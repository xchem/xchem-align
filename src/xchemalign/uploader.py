# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import re
import argparse
import tarfile
from pathlib import Path
from urllib.parse import urljoin, urlsplit
import time
import datetime
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess

from tqdm import tqdm
import yaml
from yaml.parser import ParserError
import requests
from requests.exceptions import JSONDecodeError, ConnectionError

from requests_toolbelt.multipart.encoder import (
    MultipartEncoder,
    MultipartEncoderMonitor,
)


from xchemalign import utils


# for compression
NUM_PROCESSES = max(1, os.cpu_count() - 1)

LOGIN_URL = "/accounts/login/"
UPLOAD_URL = "/api/upload_target_experiments/"
LANDING_PAGE_URL = '/viewer/react/landing/'

# used to find the url from env variables
FRAGALYSIS_URL_PREFIX = "XCHEMALIGN_FRAGALYSIS_URL_"

DEFAULT_TARBALL_TEMPLATE = "{target}_v{version}_{upload_no}_{date}.tgz"


DATA_DIR = "upload-current"


# this needs to be kept more or less up to date
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
# USER_AGENT = "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:133.0) Gecko/20100101 Firefox/133.0"

META_ALIGNER = 'meta_aligner.yaml'
CONFIG_FILE = "config.yaml"


logger = utils.Logger()


def read_yaml(yfile):
    try:
        return yaml.safe_load(yfile.read())
    except ParserError as exc:
        raise ParserError(f"{yfile.name} is not a valid YAML file") from exc


def get_validation_files_from_tarball(input_file):
    with tarfile.open(input_file, 'r') as tar:
        try:
            metafile = next(filter(re.compile(f'upload_\d+/{META_ALIGNER}').match, tar.getnames()))
        except StopIteration:
            raise FileNotFoundError(f"'{META_ALIGNER}' not found in the tarball")

        extracted_meta = tar.extractfile(metafile)
        if not extracted_meta:
            raise ValueError(f"Failed to extract '{META_ALIGNER}' from the tarball")

        try:
            configfile = next(filter(re.compile(f'upload_\d+/{CONFIG_FILE}').match, tar.getnames()))
        except StopIteration:
            raise FileNotFoundError(f"'{CONFIG_FILE}' not found in the tarball")

        extracted_config = tar.extractfile(configfile)
        if not extracted_config:
            raise ValueError(f"Failed to extract '{CONFIG_FILE}' from the tarball")

        meta = read_yaml(extracted_meta)
        config = read_yaml(extracted_config)

    return meta, config


def get_validation_files_from_fs(upload_path):
    meta_path = upload_path.joinpath(META_ALIGNER)
    if not meta_path.exists():
        raise FileNotFoundError(f"'{META_ALIGNER}' not found")

    config_path = upload_path.joinpath(CONFIG_FILE)
    if not config_path.exists():
        raise FileNotFoundError(f"'{CONFIG_FILE}' not found")

    with open(str(meta_path), "r") as meta_file, open(str(config_path), "r") as config_file:
        meta = read_yaml(meta_file)
        config = read_yaml(config_file)

    return meta, config


def get_validation_data(meta, config):
    try:
        data_version = meta["data_format_version"]
    except KeyError as exc:
        raise KeyError(f"Key 'data_format_version' missing from {META_ALIGNER}") from exc

    try:
        target_name = config["target_name"]
    except KeyError as exc:
        raise KeyError(f"Key 'target_name' missing from {CONFIG_FILE}") from exc

    return {
        'data_version': data_version,
        'target_name': target_name,
    }


def count_files(directory):
    """Count total files in directory."""
    # used to display a progress bar
    return sum(len(files) for _, _, files in os.walk(directory))


def iter_files(directory):
    """Generator that yields file paths."""
    for root, _, files in os.walk(directory):
        for file in files:
            yield os.path.join(root, file)


def detect_compression():
    """Detect the best available compression tool."""
    if shutil.which("pigz"):
        return "pigz"
    elif shutil.which("gzip"):
        return "gzip"
    else:
        return "python"


def compress_directory(upload_path, tarball_path):
    """Compress the tarball using the best available method."""
    compression_tool = detect_compression()

    parent_dir = upload_path.absolute().parent

    logger.info(f"Using {compression_tool} for compression...")

    if compression_tool in ["pigz", "gzip"]:
        # using subprocess to show a progress bar
        with tqdm(
            total=os.path.getsize(upload_path),
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=f"Compressing {upload_path}",
            dynamic_ncols=True,
        ) as pbar, open(tarball_path, "wb") as output_file:
            # due to the way gzip and pigz work, specifically, not
            # allowing to set the output file but use pipes instead, I
            # have to create 2 processes, one creates the zipped
            # stream and sends it to stdout, the other one catches it
            # and creates the file. really annoying but worth it
            # because so much faster.
            tar_process = subprocess.Popen(
                ["tar", "-C", parent_dir, "-cf", "-", upload_path.name],
                stdout=subprocess.PIPE,
            )
            compress_process = subprocess.Popen(
                [compression_tool, "-c"],
                stdin=tar_process.stdout,
                stdout=output_file,
            )
            while compress_process.poll() is None:
                pbar.update(1024 * NUM_PROCESSES)  # approximate progress step

            tar_process.stdout.close()
            tar_process.wait()
            compress_process.wait()

    else:
        # no system tools available, going for the fallback option
        total_files = count_files(upload_path)
        with tarfile.open(tarball_path, "w:gz") as tar:
            with ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:
                file_list = iter_files(upload_path)
                with tqdm(
                    total=total_files,
                    unit="file",
                    desc="Creating tarball",
                    dynamic_ncols=True,
                ) as pbar:
                    # executor.submit expects a function, hence the
                    # lambda event though it doesn't do much
                    futures = {executor.submit(lambda x: x, file): file for file in file_list}
                    for future in as_completed(futures):
                        file_path = future.result()
                        tar.add(
                            file_path,
                            arcname=os.path.relpath(
                                file_path,
                                start=os.path.dirname(upload_path),
                            ),
                        )
                        pbar.update(1)

    logger.info(f"Tarball saved at: {tarball_path}")


def get_upload_file(use_default=False, use_custom=None):
    # figure out data to upload. 3 options, create on the fly, use
    # existing or use custom
    upload_path = None
    if use_custom:
        tarball_path = use_custom
        meta, config = get_validation_files_from_tarball(tarball_path)
        validation_data = get_validation_data(meta, config)
    else:
        # no custom filename means I need to figure it out and that
        # means looking at the files in the fs
        root_path = Path(DATA_DIR)
        try:
            # get the last 'upload_x' folder
            *_, upload_path = root_path.glob("upload_*")
        except ValueError:
            raise FileNotFoundError("Upload directory not found")

        meta, config = get_validation_files_from_fs(upload_path)
        validation_data = get_validation_data(meta, config)
        upload_dir = upload_path.parts[-1]

        try:
            target_name = config["target_name"]
        except KeyError as exc:
            raise KeyError(f"'target_name' attribute missing from {CONFIG_FILE}") from exc

        tarball_path = root_path.joinpath(
            DEFAULT_TARBALL_TEMPLATE.format(
                target=target_name,
                version=validation_data["data_version"],
                upload_no=upload_dir,
                date=datetime.datetime.today().strftime('%Y-%m-%d'),
            )
        )

        if use_default:
            if not tarball_path.exists():
                raise FileNotFoundError(f"Tarball {str(tarball_path)} not found")

        # at this point, the decision is that a new tarball needs to
        # be created. but don't do it just yet, first run a validation

    return upload_path, tarball_path, validation_data


def upload(url, proposal, auth_token=None, use_default=False, use_custom=None):
    # figure out necessary urls
    splits = urlsplit(url)
    base_url = f'{splits.scheme}://{splits.netloc}'
    upload_url = urljoin(base_url, UPLOAD_URL)
    landing_page_url = urljoin(url, LANDING_PAGE_URL)

    try:
        upload_path, input_file, validation_data = get_upload_file(
            use_default=use_default,
            use_custom=use_custom,
        )
    except (FileNotFoundError, ValueError, KeyError) as exc:
        logger.error(exc.args)
        logger.error(exc.args[0])
        return

    filename = Path(input_file).name

    with requests.Session() as session:
        session.headers.update(
            {
                "User-Agent": USER_AGENT,
                "Referer": landing_page_url,
                "Referrer-policy": "same-origin",
            }
        )

        # session = check_deployment(session, auth_token)
        session.get(landing_page_url)  # sets csrftoken

        csrftoken = session.cookies.get('csrftoken', None)
        if csrftoken:
            session.headers.update(
                {
                    "X-CSRFToken": csrftoken,
                    "User-Agent": USER_AGENT,
                }
            )

        if auth_token:
            session.cookies.update(
                {
                    "sessionid": auth_token,
                }
            )

        validation_data["target_access_string"] = proposal

        logger.info("Checking data version...")
        validation_result = session.post(
            upload_url,
            data=validation_data,
        )
        if validation_result.url.find("keycloak") > 0:
            logger.error("You are not logged in to Fragalysis")
            return

        result_json = validation_result.json()
        if validation_result.ok:
            # data validation errors
            if "success" in result_json.keys():
                if not result_json["success"]:
                    logger.error(result_json["message"])
                    return
            elif "detail" in result_json.keys():
                logger.error(result_json["detail"])
                return
            else:
                if result_json["message"]:
                    logger.warn(result_json["message"])
        else:
            # django validation errors. The only think that can be
            # here is propsal (which on the server is called
            # target_access_string
            tas_error = result_json.pop("target_access_string")
            if tas_error:
                logger.error(f"Proposal number error: {tas_error[0]}")

            # but just in case, if there's anything else, print them
            for field, errors in result_json.items():
                # because comes as a list
                for e in errors:
                    logger.error(f"{field}: {e}")
            return

        # validation passed, attempt file upload. create if necessary
        if upload_path and upload_path.exists():
            compress_directory(upload_path, input_file)

        encoder = MultipartEncoder(
            fields={
                "target_access_string": proposal,
                "file": (
                    str(input_file),
                    open(input_file, "rb"),
                    "application/octet-stream",
                ),
            }
        )

        with open(input_file, "rb") as f:
            file_size = int(f.seek(0, 2))
            f.seek(0)

            def callback(monitor):
                pbar.update(monitor.bytes_read - pbar.n)

            monitor = MultipartEncoderMonitor(encoder, callback)

            session.headers.update({"Content-Type": monitor.content_type})

            with tqdm(total=file_size, unit="B", unit_scale=True, desc='Uploading') as pbar:
                try:
                    response = session.post(
                        upload_url,
                        data=monitor,
                    )
                except ConnectionError as exc:
                    logger.info('Connection closed by server')
                    return

        if response.ok:
            # not quite there yet
            try:
                jresp = response.json()
                # and once more
                try:
                    task_status_url = jresp["task_status_url"]
                except KeyError:
                    logger.error('Server response does not contain task url')
                    return
            except JSONDecodeError:
                # I have occasionally observed it, but don't know how
                # to force it for debugging
                logger.error('Server response does not contain json payload')
                return

            # should get back "task_status_url" in response. keep polling for updates
            task_status = {}
            count = 0
            status_url = urljoin(base_url, task_status_url)
            finished = False
            errors = []

            logger.info(f'{filename} uploaded, processing:')
            while True:
                # ping the task url endpoint for a while and print the messages
                status = session.get(status_url)
                # print('ping', status.json())

                count = count + 1
                try:
                    statjson = status.json()
                except JSONDecodeError:
                    continue

                if task_status != statjson:
                    task_status = statjson
                    msg = task_status.get("messages", None)
                    if isinstance(msg, list):
                        for m in msg:
                            logger.log(m, level=-1)
                            if m.startswith('ERROR'):
                                errors.append(m)
                    else:
                        logger.info(f'{task_status.get("status", "")}: {msg}')
                    if task_status in ("SUCCESS", "FAILED", "CANCELED", "FATAL"):
                        finished = True
                        break

                if statjson.get("ready", False) is True:
                    # if task complete, quit
                    finished = True
                    break
                if statjson.get("status", "") in ("SUCCESS", "FAILED", "CANCELED", "FATAL"):
                    # if task complete, quit
                    finished = True
                    break
                if count > 100:
                    # if no changes for {randomly selected number of pings}, quit
                    break

                time.sleep(1)

            if finished:
                if errors:
                    logger.log(f'The following errors were encountered when processing {filename}:', level=-1)
                    for e in errors:
                        logger.log(e, level=-1)

        else:
            # upload response else, failure
            logger.error("Upload failed: error {}, {}".format(response.status_code, response.reason))


def main():
    # grab the predefined urls from the environment
    fragalysis_urls = {
        key[len(FRAGALYSIS_URL_PREFIX) :].lower(): value
        for key, value in os.environ.items()
        if key.startswith(FRAGALYSIS_URL_PREFIX)
    }

    # compile url help text, listing, if found, the keyword urls
    url_help = "Upload url."
    if fragalysis_urls:
        url_help = (
            url_help
            + " Use any of the following choice keywords or a custom url:"
            + os.linesep
            + os.linesep.join([f"{key} ({value})" for key, value in fragalysis_urls.items()])
        )

    parser = argparse.ArgumentParser(
        description="Upload files to fragalysis",
        formatter_class=argparse.RawTextHelpFormatter,  # preserves newlines
    )
    parser.add_argument("-u", "--url", metavar="url", required=True, help=url_help)
    parser.add_argument("-p", "--proposal", metavar="proposal number", required=True)
    parser.add_argument(
        "-t",
        "--token",
        default=None,
        required=False,
        metavar="authentication token",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-d",
        "--use-default",
        action="store_true",
        help="Upload previously created tarball instead of creating a new one",
    )
    group.add_argument(
        "-c",
        "--use-custom",
        metavar="tarball",
        required=False,
        help="Upload custom tarball",
    )

    args = parser.parse_args()

    if args.url in fragalysis_urls.keys():
        args.url = fragalysis_urls[args.url]

    upload(
        args.url,
        args.proposal,
        auth_token=args.token,
        use_default=args.use_default,
        use_custom=args.use_custom,
    )


if __name__ == "__main__":
    main()
