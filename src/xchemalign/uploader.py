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

import re
import argparse
import tarfile
from pathlib import Path
from urllib.parse import urljoin, urlsplit
import time

from tqdm import tqdm
import yaml
import requests
from requests.exceptions import JSONDecodeError, ConnectionError

from requests_toolbelt.multipart.encoder import (
    MultipartEncoder,
    MultipartEncoderMonitor,
)


from xchemalign import utils


LOGIN_URL = "/accounts/login/"
UPLOAD_URL = "/api/upload_target_experiments/"
LANDING_PAGE_URL = '/viewer/react/landing/'


# this needs to be kept more or less up to date
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
# USER_AGENT = "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:133.0) Gecko/20100101 Firefox/133.0"

META_ALIGNER = 'meta_aligner.yaml'
CONFIG_FILE = "config.yaml"


logger = utils.Logger()


def get_validation_data(input_file):
    with tarfile.open(input_file, 'r') as tar:
        # Check if the specified YAML file exists in the tarball
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

        # Read the file content and parse it as YAML
        meta = yaml.safe_load(extracted_meta.read())
        config = yaml.safe_load(extracted_config.read())

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


def upload(input_file, url, proposal, auth_token=None):
    splits = urlsplit(url)
    base_url = f'{splits.scheme}://{splits.netloc}'
    upload_url = urljoin(base_url, UPLOAD_URL)

    filename = Path(input_file).name

    with requests.Session() as session:
        session.headers.update(
            {
                "User-Agent": USER_AGENT,
                "Referer": urljoin(url, LANDING_PAGE_URL),
                "Referrer-policy": "same-origin",
            }
        )

        # session = check_deployment(session, auth_token)
        session.get(url)  # sets csrftoken

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

        try:
            validation_data = get_validation_data(input_file)
        except (FileNotFoundError, ValueError, KeyError) as exc:
            logger.error(exc.args)
            logger.error(exc.args[0])
            return

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
        if not result_json["success"]:
            logger.error(result_json["message"])
            return
        else:
            if result_json["message"]:
                logger.warn(result_json["message"])

        # validation passed, attempt file upload
        encoder = MultipartEncoder(
            fields={
                "target_access_string": proposal,
                "file": (
                    input_file,
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
    parser = argparse.ArgumentParser(description="Upload files to fragalysis")
    parser.add_argument("-i", "--input", required=True, help="Input")
    parser.add_argument("-u", "--url", required=True, help="Url to upload")
    parser.add_argument("-p", "--proposal", required=True, help="Project/proposal")
    parser.add_argument("-t", "--token", default=None, required=False, help="Authentication token")

    args = parser.parse_args()

    upload(args.input, args.url, args.proposal, auth_token=args.token)


if __name__ == "__main__":
    main()
