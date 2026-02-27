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
import hashlib
from pathlib import Path
from urllib.parse import urljoin, urlsplit
import time
import datetime
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
from tempfile import TemporaryDirectory
from typing import Any
from abc import ABC, abstractmethod

from tqdm import tqdm
import boto3
from boto3.s3.transfer import TransferConfig
from botocore.config import Config
from botocore.exceptions import ParamValidationError
import yaml
from yaml.parser import ParserError
import requests
from requests.exceptions import JSONDecodeError, ConnectionError

from requests_toolbelt.multipart.encoder import (
    MultipartEncoder,
    MultipartEncoderMonitor,
)


from xchemalign import utils
from .copier import handle_inputs

# from xchemalign.utils import Constants


# for compression
NUM_PROCESSES = max(1, os.cpu_count() - 1)
AWS_BUCKET_NAME = os.environ.get("AWS_BUCKET_NAME", '')


# this needs to be kept more or less up to date
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"


logger = utils.LOG


def compile_stack_urls(url_prefix: str) -> dict[str, str]:
    """Extract available stack urls from env variables"""
    return {key[len(url_prefix) :].lower(): value for key, value in os.environ.items() if key.startswith(url_prefix)}


def calculate_sha256(filepath) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read the file in chunks of 4096 bytes
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def read_yaml(yfile):
    try:
        return yaml.safe_load(yfile.read())
    except ParserError as exc:
        raise ParserError(f"{yfile.name} is not a valid YAML file") from exc


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


def upload_to_s3(tarball_path, object_key=""):
    logger.info("Uploading to s3")

    filesize = os.path.getsize(tarball_path)

    # workaround to a recent boto3 issue: https://github.com/boto/boto3/issues/4398
    config = Config(
        request_checksum_calculation="WHEN_REQUIRED",
        response_checksum_validation="WHEN_REQUIRED",
    )
    s3 = boto3.client("s3", config=config)

    progress = tqdm(
        total=filesize, unit="B", unit_scale=True, unit_divisor=1024, desc=f"Uploading {tarball_path.name}"
    )

    def progress_callback(bytes_amount):
        progress.update(bytes_amount)

    with open(str(tarball_path), "rb") as f:
        try:
            s3.upload_fileobj(f, AWS_BUCKET_NAME, object_key, Callback=progress_callback)
        except ParamValidationError as exc:
            logger.error(exc.args[0])
            logger.error("Do you have the AWS_BUCKET_NAME environment variable set?")

    progress.close()


class AuthenticationError(Exception):
    pass


class ValidationError(Exception):
    pass


class FileTransferError(Exception):
    pass


class RetryLimitReachedError(Exception):
    pass


class StatusError(Exception):
    pass


class DataSource(ABC):
    META_ALIGNER_FILE = 'meta_aligner.yaml'
    CONFIG_FILE = "config.yaml"

    DEFAULT_TARBALL_TEMPLATE = "{target}_v{version}_{upload_no}_{date}.tgz"
    DEFAULT_INPUTS_TEMPLATE = "{target}_v{version}_{upload_no}_{date}_inputs.tgz"

    def __init__(self) -> None:
        self._config = None
        self._meta = None
        self._data_version = None
        self._upload_version = None
        self._target_name = None
        self._tarball_path = None
        self._data_source = None

    @property
    @abstractmethod
    def data_source(self) -> Path:
        pass

    @property
    def config(self) -> dict[str, Any]:
        if self._config is None:
            self._config = self._load_yaml_data(self.CONFIG_FILE)
        return self._config

    @property
    def meta(self) -> dict[str, Any]:
        if self._meta is None:
            self._meta = self._load_yaml_data(self.META_ALIGNER_FILE)
        return self._meta

    @property
    def data_version(self) -> str:
        if self._data_version is None:
            try:
                self._data_version = self.meta["data_format_version"]
            except KeyError as exc:
                raise KeyError(
                    f"Key 'data_format_version' missing from {self.META_ALIGNER_FILE}",
                ) from exc
        return self._data_version

    @property
    def upload_version(self) -> str:
        if self._upload_version is None:
            try:
                self._upload_version = self.meta["version_number"]
            except KeyError as exc:
                raise KeyError(
                    f"Key 'version_number' missing from {self.META_ALIGNER_FILE}",
                ) from exc
        return self._upload_version

    @property
    def target_name(self) -> str:
        if self._target_name is None:
            try:
                self._target_name = self.config["target_name"]
            except KeyError as exc:
                raise KeyError(
                    f"Key 'target_name' missing from {self.CONFIG_FILE}",
                ) from exc
        return self._target_name

    @property
    def tarball_path(self) -> Path:
        if self._tarball_path is None:
            self._tarball_path = self._get_tarball_path()
        return self._tarball_path

    @abstractmethod
    def _load_yaml_data(self, filename):
        raise NotImplementedError

    @abstractmethod
    def _get_tarball_path(self) -> Path:
        raise NotImplementedError

    def get_validation_payload(self) -> dict[str, str]:
        return {
            'data_version': self.data_version,
            'target_name': self.target_name,
            'upload_version': self.upload_version,
        }

    def checksum(self):
        return calculate_sha256(self.tarball_path)


class XCADataUpload:
    """Handle validation, upload, and status check of data tarball."""

    LOGIN_URL = "/accounts/login/"
    VALIDATE_URL = "/api/validate_target_experiments/"
    UPLOAD_URL = "/api/upload_target_experiments/"
    LANDING_PAGE_URL = '/viewer/react/landing/'

    def __init__(
        self,
        url: str,
        data_source: DataSource,
        proposal: str,
        auth_token: str | None = None,
        retries: int = 3,
    ) -> None:
        self._data_source = data_source
        self._proposal = proposal
        self._auth_token = auth_token
        self._retries = retries

        splits = urlsplit(url)
        self._base_url = f'{splits.scheme}://{splits.netloc}'
        self._validate_url = urljoin(self._base_url, self.VALIDATE_URL)
        self._upload_url = urljoin(self._base_url, self.UPLOAD_URL)
        self._landing_page_url = urljoin(url, self.LANDING_PAGE_URL)

        self._session = requests.Session()

        # the way fragalysis task works, all messages are retained and
        # retrieved by task status checker. potential for
        # improvement. for now, keep track of how many messages have
        # been printed and don't show the old ones
        self._progress_message_count = 0
        self._progress_errors = []

        self._futile_ping_count = 100

    def _init_session(self):
        """Start session and populate necessary params"""
        self._session.headers.update(
            {
                "User-Agent": USER_AGENT,
                "Referer": self._landing_page_url,
                "Referrer-policy": "same-origin",
            }
        )

        # session = check_deployment(session, auth_token)
        self._session.get(self._landing_page_url)  # sets csrftoken

        csrftoken = self._session.cookies.get('csrftoken', None)
        if csrftoken:
            self._session.headers.update(
                {
                    "X-CSRFToken": csrftoken,
                    "User-Agent": USER_AGENT,
                }
            )

        if self._auth_token:
            self._session.cookies.update(
                {
                    "sessionid": self._auth_token,
                }
            )

    def _validate(self):
        """Validate upload data.

        Checks data version, proposal membership, authentication (if given).
        """
        logger.info('Validating upload')
        validation_data = self._data_source.get_validation_payload()
        logger.info(f'Local data version: {validation_data["data_version"]}')
        logger.info(f'Local upload version: {validation_data["upload_version"]}')

        validation_data["target_access_string"] = self._proposal

        validation_result = self._session.post(
            self._validate_url,
            data=validation_data,
        )
        if validation_result.url.find("keycloak") > 0:
            raise AuthenticationError("You are not logged in to Fragalysis")

        result_json = validation_result.json()

        if validation_result.ok:
            # data validation errors
            if "success" in result_json.keys():
                if not result_json["success"]:
                    errors = False
                    for msg in result_json["message"]:
                        errors = True
                        logger.error(msg)
                    if errors:
                        # need to raise error.. but message superfluous?
                        raise ValidationError('Data validation errors, quitting')

            elif "detail" in result_json.keys():
                raise ValidationError(result_json["detail"])

            else:
                if result_json["message"]:
                    for msg in result_json["message"]:
                        logger.warn(msg)
        else:
            # django validation errors. The only think that can be
            # here is propsal (which on the server is called
            # target_access_string
            tas_error = result_json.pop("target_access_string")
            if tas_error:
                raise ValidationError(f"Proposal number error: {tas_error[0]}")

            # but just in case, if there's anything else, print them
            for field, errors in result_json.items():
                # because comes as a list
                if errors:
                    for e in errors:
                        logger.error(f"{field}: {e}")

                    raise ValidationError('Server validation errors')

    def _upload(self) -> str:
        """File upload attempt"""
        checksum = self._data_source.checksum()
        encoder = MultipartEncoder(
            fields={
                "target_access_string": self._proposal,
                "sha256checksum": checksum,
                "file": (
                    str(self._data_source.tarball_path),
                    open(self._data_source.tarball_path, "rb"),
                    "application/octet-stream",
                ),
            }
        )
        with open(self._data_source.tarball_path, "rb") as f:
            file_size = int(f.seek(0, 2))
            f.seek(0)

            def callback(monitor):
                pbar.update(monitor.bytes_read - pbar.n)

            monitor = MultipartEncoderMonitor(encoder, callback)

            self._session.headers.update({"Content-Type": monitor.content_type})

            with tqdm(total=file_size, unit="B", unit_scale=True, desc='Uploading') as pbar:
                try:
                    response = self._session.post(
                        self._upload_url,
                        data=monitor,
                    )
                except ConnectionError as exc:
                    raise ConnectionError('Connection closed by server') from exc

        if response.ok:
            try:
                response_json = response.json()
                # and once more
            except JSONDecodeError as exc:
                # I have occasionally observed it, but don't know how
                # to force it for debugging
                logger.error('Server response does not contain json payload')
                raise JSONDecodeError from exc

            try:
                response_status_url = response_json["task_status_url"]
                task_status_url = urljoin(self._base_url, response_status_url)
                self._task_status_url = task_status_url
                logger.info(f"task_status_url={task_status_url}")
            except KeyError:
                logger.error('Server response does not contain task url')
                raise Exception

            return task_status_url

        else:
            try:
                response_json = response.json()
            except JSONDecodeError as exc:
                # try something else?
                msg = 'Response body does not contain JSON payload'
                logger.error(msg)
                logger.error(response.text)
                raise JSONDecodeError(msg) from exc
            except Exception as exc:
                # other errors, no idea what they might be
                msg = "Upload failed: error {}, {}".format(response.status_code, response.reason)
                logger.error(msg)
                logger.error(response.text)
                raise Exception(msg) from exc

            # no network errors, response resoved correctly, but there might
            # see if the server raises any issues
            if fname_errors := response_json.get('filename', ''):
                # comes as a list, have to go through everything
                for e in fname_errors:
                    if e.find('checksum'):
                        raise FileTransferError
                    else:
                        raise Exception(f'Upload file error: {e}')
            else:
                raise Exception

    def _check_status(self, url: str):
        """Check task status in Fragalysis.

        Return new messaages about processing progress.
        """
        status = self._session.get(url)
        finished = False

        status_json = status.json()

        if status_json.get('status', None) in ("SUCCESS", "FAILED", "CANCELED", "FATAL"):
            finished = True

        if status_json.get("ready", False) is True:
            finished = True

        if error := status_json.get("error", None):
            raise StatusError(error)

        msg = status_json.get("messages", None)

        if isinstance(msg, list):
            messages = msg
        else:
            messages = [msg]

        return finished, messages

    def _close(self):
        self._session.close()

    # def _upload_to_s3(self):
    #     copy_inputs(
    #         base_dir, inputs, ref_datasets, inputs_path  # pylint: disable=undefined-variable
    #     )  # pylint: disable=used-before-assignment
    #     object_key = f"{validation_data['target_name']}/{str(inputs_path)}"
    #     upload_to_s3(inputs_path, object_key=object_key)

    def upload(self):
        logger.info(f'Uploading {self._data_source.data_source}')

        try:
            self._init_session()
            self._validate()
            logger.info('Validation successful, starting file upload')

            for i in range(self._retries):
                try:
                    task_status_url = self._upload()
                    break
                except FileTransferError:
                    logger.info(
                        'Uploaded file checksum does not match the calculated checksum, '
                        + f'file was likely corrupt during the transfer. Retrying {i + 1}',
                    )
            else:
                raise RetryLimitReachedError(f'Upload retry limit ({self._retries}) reached, quitting')

            logger.info('File uploaded, checking progress...')
            count = 0
            missing_proj_msg_already_shown = False
            while True:
                time.sleep(1)
                try:
                    finished, messages = self._check_status(url=task_status_url)
                except JSONDecodeError:
                    # this is fine, task is in such an early state
                    # that the message body has not been initiated yet
                    continue
                except StatusError as exc:
                    # this is an error raised by server, see what it says
                    if exc.args[0] == f'Proposal {self._proposal} not found':
                        # common in staging which gets wiped often and
                        # users upload to projects that don't exist
                        # yet, but is being created by this very
                        # upload. safe to ignore right now, remove
                        # later when project creation method is
                        # changed
                        if not missing_proj_msg_already_shown:
                            # printing this once is enough
                            missing_proj_msg_already_shown = True
                            logger.info(f'Proposal {self._proposal} not found!')
                            logger.info(
                                'This indicates that the proposal does not yet exist and '
                                + 'will now be created. It is not possible to display live '
                                + 'updates of the target loading process; all collected '
                                + 'messages will be presented once the process has completed. '
                                + 'Depending on the size of the tarball, this may take some time.'
                            )

                        # reset counter. Not ideal if gets stuck but nothing to do here
                        count = 0
                        continue
                    else:
                        raise StatusError from exc

                if self._progress_message_count == 0 or self._progress_message_count != len(messages):
                    # no old messages, everything or at least
                    # some incoming messages are new. reset counter,
                    # don't restart loop, go to printing
                    count = 0
                else:
                    # no new messages, nothing to print, start new
                    # cycle and continue checking
                    count = count + 1
                    continue

                if count > self._futile_ping_count:
                    # if no changes for {randomly selected number of pings}, quit
                    logger.info(f'No changes in {self._futile_ping_count} pings, quitting')
                    break

                # db message buffer contains all messages that have
                # been written from the beginning of the process. only
                # print the ones that have not been printed yet. if
                # error, store for showing later
                for m in messages[self._progress_message_count :]:
                    logger.log(m, level=-1)
                    if m.startswith('ERROR'):
                        self._progress_errors.append(m)

                self._progress_message_count = len(messages)

                if finished:
                    break

        finally:
            # at the end, print all the errors once more
            if self._progress_errors:
                logger.log(
                    'The following errors were encountered when processing ' + 'f{self._data_source.tarball_path}:',
                    level=-1,
                )
                for e in self._progress_errors:
                    logger.log(e, level=-1)
            self._close()


class TarballSource(DataSource):
    def __init__(self, tarball_path: str) -> None:
        super().__init__()
        self._tarball_path: Path = Path(tarball_path)

    @property
    def data_source(self) -> Path:
        return self._tarball_path

    def _get_tarball_path(self) -> Path:
        """Override: validate exist"""
        if not self._tarball_path.exists():
            raise FileNotFoundError(f"Tarball {str(self._tarball_path)} not found")
        return self._tarball_path

    def _load_yaml_data(self, filename):
        """Override: extract yaml files from the tarball"""
        logger.info(f'Extracting {filename} from tarball')
        with tarfile.open(self.tarball_path, 'r') as tar:
            try:
                yaml_file = next(
                    filter(
                        re.compile(f'upload_\d+/{filename}').match,
                        tar.getnames(),
                    ),
                )
            except StopIteration:
                raise FileNotFoundError(f"'{filename}' not found in the tarball")

            extracted = tar.extractfile(yaml_file)
            if not extracted:
                raise ValueError(f"Failed to extract '{self.META_ALIGNER_FILE}' from the tarball")

            parsed = read_yaml(extracted)
        return parsed


class FilesystemSource(DataSource):
    DATA_DIR = "upload-current"

    def __init__(self, use_default: bool = False) -> None:
        super().__init__()
        self.use_default: bool = use_default
        self._upload_path: Path | None = None
        self._root_path: Path | None = None

    @property
    def data_source(self) -> Path:
        if self.use_default:
            return self._tarball_path
        else:
            return self.upload_path

    @property
    def tarball_path(self) -> Path:
        """Override: compress directory to tarball if not yet done"""
        if self._tarball_path is None:
            self._tarball_path = self._get_tarball_path()

        if not self._tarball_path.exists():
            compress_directory(self.upload_path, self._tarball_path)
        return self._tarball_path

    @property
    def upload_path(self) -> Path:
        """Data directory, uncompressed"""
        if self._upload_path is None:
            # uploader can be run from 3 places
            root_path = Path(self.DATA_DIR).absolute()
            if not root_path.exists():
                if root_path.parent.name == self.DATA_DIR:
                    root_path = root_path.parent
                elif root_path.parent.parent.name == self.DATA_DIR:
                    root_path = root_path.parent.parent

            if not root_path.exists():
                raise FileNotFoundError("xchemalign.uploader run from unexpected location")

            self._root_path = root_path

            try:
                # get the last 'upload_x' folder
                *_, upload_path = root_path.glob("upload_*")
            except ValueError:
                raise FileNotFoundError("Upload directory not found")

            self._upload_path = upload_path

        return self._upload_path

    def upload_inputs(self):
        # NB! method not tested
        validation_data = self.get_validation_payload()

        # TODO: a possibility to check with the fragalysis whether the
        # target upload was successful. This will require target name and
        # upload version from validation_data (returned by
        # get_upload_file()). what I dont't have here atm, is the stack
        # url, if going to check this, this parameter must be made
        # available from cli

        date = datetime.datetime.today().strftime('%Y-%m-%d')
        upload_dir = self.upload_path.parts[-1]

        inputs_path = self._root_path.joinpath(
            self.DEFAULT_INPUTS_TEMPLATE.format(
                target=self.target_name,
                version=self.data_version,
                upload_no=upload_dir,
                date=date,
            )
        )

        object_key = f"{validation_data['target_name']}/{str(inputs_path)}"

        if not inputs_path.is_file():
            base_dir = self._get_copy_inputs(self.config)
            self._copy_inputs(base_dir, self.config, inputs_path)

        upload_to_s3(inputs_path, object_key=object_key)

    def _get_tarball_path(self) -> Path:
        """Override: compile tarball path"""
        date = datetime.datetime.today().strftime('%Y-%m-%d')
        upload_dir = self.upload_path.parts[-1]

        tarball_path = self._root_path.joinpath(
            self.DEFAULT_TARBALL_TEMPLATE.format(
                target=self.target_name,
                version=self.data_version,
                upload_no=upload_dir,
                date=date,
            )
        )

        if self.use_default:
            if not self.tarball_path.exists():
                raise FileNotFoundError(f"Tarball {str(self.tarball_path)} not found")

        return tarball_path

    def _load_yaml_data(self, filename):
        """Override: find yaml files from the filesystem"""
        path = self.upload_path.joinpath(filename)
        if not path.exists():
            raise FileNotFoundError(f"'{filename}' not found")

        with open(str(path), "r") as input_file:
            contents = read_yaml(input_file)

        return contents

    @classmethod
    def _get_copy_inputs(cls, config):
        try:
            base_dir = config["base_dir"]
        except KeyError as exc:
            raise KeyError(f"Key 'base_dir' missing from {cls.CONFIG_FILE}") from exc

        return base_dir

    @classmethod
    def _copy_inputs(cls, base_dir, config, inputs_path):
        with TemporaryDirectory() as tempdir:
            compressible_tempdir = Path(tempdir).joinpath("inputs")
            try:
                handle_inputs(base_dir, config, str(compressible_tempdir), logger)
            except Exception as exc:
                logger.error(f"Error copying inputs: {exc.args}")
            compress_directory(compressible_tempdir, inputs_path)


class XCAUploadManager:
    # prefix used in fragalysis urls in env variables
    _FRAGALYSIS_URL_PREFIX = "XCHEMALIGN_FRAGALYSIS_URL_"

    # dict (short name: url) of available fragalysis stack URLs
    # defined in user's environment (staging, production, local dev,
    # etc)
    _FRAGALYSIS_URLS = compile_stack_urls(_FRAGALYSIS_URL_PREFIX)

    def __init__(
        self,
        *,
        url: str,
        proposal: str,
        use_default: bool = False,
        use_custom: str | None = None,
        no_copy: bool = False,
        retries: int = 3,
        auth_token: str | None = None,
    ) -> None:
        self._url = self.get_stack_url(url)
        self._proposal = proposal
        self._auth_token = auth_token
        self._no_copy = no_copy
        self._no_copy = False  # no copying atm
        self._retries = retries

        if use_custom:
            # user gave the path to the tarball they want to upload
            data_source = TarballSource(tarball_path=use_custom)
        else:
            # tarball not given, that means looking for the files in
            # the fs. 2 possibilities:
            # 1. compress the input directory and upload the tarball
            # 2. this is a retry and the tarball was already created in 1.
            data_source = FilesystemSource(use_default=use_default)

        self.data_source: DataSource = data_source

    def upload_target(self):
        uploader = XCADataUpload(
            url=self._url,
            data_source=self.data_source,
            proposal=self._proposal,
            auth_token=self._auth_token,
            retries=self._retries,
        )

        try:
            uploader.upload()
        except KeyboardInterrupt:
            logger.info(
                'Quitting. If the tarball was already uploaded, this will not ' + 'stop the processing in Fragalys.'
            )
        except (
            AuthenticationError,
            ValidationError,
            ConnectionError,
            RetryLimitReachedError,
            FileNotFoundError,
            StatusError,
        ) as exc:
            logger.error(exc)

    @staticmethod
    def upload_inputs():
        data_source = FilesystemSource()

        try:
            data_source.upload_inputs()
        except KeyboardInterrupt:
            logger.info(
                'Quitting. If the tarball was already uploaded, this will not ' + 'stop the processing in Fragalys.'
            )
        except (FileNotFoundError,) as exc:
            logger.error(exc)

    @classmethod
    def get_stack_url(cls, url_name) -> str:
        """Return a valid url based on user-supplied short name

        If url_name not in precompiled urls, assume valid URL and
        continue.
        Possble TODO: add validation
        """
        return cls._FRAGALYSIS_URLS.get(url_name, url_name)

    @classmethod
    def get_url_help_text(cls) -> str:
        """Compile url help text, listing, if found, the keyword urls"""
        url_help = "Upload url."
        if cls._FRAGALYSIS_URLS:
            url_help = (
                url_help
                + " Use any of the following choice keywords or a custom url:"
                + os.linesep
                + os.linesep.join(
                    [f"{key} ({value})" for key, value in cls._FRAGALYSIS_URLS.items()],
                )
            )

        return url_help


def main():
    url_help = XCAUploadManager.get_url_help_text()

    parser = argparse.ArgumentParser(
        description="Upload files to fragalysis",
        formatter_class=argparse.RawTextHelpFormatter,  # preserves newlines
    )

    # two mutually exclusive groups of input args:
    # 1) upload data to fragalysis (with it's own mutually exclusive args)
    # 2) upload source data
    # since argparse doesn't allow nested groups, have to validate
    # manually

    parser.add_argument("-u", "--url", metavar="url", help=url_help)
    parser.add_argument("-p", "--proposal", metavar="proposal number")
    parser.add_argument(
        "-t",
        "--token",
        default=None,
        required=False,
        metavar="authentication token",
    )

    parser.add_argument(
        "-d",
        "--use-default",
        action="store_true",
        help="Upload previously created tarball instead of creating a new one",
    )
    parser.add_argument(
        "-c",
        "--use-custom",
        metavar="tarball",
        required=False,
        help="Upload custom tarball",
    )
    parser.add_argument(
        "--no-copy",
        required=False,
        action="store_true",
        help="Do not copy input data to destination",
    )
    parser.add_argument(
        "-s",
        "--upload-inputs",
        action="store_true",
        help="Skip the target upload and only upload XCA input data",
    )
    parser.add_argument(
        "--max-retries",
        default=3,
        type=int,
        required=False,
        metavar="N",
        help="Retry broken download N times",
    )
    args = parser.parse_args()

    if args.upload_inputs:
        logger.info(
            "'upload_inputs' argument given, ignoring any others " + "and proceeding to upload XCA input data",
        )
        XCAUploadManager.upload_inputs()
    else:
        uploader = XCAUploadManager(
            url=args.url,
            proposal=args.proposal,
            auth_token=args.token,
            use_default=args.use_default,
            use_custom=args.use_custom,
            no_copy=args.no_copy,
            retries=args.max_retries,
        )
        uploader.upload_target()


if __name__ == "__main__":
    main()
