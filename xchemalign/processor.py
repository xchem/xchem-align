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

from . import utils


METADATA_FILENAME = 'metadata.yaml'
CONFIG_FILENAME = 'config.yaml'
VERSION_DIR_PREFIX = 'upload_'


class Processor:
    """
    Base class that understands versioning and basic configuration
    """

    def __init__(self, config_file, logger=None):
        self.config_file = config_file

        config = utils.read_config_file(config_file)

        self.base_dir = config['base_dir']
        self.input_dirs = config['input_dirs']
        self.output_dir = config['output_dir']
        self.target_name = config['target_name']
        self.config = config
        self.version_dir = None
        self.meta_history = []
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = os.path.join(self.output_dir, VERSION_DIR_PREFIX + str(version))
            if os.path.isdir(v_dir):
                version += 1
            else:
                break
        if version == 1:
            self.logger.error('No version directory found. Please create one named upload_1')
            return None

        # the working version dir is one less than the current value
        version -= 1
        self.logger.info('Version is {}'.format(version))
        v_dir = os.path.join(self.output_dir, VERSION_DIR_PREFIX + str(version))

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(1, version):
                self.logger.info('Reading metadata for version {}'.format(v))
                dir_name = os.path.join(self.output_dir, VERSION_DIR_PREFIX + str(v))
                meta = self.read_metadata(dir_name)
                self.meta_history.append(meta)

        self.logger.info('Setting version dir to {}'.format(v_dir))
        self.version_dir = v_dir

        num_old_metas = len(self.meta_history)
        if num_old_metas:
            self.logger.info('Found {} metadata files from previous versions'.format(num_old_metas))

        return v_dir

    def read_metadata(self, version_dir):
        self.logger.info('Reading metadata for version {}'.format(version_dir))
        meta_file = os.path.join(version_dir, METADATA_FILENAME)

        meta = utils.read_config_file(meta_file)
        return meta
