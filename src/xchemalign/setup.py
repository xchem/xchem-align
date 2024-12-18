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

import argparse
import math
import os
from pathlib import Path


from xchemalign import utils


class Setup:
    def __init__(self, working_dir: Path):
        if not working_dir:
            self.working_dir = Path.cwd()
        else:
            self.working_dir = Path(working_dir)
        self.version_dir = Path('upload-v' + str(math.floor(utils.DATA_FORMAT_VERSION)))
        self.sym_dir = Path('upload-current')
        self.logger = utils.Logger(logfile=self.working_dir / 'setup.log')
        utils.LOG = self.logger

    def run(self):
        self.logger.info("Using {} as working dir".format(self.working_dir))
        self.logger.info("Current data format version is", utils.DATA_FORMAT_VERSION)

        if self._pre_flight_checks() > 1:
            exit(1)
        self._create()
        self._inform()

    def _pre_flight_checks(self):
        if (self.working_dir / self.sym_dir).exists():
            self.logger.error(
                "Working directory",
                self.working_dir,
                "already seems to have been setup. Please empty it and try again",
            )
            return 1
        return 0

    def _create(self):
        os.mkdir(self.working_dir / self.version_dir)
        os.mkdir(self.working_dir / self.version_dir / 'upload_1')
        open(self.working_dir / self.version_dir / "config.yaml", 'w').close()
        open(self.working_dir / self.version_dir / "assemblies.yaml", 'w').close()
        cwd = Path.cwd()
        os.chdir(self.working_dir)
        os.symlink(self.version_dir, self.sym_dir, target_is_directory=True)
        os.chdir(cwd)

    def _inform(self):
        self.logger.info("Your working environment has been set up in", (self.working_dir / self.sym_dir))
        self.logger.info("In there you will find dummy config.yaml and assemblies.yaml files")
        self.logger.info("You will need to edit these as described in the User Guide.")
        self.logger.info(
            "Then you can run collator like this (omit the -d argument if you are already in "
            + "that directory):\n      python -m xchemalign.collator -d",
            self.working_dir,
        )


def main():
    parser = argparse.ArgumentParser(description="setup")

    parser.add_argument("-d", "--dir", help="Working directory")

    args = parser.parse_args()

    s = Setup(args.dir)
    s.run()


if __name__ == "__main__":
    main()
