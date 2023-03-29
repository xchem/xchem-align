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

from . import utils, processor, validator


class Aligner(processor.Processor):

    def validate(self):
        v = validator.Validator(self.base_dir, self.input_dirs, self.output_dir, self.target_name, logger=self.logger)
        warnings, errors = v.validate_paths()
        return warnings, errors

    def run(self):
        self.logger.info('Running aligner...')
        v_dir = self.read_versions()
        if not v_dir:
            self.logger.error('Error with version dir. Please fix and try again.')
            return None
        self.logger.info('Using version dir {}'.format(v_dir))

        meta = self.read_metadata(v_dir)

        self._perform_alignments(self.config, meta, self.meta_history, [], self.output_dir)

    def _perform_alignments(self, config, meta, meta_history, panddas_dirs, output_dir):
        self.logger.info('Performing alignments (well, not really)')


def main():

    parser = argparse.ArgumentParser(description='aligner')

    parser.add_argument('-c', '--config-file', default='config.yaml', help="Configuration file")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("aligner: ", args)

    logger = utils.Logger(logfile=args.log_file, level=args.log_level)

    a = Aligner(args.config_file, logger=logger)

    warnings, errors = a.validate()

    if not args.validate:
        if errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            a.run()


if __name__ == "__main__":
    main()

