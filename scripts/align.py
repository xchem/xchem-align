import argparse

from xchemalign import utils
from xchemalign.utils import Constants
from xchemalign.aligner import Aligner


def main():
    parser = argparse.ArgumentParser(description="aligner")

    parser.add_argument("-d", "--version-dir", required=True, help="Path to version dir")
    parser.add_argument("-m", "--metadata_file", default=Constants.METADATA_XTAL_FILENAME, help="Metadata YAML file")
    parser.add_argument("-x", "--xtalforms", help="Crystal forms YAML file")
    parser.add_argument("-a", "--assemblies", help="Crystal forms YAML file")

    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()
    print("aligner: ", args)

    logger = utils.Logger(logfile=args.log_file, level=args.log_level)

    a = Aligner(args.version_dir, args.metadata_file, args.xtalforms, args.assemblies, logger=logger)
    num_errors, num_warnings = a.validate()

    if not args.validate:
        if num_errors:
            print("There are errors, cannot continue")
            exit(1)
        else:
            a.run()


if __name__ == "__main__":
    main()
