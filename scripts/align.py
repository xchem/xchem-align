import argparse
import shutil

from pathlib import Path

from xchemalign import utils
from xchemalign.utils import Constants
from xchemalign.aligner import Aligner


def main():
    parser = argparse.ArgumentParser(description="aligner")

    parser.add_argument("-d", "--version-dir", required=True, help="Path to version dir")
    parser.add_argument(
        "-m", "--metadata_file", default=Constants.METADATA_XTAL_FILENAME.format(""), help="Metadata YAML file"
    )
    parser.add_argument("-a", "--assemblies", help="Assemblies YAML file")

    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()

    if args.log_file:
        log = args.log_file
    else:
        log = str(Path(args.version_dir).parent / 'aligner.log')
        print("Using {} for log file".format(str(log)))

    logger = utils.Logger(logfile=log, level=args.log_level)
    logger.info("aligner: ", args)

    a = Aligner(args.version_dir, args.metadata_file, args.assemblies, logger=logger)
    num_errors, num_warnings = a.validate()

    if not args.validate:
        if num_errors:
            logger.error("There are errors, cannot continue")
            exit(1)
        else:
            a.run()
            # write a summary of errors and warnings
            logger.report()
            logger.close()
            if logger.logfilename:
                to_path = a.version_dir / 'aligner.log'
                print("copying log file", logger.logfilename, "to", to_path)
                f = shutil.copy2(logger.logfilename, to_path)
                if not f:
                    print("Failed to copy log file {} to {}".format(logger.logfilename, to_path))


if __name__ == "__main__":
    main()
