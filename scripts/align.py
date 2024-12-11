import argparse
import shutil

from pathlib import Path

from xchemalign import utils
from xchemalign.utils import Constants
from xchemalign.aligner import Aligner


def main():
    parser = argparse.ArgumentParser(description="aligner")

    parser.add_argument("-d", "--dir", help="Working directory")

    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()

    if args.dir:
        log = str(Path(args.dir) / 'aligner.log')
    else:
        log = 'aligner.log'

    logger = utils.Logger(logfile=log, level=args.log_level)
    logger.info("aligner: ", args)
    logger.info("Using {} for log file".format(str(log)))
    utils.LOG = logger

    a = Aligner(args.dir, logger=logger)
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
