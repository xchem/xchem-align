import argparse
import shutil
from pathlib import Path

from xchemalign import utils
from xchemalign import setup
from xchemalign.collator import Collator, _check_working_dir


def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-d", "--dir", help="Working directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("-v", "--validate", action="store_true", help="Only perform validation")
    parser.add_argument("--no-git-info", action="store_false", help="Don't add GIT info to metadata")

    args = parser.parse_args()

    if args.dir:
        working_dir = Path(args.dir)
    else:
        working_dir = Path.cwd()

    if _check_working_dir(working_dir):
        print("Working dir does not seem to have been initialised - missing 'upload_current' symlink")
        inp = input("Do you want the working dir to be initialised? (Y/N)")
        if inp == "Y" or inp == "y":
            print("Initialising working dir")
            s = setup.Setup(args.dir)
            s.run()
        exit(1)

    c = Collator(working_dir, log_file=args.log_file, log_level=args.log_level, include_git_info=args.no_git_info)

    logger = c.logger
    logger.info("collator: ", str(args))
    utils.LOG = logger

    meta = c.validate()

    if not args.validate:
        if meta is None or len(logger.errors) > 0:
            logger.error("There are errors, cannot continue")
            logger.report()
            logger.close()
            exit(1)
        else:
            c.run(meta)
            # write a summary of errors and warnings
            logger.report()
            logger.close()
            if logger.logfilename:
                to_path = c.output_path / c.version_dir / 'collator.log'
                print("copying log file", logger.logfilename, "to", to_path)
                f = shutil.copy2(logger.logfilename, to_path)
                if not f:
                    print("Failed to copy log file {} to {}".format(logger.logfilename, to_path))


if __name__ == "__main__":
    main()
