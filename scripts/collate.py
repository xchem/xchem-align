import argparse
import shutil

from xchemalign import utils
from xchemalign.collator import Collator


def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-c", "--config-file", default="config.yaml", help="Configuration file")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("-v", "--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()

    c = Collator(args.config_file, log_file=args.log_file, log_level=args.log_level)
    logger = c.logger
    logger.info("collator: ", str(args))

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
