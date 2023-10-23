import argparse

from xchemalign import utils
from xchemalign.collator import Collator


def main():
    parser = argparse.ArgumentParser(description="collator")

    parser.add_argument("-c", "--config-file", default="config.yaml", help="Configuration file")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level")
    parser.add_argument("-v", "--validate", action="store_true", help="Only perform validation")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("collator: ", args)

    c = Collator(args.config_file, logger=logger)

    meta, num_errors, num_warnings = c.validate()

    if not args.validate:
        if meta is None or num_errors:
            print("There are errors, cannot continue")
            exit(1)
        else:
            c.run(meta)


if __name__ == "__main__":
    main()
