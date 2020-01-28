#!/usr/bin/env python3

import logging
import sys
import os


from dmdpy.utility import utilities
import dmdpy.calculation


def main():
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .dmdpy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    # TODO have it receive argument for the number of cores, time to run, and a scratch directory

    try:
        logger.debug("Attempting to begin the calculation")
        c = dmdpy.calculation(cores=8)

    except:
        logger.exception("Check the error")
        logger.error("Error on the calculation")
        sys.exit(1)


if __name__ == "__main__":
    main()
