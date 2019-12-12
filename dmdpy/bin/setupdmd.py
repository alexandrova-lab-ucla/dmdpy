#!/usr/bin/env python3

import logging
import sys


import dmdpy.utility.utilities as utilities
from dmdpy.setupjob import setupDMDjob


def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .dmdpy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    sdj = setupDMDjob()

    try:
        pass

    except:
        logger.critical("Exception encountered, quiting")
        sys.exit(1)


if __name__ == "__main__":
    main()

