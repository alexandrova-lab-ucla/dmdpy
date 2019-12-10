#!/usr/bin/env python3

import logging
import sys


import dmdpy.utilities.utilities as utilities
from dmdpy.setupjob import setupDMDjob

def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .turbopy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    try:
        sdj = setupDMDjob()

    except:
        logger.critical("Exception encountered, quiting")
        raise


if __name__ == "__main__":
    main()

