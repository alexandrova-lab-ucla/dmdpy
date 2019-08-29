#!/usr/bin/env python3

import logging


import dmdpy.utilities.utilities as utilities
from dmdpy.setupdmdjob import setupDMDjob


if __name__ == "__main__":

    utilities.load_logger_config()
    logger = logging.getLogger(__name__)

    try:
        sdj = setupDMDjob()

    except:
        logger.critical("Exception encountered, quiting")
        raise
