#!/usr/bin/env python3

import logging
import sys
import os
import json


from dmdpy.utility import utilities
import dmdpy.calculation

def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .dmdpy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    # commands = {
    #     "1" : {
    #         "Time" : 50
    #     }
    # }

    c = dmdpy.calculation(cores=8)


if __name__ == "__main__":
    main()