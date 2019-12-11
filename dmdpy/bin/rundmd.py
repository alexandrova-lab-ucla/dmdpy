#!/usr/bin/env python3

import logging
import sys


from dmdpy.utility import utilities
import dmdpy.calculation

def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .dmdpy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    commands = {
        "1" : {
            "continue": True,
            "Start time": 0,
            "Max time": 50,
        }
    }

    c = dmdpy.calculation(8, commands)


if __name__ == "__main__":
    main()