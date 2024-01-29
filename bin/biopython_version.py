#!/usr/bin/env python

import Bio
import sys


def biopython_version():
    return Bio.__version__


def main(args=None):
    print(biopython_version())


if __name__ == "__main__":
    sys.exit(main())
