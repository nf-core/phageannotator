#!/usr/bin/env python

import numpy as np
import sys


def numpy_version():
    return np.__version__


def main(args=None):
    print(numpy_version())


if __name__ == "__main__":
    sys.exit(main())
