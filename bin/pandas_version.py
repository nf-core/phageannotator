#!/usr/bin/env python

import pandas as pd
import sys


def pandas_version():
    return pd.__version__


def main(args=None):
    print(pandas_version())


if __name__ == "__main__":
    sys.exit(main())
