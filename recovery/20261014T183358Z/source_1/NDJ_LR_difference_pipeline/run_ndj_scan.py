#!/usr/bin/env python3
import sys
from ndj_scan.cli import main

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr, flush=True)
        raise
