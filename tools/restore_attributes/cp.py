#!/usr/bin/env python

import shutil
import sys

if len(sys.argv) != 3:
    print >> sys.stderr, 'Usage: %s <src> <dst>' % sys.argv[0]
    sys.exit(1)

src, dst = sys.argv[1:3]
shutil.copy2(src, dst)
