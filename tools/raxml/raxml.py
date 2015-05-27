#!/usr/bin/env python

import random
import sys
import shutil
import gd_util

################################################################################

if len(sys.argv) != 3:
    gd_util.die('Usage')

input, output = sys.argv[1:]
random.seed()

################################################################################

prog = 'raxmlHPC'

args = [ prog ]

## required: -s sequenceFileName -n outputFileName -m substitutionModel
## we supply -s, -n (they are not allowed from user)

args.append('-s')           # name of the alignment data file in PHYLIP format
args.append(input)

args.append('-n')           # name of the output file
args.append('fake')

## default options
args.append('-m')           # substitutionModel
args.append('GTRGAMMA')     # GTR + Optimization of substitution rates + GAMMA model of rate
                            # heterogeneity (alpha parameter will be estimated)

args.append('-N')           # number of alternative runs on distinct starting trees
args.append(1000)

args.append('-f')           # select algorithm
args.append('a')            # rapid Bootstrap analysis and search for
                            # best-scoring ML tree in one program run

args.append('-x')           # integer random seed and turn on rapid bootstrapping
args.append(random.randint(0,100000000000000))

args.append('-p')           # random seed for parsimony inferences
args.append(random.randint(0,100000000000000))

gd_util.run_program(prog, args)
shutil.copy2('RAxML_bipartitions.fake', output)
sys.exit(0)
