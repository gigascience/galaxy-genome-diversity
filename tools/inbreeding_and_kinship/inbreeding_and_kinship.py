#!/usr/bin/env python

import sys
import gd_util

################################################################################

if len(sys.argv) != 6:
    gd_util.die('Usage')

ped_input, ind_input, computed_value, output, kinship_input = sys.argv[1:]

################################################################################

prog = 'inbreed'

args = [ prog ]
args.append(ped_input)      # pedigree
args.append(ind_input)      # specified individuals (e.g.,,potential breeding population)
args.append(kinship_input)  # kinships of founders
args.append(computed_value) # 0 = inbreedng coefficients, 1 = kinships, 2 = mean kinships

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

sys.exit(0)

