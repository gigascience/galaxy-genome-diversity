#!/usr/bin/env python

import sys
import gd_util

from Population import Population

################################################################################

if len(sys.argv) != 7:
    gd_util.die('Usage')

input, input_type, ind_arg, p1_input, p2_input, output = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(p1_input)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the first population that is not in the SNP table')

p2 = Population()
p2.from_population_file(p2_input)
if not p_total.is_superset(p2):
    gd_util.die('There is an individual in the second population that is not in the SNP table')

################################################################################

prog = 'offspring_heterozygosity'

args = [ prog ]
args.append(input)  # a Galaxy SNP table

for tag in p1.tag_list():
    column, name = tag.split(':')

    if input_type == 'gd_genotype':
        column  = int(column) - 2

    tag = '{0}:{1}:{2}'.format(column, 0, name)
    args.append(tag)

for tag in p2.tag_list():
    column, name = tag.split(':')

    if input_type == 'gd_genotype':
        column = int(column) - 2

    tag = '{0}:{1}:{2}'.format(column, 1, name)
    args.append(tag)

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

sys.exit(0)

