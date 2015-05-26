#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 6:
    gd_util.dir('Usage')

input, p1_input, output, input_type, ind_arg  = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(p1_input)

if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population that is not in the SNP table')

################################################################################

prog = 'aggregate'

args = [ prog ]
args.append(input)

if input_type == 'gd_snp':
    args.append(1)
elif input_type == 'gd_genotype':
    args.append(0)
else:
    die('unknown input type: {0}'.format(input_type))

columns = p1.column_list()

for column in sorted(columns):
    if input_type == 'gd_genotype':
        column = str(int(column) - 2)
    args.append(column)

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

