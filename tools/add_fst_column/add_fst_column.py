#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 13:
    gd_util.die('Usage')

input, p1_input, p2_input, input_type, data_source, min_reads, min_qual, retain, discard_fixed, biased, output, ind_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(p1_input)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in population 1 that is not in the SNP table')

p2 = Population()
p2.from_population_file(p2_input)
if not p_total.is_superset(p2):
    gd_util.die('There is an individual in population 2 that is not in the SNP table')

################################################################################

prog = 'Fst_column'

args = [ prog ]
args.append(input)
args.append(data_source)
args.append(min_reads)
args.append(min_qual)
args.append(retain)
args.append(discard_fixed)
args.append(biased)

columns = p1.column_list()
for column in columns:
    if input_type == 'gd_genotype':
        column = int(column) - 2
    args.append('{0}:1'.format(column))

columns = p2.column_list()
for column in columns:
    if input_type == 'gd_genotype':
        column = int(column) - 2
    args.append('{0}:2'.format(column))

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

