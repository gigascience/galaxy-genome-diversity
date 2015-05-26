#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 12:
    gd_util.die('Usage')

input, p1_input, p2_input, input_type, data_source, min_total_count, discard_fixed, output, shuffles, p0_input, ind_arg = sys.argv[1:]

try:
    shuffle_count = int(shuffles)
except:
    shuffle_count = 0

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

p0 = None
if shuffle_count > 0:
    p0 = Population()
    p0.from_population_file(p0_input)
    if not p_total.is_superset(p0):
        gd_util.die('There is an individual in population 0 that is not in the SNP table')

################################################################################

prog = 'Fst_ave'

args = [ prog ]
args.append(input)
args.append(data_source)
args.append(min_total_count)
args.append(discard_fixed)
args.append(shuffles)

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

if p0 is not None:
    columns = p0.column_list()
    for column in columns:
        if input_type == 'gd_genotype':
            column = int(column) - 2
        args.append('{0}:0'.format(column))

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

