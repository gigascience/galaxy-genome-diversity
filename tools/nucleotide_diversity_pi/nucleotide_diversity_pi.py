#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 7:
    gd_util.die('Usage')

gd_saps_file, gd_snps_file, covered_intervals_file, gd_indivs_file, output_file, ind_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(gd_indivs_file)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population individuals that is not in the SNP table')

################################################################################

prog = 'get_pi'

args = [ prog ]
args.append(gd_saps_file)
args.append(gd_snps_file)
args.append(covered_intervals_file)

columns = p1.column_list()
for column in columns:
    args.append(column)

with open(output_file, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

