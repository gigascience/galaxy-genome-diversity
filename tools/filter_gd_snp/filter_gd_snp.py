#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

def convert_percent(string_value):
    if string_value.endswith('%'):
        val = convert_non_negative_int(string_value[:-1])
        if val > 100:
            print >> sys.stderr, 'percentage: "%d" > 100' % val
            sys.exit(1)
        val = val * -1
    else:
        val = convert_non_negative_int(string_value)

    return str(val)

def convert_non_negative_int(string_value):
    try:
        val = int(string_value)
    except:
        print >> sys.stderr, '"%s" is not an integer' % string_value
        sys.exit(1)

    if val < 0:
        print >> sys.stderr, '"%d" is negative' % val
        sys.exit(1)

    return val

################################################################################

if len(sys.argv) != 13:
    gd_util.die('Usage')

input, output, ref_chrom_col, min_spacing, lo_genotypes, p1_input, input_type, lo_coverage, hi_coverage, low_ind_cov, low_quality, ind_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(p1_input)

if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population that is not in the SNP table')

lo_coverage = convert_percent(lo_coverage)
hi_coverage = convert_percent(hi_coverage)

if input_type == 'gd_snp':
    type_arg = 1
elif input_type == 'gd_genotype':
    type_arg = 0
else:
    gd_util.die('unknown input_type: {0}'.format(input_type))

################################################################################

prog = 'filter_snps'

args = [ prog ]
args.append(input)          # file containing a Galaxy table
args.append(type_arg)       # 1 for a gd_snp file, 0 for gd_genotype
args.append(lo_coverage)    # lower bound on total coverage (< 0 means interpret as percentage)
args.append(hi_coverage)    # upper bound on total coveraae (< 0 means interpret as percentage)
args.append(low_ind_cov)    # lower bound on individual coverage
args.append(low_quality)    # lower bound on individual quality value
args.append(lo_genotypes)   # lower bound on the number of defined genotypes
args.append(min_spacing)    # lower bound on the spacing between SNPs
args.append(ref_chrom_col)  # reference-chromosome column (base-1); ref position in next column

columns = p1.column_list()
for column in sorted(columns):
    args.append(column)     # the starting columns (base-1) for the chosen individuals

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

