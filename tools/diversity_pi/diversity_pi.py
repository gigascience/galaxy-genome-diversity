#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

def load_pop(file, wrapped_dict):
    if file == '/dev/null':
        pop = None
    else:
        pop = Population()
        pop.from_wrapped_dict(wrapped_dict)
    return pop

def append_tags(the_list, p, p_type, val):
    if p is None:
        return
    for tag in p.tag_list():
        column, name = tag.split(':')
        if p_type == 'gd_genotype':
            column = int(column) - 2
        the_list.append('{0}:{1}:{2}'.format(val, column, name))

################################################################################

if len(sys.argv) != 11:
    gd_util.die('Usage')

snp_input, snp_ext, snp_arg, cov_input, cov_ext, cov_arg, indiv_input, min_coverage, req_thresh, output = sys.argv[1:]

p_snp = load_pop(snp_input, snp_arg)
p_cov = load_pop(cov_input, cov_arg)

p_ind = Population()
p_ind.from_population_file(indiv_input)

if not p_snp.is_superset(p_ind):
  gd_util.die('There is an individual in the population individuals that is not in the SNP/Genotype table')

if p_cov is not None and (not p_cov.is_superset(p_ind)):
  gd_util.die('There is an individual in the population individuals that is not in the Coverage table')

################################################################################

prog = 'mito_pi'

args = [ prog ]
args.append(snp_input)
args.append(cov_input)
args.append(min_coverage)
args.append(req_thresh)

append_tags(args, p_ind, 'gd_indivs', 0)
append_tags(args, p_snp, snp_ext, 1)
append_tags(args, p_cov, cov_ext, 2)

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

