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


snp_file, snp_ext, snp_arg, indiv_input, annotation_input, cov_file, cov_ext, cov_arg, min_coverage, output = sys.argv[1:]

p_snp = load_pop(snp_file, snp_arg)
p_cov = load_pop(cov_file, cov_arg)

if indiv_input == '/dev/null':
    if p_snp is not None:
        p_ind = p_snp
    elif p_cov is not None:
        p_ind = p_cov
    else:
        p_ind = None
    order_p_ind = True
else:
    p_ind = Population()
    p_ind.from_population_file(indiv_input)
    order_p_ind = False

## p ind must be from either p_snp or p_cov
if p_snp is not None and p_cov is not None:
    if not (p_snp.is_superset(p_ind) or p_cov.is_superset(p_ind)):
        gd_util.die('There is an individual in the population individuals that is not in the SNP/Genotype or Coverage table')
elif p_snp is not None:
    if not p_snp.is_superset(p_ind):
        gd_util.die('There is an individual in the population individuals that is not in the SNP/Genotype table')
elif p_cov is not None:
    if not p_cov.is_superset(p_ind):
        gd_util.die('There is an individual in the population individuals that is not in the Coverage table')


################################################################################

prog = 'mito_draw'

args = [ prog ]
args.append(snp_file)
args.append(cov_file)
args.append(annotation_input)
args.append(min_coverage)

if order_p_ind:
    for column in sorted(p_ind.column_list()):
        individual = p_ind.individual_with_column(column)
        name = individual.name.split()[0]
        args.append('{0}:{1}:{2}'.format(0, column, name))
else:
    append_tags(args, p_ind, 'gd_indivs', 0)

append_tags(args, p_snp, snp_ext, 1)
append_tags(args, p_cov, cov_ext, 2)

with open('Ji.spec', 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'varplot'

args = [ prog ]
args.append('-w')
args.append(3)
args.append('-s')
args.append(0.3)
args.append('-g')
args.append(0.2)
args.append('Ji.spec')

with open('Ji.svg', 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'convert'

args = [ prog ]
args.append('-density')
args.append(200)
args.append('-resize')
args.append('140%')
args.append('Ji.svg')
args.append('-compress')
args.append('zip')
args.append('tiff:{0}'.format(output))

gd_util.run_program(prog, args)
sys.exit(0)
