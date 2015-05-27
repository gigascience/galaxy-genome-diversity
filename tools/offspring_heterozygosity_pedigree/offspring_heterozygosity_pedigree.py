#!/usr/bin/env python

import sys
import gd_util

from Population import Population

def load_and_check_pop(file, total_pop, name):
        p = Population()
        p.from_population_file(file)
        if not total_pop.is_superset(p):
            gd_util.die('There is an individual in the {0} that is not in the SNP table'.format(name))
        return p

def append_breeders_from_file(the_list, filename, kind):
        with open(filename) as fh:
            for line in fh:
                elems = line.split()
                breeder = elems[0].rstrip('\r\n')
                the_list.append('{0}:{1}'.format(kind, breeder))

################################################################################

if len(sys.argv) != 9:
    gd_util.die('Usage')

input, input_type, pedigree, ind_arg, founders, b1_input, b2_input, output = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

f1 = load_and_check_pop(founders, p_total, 'founders')

################################################################################

prog = 'offspring_heterozygosity2'

args = [ prog ]
args.append(input)      # a Galaxy SNP table
args.append(pedigree)   # a pedigree, where the SNP table is for the founders

for tag in f1.tag_list():
    column, name = tag.split(':')
    if type == 'gd_genotype':
        column = int(column) - 2
    tag = 'founder:{0}:{1}'.format(column, name)
    args.append(tag)

append_breeders_from_file(args, b1_input, 0)
append_breeders_from_file(args, b2_input, 1)

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

sys.exit(0)

