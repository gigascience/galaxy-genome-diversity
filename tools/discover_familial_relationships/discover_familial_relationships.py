#!/usr/bin/env python

import sys
import gd_util

from Population import Population

################################################################################

if len(sys.argv) != 6:
    gd_util.die('Usage')

input, input_type, ind_arg, pop_input, output = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(pop_input)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population that is not in the SNP table')

################################################################################

prog = 'kinship_prep'

args = [ prog ]
args.append(input)  # a Galaxy SNP table
args.append(0)      # required number of reads for each individual to use a SNP
args.append(0)      # required genotype quality for each individual to use a SNP
args.append(0)      # minimum spacing between SNPs on the same scaffold

for tag in p1.tag_list():
    if input_type == 'gd_genotype':
        column, name = tag.split(':')
        tag = '{0}:{1}'.format(int(column) - 2, name)
    args.append(tag)

gd_util.run_program(prog, args)

# kinship.map
# kinship.ped
# kinship.dat

################################################################################

prog = 'king'

args = [ prog ]
args.append('-d')
args.append('kinship.dat')
args.append('-p')
args.append('kinship.ped')
args.append('-m')
args.append('kinship.map')
args.append('--kinship')

gd_util.run_program(prog, args)

# king.kin

################################################################################

valid_header = 'FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n'

with open('king.kin') as fh:
    header = fh.readline()
    if header != valid_header:
        gd_util.die('crap')

    with open(output, 'w') as ofh:

        for line in fh:
            elems = line.split('\t')
            if len(elems) != 10:
                gd_util.die('crap')

            x = elems[1]
            y = elems[2]
            z = elems[8]

            f = float(z)

            message = ''

            if f > 0.354:
                message = 'duplicate or MZ twin'
            elif f >= 0.177:
                message = '1st degree relatives'
            elif f >= 0.0884:
                message = '2nd degree relatives'
            elif f >= 0.0442:
                message = '3rd degree relatives'

            print >> ofh, '\t'.join([x, y, z, message])

################################################################################

sys.exit(0)

