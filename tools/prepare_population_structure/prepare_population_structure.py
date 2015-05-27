#!/usr/bin/env python

import gd_util
import os
import shutil
import sys
from Population import Population
import gd_composite

################################################################################

def do_import(filename, files_path, min_reads, min_qual, min_spacing, using_info, population_list):
    info_page = gd_composite.InfoPage()
    info_page.set_title('Prepare to look for population structure Galaxy Composite Dataset')

    display_file = gd_composite.DisplayFile()
    display_value = gd_composite.DisplayValue()

    out_ped = gd_composite.Parameter(name='admix.ped', value='admix.ped', display_type=display_file)
    out_map = gd_composite.Parameter(name='admix.map', value='admix.map', display_type=display_file)
    out_use = gd_composite.Parameter(description=using_info, display_type=display_value)

    info_page.add_output_parameter(out_ped)
    info_page.add_output_parameter(out_map)
    info_page.add_output_parameter(out_use)

    in_min_reads = gd_composite.Parameter(description='Minimum reads covering a SNP, per individual', value=min_reads, display_type=display_value)
    in_min_qual = gd_composite.Parameter(description='Minimum quality value, per individual', value=min_qual, display_type=display_value)
    in_min_spacing = gd_composite.Parameter(description='Minimum spacing between SNPs on the same scaffold', value=min_spacing, display_type=display_value)

    info_page.add_input_parameter(in_min_reads)
    info_page.add_input_parameter(in_min_qual)
    info_page.add_input_parameter(in_min_spacing)

    misc_populations = gd_composite.Parameter(name='Populations', value=population_list, display_type=gd_composite.DisplayPopulationList())
    info_page.add_misc(misc_populations)

    with open(filename, 'w') as ofh:
        print >> ofh, info_page.render()

################################################################################

if len(sys.argv) < 10:
    gd_util.die('Usage')

# parse command line
input_snp_filename, input_type, min_reads, min_qual, min_spacing, output_filename, output_files_path, ind_arg = sys.argv[1:9]
args = sys.argv[9:]

population_files = []
all_individuals = False

for arg in args:
    if arg == 'all_individuals':
        all_individuals = True
    elif len(arg) > 11 and arg[:11] == 'population:':
        file, name = arg[11:].split(':', 1)
        population_files.append((file, name))

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

individual_population = {}
population_list = []

if all_individuals:
    p1 = p_total
    p1.name = 'All Individuals'
    population_list.append(p1)
else:
    p1 = Population()
    for file, name in population_files:
        this_pop = Population(name)
        this_pop.from_population_file(file)
        population_list.append(this_pop)

        for tag in this_pop.tag_list():
            if tag not in individual_population:
                individual_population[tag] = name

        # add individuals from this file to p1
        p1.from_population_file(file)


if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population that is not in the SNP table')

################################################################################

prog = 'admix_prep'

args = [ prog ]
args.append(input_snp_filename)
args.append(min_reads)
args.append(min_qual)
args.append(min_spacing)

for tag in p1.tag_list():
    if input_type == 'gd_genotype':
        column, name = tag.split(':', 1)
        tag = '{0}:{1}'.format(int(column) - 2, name)
    args.append(tag)

stdoutdata, stderrdata = gd_util.run_program(prog, args)
using_info = stdoutdata.rstrip('\r\n')

################################################################################

gd_util.mkdir_p(output_files_path)

output_ped_filename = os.path.join(output_files_path, 'admix.ped')
output_map_filename = os.path.join(output_files_path, 'admix.map')
shutil.copy2('admix.ped', output_ped_filename)
shutil.copy2('admix.map', output_map_filename)

do_import(output_filename, output_files_path, min_reads, min_qual, min_spacing, using_info, population_list)
sys.exit(0)

