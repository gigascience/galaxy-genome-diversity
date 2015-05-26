#!/usr/bin/env python

import gd_util
import os
import sys
from Population import Population
import gd_composite

################################################################################

if len(sys.argv) < 7:
    gd_util.die('Usage')

input, data_source, output, extra_files_path, ind_arg = sys.argv[1:6]

population_info = []
p1_input = None
all_individuals = False

for arg in sys.argv[6:]:
    if arg == 'all_individuals':
        all_individuals = True
    elif len(arg) > 12 and arg[:12] == 'individuals:':
        p1_input = arg[12:]
    elif len(arg) > 11 and arg[:11] == 'population:':
        file, name = arg[11:].split(':', 1)
        population_info.append((file, name))

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

################################################################################

gd_util.mkdir_p(extra_files_path)

################################################################################

prog = 'coverage'

args = [ prog ]
args.append(input)
args.append(data_source)

user_coverage_file = os.path.join(extra_files_path, 'coverage.txt')
args.append(user_coverage_file)

population_list = []

if all_individuals:
    tags = p_total.tag_list()
elif p1_input is not None:
    p1 = Population()
    this_pop = Population()
    this_pop.from_population_file(p1_input)
    population_list.append(this_pop)
    p1.from_population_file(p1_input)
    if not p_total.is_superset(p1):
        gd_util.die('There is an individual in the population that is not in the SNP table')
    tags = p1.tag_list()
else:
    tags = []
    for population_file, population_name in population_info:
        population = Population()
        this_pop = Population()
        this_pop.from_population_file(population_file)
        population_list.append(this_pop)
        population.from_population_file(population_file)
        if not p_total.is_superset(population):
            gd_util.die('There is an individual in the {} population that is not in the SNP table'.format(population_name))
        columns = population.column_list()
        for column in columns:
            tags.append('{0}:{1}'.format(column, population_name))

for tag in tags:
    args.append(tag)

## text output
coverage_file = 'coverage.txt'
with open(coverage_file, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

## graphical output
coverage2_file = 'coverage2.txt'
with open(coverage_file) as fh, open(coverage2_file, 'w') as ofh:
    for line in fh:
        line = line.rstrip('\r\n')
        elems = line.split('\t')
        name = elems.pop(0)
        values = [ elems[0] ]
        for idx in range(1, len(elems)):
            val = str(float(elems[idx]) - float(elems[idx-1]))
            values.append(val)
        print >> ofh, '{0}\t{1}'.format(name, '\t'.join(values))

################################################################################

prog = 'Rscript'

args = [ prog ]

_realpath = os.path.realpath(__file__)
_script_dir = os.path.dirname(_realpath)
r_script_file = os.path.join(_script_dir, 'coverage_plot.r')
args.append(r_script_file)

pdf_file = os.path.join(extra_files_path, 'coverage.pdf')
args.append(pdf_file)

gd_util.run_program(prog, args)

################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('Coverage distributions Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='coverage.pdf', value='coverage.pdf', display_type=display_file)
out_txt = gd_composite.Parameter(name='coverage.txt', value='coverage.txt', display_type=display_file)

info_page.add_output_parameter(out_pdf)
info_page.add_output_parameter(out_txt)

if data_source == '0':
    data_source_value = 'sequence coverage'
elif data_source == '1':
    data_source_value = 'estimated genotype'

in_data_source = gd_composite.Parameter(description='Data source', value=data_source_value, display_type=display_value)

info_page.add_input_parameter(in_data_source)

if population_list:
    misc_populations =  gd_composite.Parameter(name='Populations', value=population_list, display_type=gd_composite.DisplayPopulationList())
    info_page.add_misc(misc_populations)
else:
    misc_individuals = gd_composite.Parameter(name='Individuals', value=tags, display_type=gd_composite.DisplayTagList())
    info_page.add_misc(misc_individuals)

with open (output, 'w') as ofh:
    print >> ofh, info_page.render()

sys.exit(0)

