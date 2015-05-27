#!/usr/bin/env python

import gd_util
import os
import sys
from Population import Population
import gd_composite

################################################################################

if len(sys.argv) != 12:
    gd_util.die('Usage')

input, output, extra_files_path, input_type, data_source_arg, minimum_coverage, minimum_quality, p1_input, dbkey, draw_tree_options, ind_arg = sys.argv[1:]

if input_type == 'gd_snp':
    if data_source_arg == 'sequence_coverage':
        data_source = 0
    elif data_source_arg == 'estimated_genotype':
        data_source = 1
    else:
        gd_util.die('Unsupported data_source: {0}'.format(data_source_arg))
elif input_type == 'gd_genotype':
        data_source = 1
        minimum_coverage = 0
        minimum_quality = 0
else:
    gd_util.die('Unsupported input_type:: {0}'.format(input_type))

# note: TEST THIS
if dbkey in ['', '?', 'None']:
    dbkey = 'none'

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

if p1_input == "all_individuals":
    tags = p_total.tag_list()
else:
    p1 = Population()
    p1.from_population_file(p1_input)
    if not p_total.is_superset(p1):
        gd_util.die('There is an individual in the population that is not in the SNP table')
    tags = p1.tag_list()

################################################################################

gd_util.mkdir_p(extra_files_path)
phylip_outfile = os.path.join(extra_files_path, 'distance_matrix.phylip')
newick_outfile = os.path.join(extra_files_path, 'phylogenetic_tree.newick')
ps_outfile = 'tree.ps'
pdf_outfile = os.path.join(extra_files_path, 'tree.pdf')
informative_snp_file = os.path.join(extra_files_path, 'informative_snps.txt')
mega_distance_matrix_file = os.path.join(extra_files_path, 'mega_distance_matrix.txt')

################################################################################

prog = 'dist_mat'

args = [ prog ]
args.append(input)
args.append(minimum_coverage)
args.append(minimum_quality)
args.append(dbkey)
args.append(data_source)
args.append(informative_snp_file)
args.append(mega_distance_matrix_file)

for tag in tags:
    if input_type == 'gd_genotype':
        column, name = tag.split(':')
        tag = '{0}:{1}'.format(int(column) - 2, name)
    args.append(tag)

with open(phylip_outfile, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'quicktree'

args = [ prog ]
args.append('-in')
args.append('m')
args.append('-out')
args.append('t')
args.append(phylip_outfile)

with open(newick_outfile, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'draw_tree'

args = [ prog ]

if draw_tree_options:
    args.append(draw_tree_options)

args.append(newick_outfile)

with open(ps_outfile, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'ps2pdf'

args = [ prog ]
args.append('-dPDFSETTINGS=/prepress')
args.append(ps_outfile)
args.append(pdf_outfile)

gd_util.run_program(prog, args)

################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('Phylogenetic tree Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='tree.pdf', value='tree.pdf', display_type=display_file)
out_newick = gd_composite.Parameter(value='phylogenetic_tree.newick', name='phylogenetic tree (newick)', display_type=display_file)
out_phylip = gd_composite.Parameter(value='distance_matrix.phylip', name='Phylip distance matrix', display_type=display_file)
out_mega = gd_composite.Parameter(value='mega_distance_matrix.txt', name='Mega distance matrix', display_type=display_file)
out_snps = gd_composite.Parameter(value='informative_snps.txt', name='informative SNPs', display_type=display_file)

info_page.add_output_parameter(out_pdf)
info_page.add_output_parameter(out_newick)
info_page.add_output_parameter(out_phylip)
info_page.add_output_parameter(out_mega)
info_page.add_output_parameter(out_snps)

in_min_cov = gd_composite.Parameter(description='Minimum coverage', value=minimum_coverage, display_type=display_value)
in_min_qual = gd_composite.Parameter(description='Minimum quality', value=minimum_quality, display_type=display_value)

include_ref_value = 'no'
if dbkey != 'none':
    include_ref_value = 'yes'

in_include_ref = gd_composite.Parameter(description='Include reference sequence', value=include_ref_value, display_type=display_value)

if data_source == 0:
    data_source_value = 'sequence coverage'
elif data_source == 1:
    data_source_value = 'estimated genotype'

in_data_source = gd_composite.Parameter(description='Data source', value=data_source_value, display_type=display_value)

branch_type_value = 'square'
if 'd' in draw_tree_options:
    branch_type_value = 'diagonal'

in_branch_type = gd_composite.Parameter(description='Branch type', value=branch_type_value, display_type=display_value)

branch_scale_value = 'yes'
if 's' in draw_tree_options:
    branch_scale_value = 'no'

in_branch_scale = gd_composite.Parameter(description='Draw branches to scale', value=branch_scale_value, display_type=display_value)

branch_length_value = 'yes'
if 'b' in draw_tree_options:
    branch_length_value = 'no'

in_branch_length = gd_composite.Parameter(description='Show branch lengths', value=branch_length_value, display_type=display_value)

tree_layout_value = 'horizontal'
if 'v' in draw_tree_options:
    tree_layout_value = 'vertical'

in_tree_layout = gd_composite.Parameter(description='Tree layout', value=tree_layout_value, display_type=display_value)

info_page.add_input_parameter(in_min_cov)
info_page.add_input_parameter(in_min_qual)
info_page.add_input_parameter(in_include_ref)
info_page.add_input_parameter(in_data_source)
info_page.add_input_parameter(in_branch_type)
info_page.add_input_parameter(in_branch_scale)
info_page.add_input_parameter(in_branch_length)
info_page.add_input_parameter(in_tree_layout)

misc_individuals = gd_composite.Parameter(name='Individuals', value=tags, display_type=gd_composite.DisplayTagList())

info_page.add_misc(misc_individuals)


with open(output, 'w') as ofh:
    print >> ofh, info_page.render()

################################################################################

sys.exit(0)

