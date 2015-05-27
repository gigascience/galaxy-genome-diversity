#!/usr/bin/env python

import gd_util
import os
import re
import shutil
import sys
from BeautifulSoup import BeautifulSoup
import gd_composite

################################################################################

def do_ped2geno(input, output):
    lines = []
    with open(input) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            lines.append(line.split())

    pair_map = {
        '0':{ '0':'9', '1':'9', '2':'9' },
        '1':{ '0':'1', '1':'2', '2':'1' },
        '2':{ '0':'1', '1':'1', '2':'0' }
    }
    with open(output, 'w') as ofh:
        for a_idx in xrange(6, len(lines[0]), 2):
            b_idx = a_idx + 1
            print >> ofh, ''.join(map(lambda line: pair_map[line[a_idx]][line[b_idx]], lines))

def do_map2snp(input, output):
    with open(output, 'w') as ofh:
        with open(input) as fh:
            for line in fh:
                elems = line.split()
                print >> ofh, '  {0} 11 0.002 2000 A T'.format(elems[1])

def make_ind_file(ind_file, input):
    pops = []
    name_map = []
    name_idx = 0

    ofh = open(ind_file, 'w')

    with open(input) as fh:
        soup = BeautifulSoup(fh)
        misc = soup.find('div', {'id': 'gd_misc'})
        populations = misc('ul')[0]

        i = 0
        for entry in populations:
            if i % 2 == 1:
                population_name = entry.contents[0].encode('utf8').strip().replace(' ', '_')
                pops.append(population_name)
                individuals = entry.ol('li')
                for individual in individuals:
                    individual_name = individual.string.encode('utf8').strip()
                    name_map.append(individual_name)
                    print >> ofh, 'ind_%s' % name_idx, 'M', population_name
                    name_idx += 1
            i += 1

    ofh.close()
    return pops, name_map

def make_par_file(par_file, geno_file, snp_file, ind_file, evec_file, eval_file):
    with open(par_file, 'w') as fh:
        print >> fh, 'genotypename: {0}'.format(geno_file)
        print >> fh, 'snpname: {0}'.format(snp_file)
        print >> fh, 'indivname: {0}'.format(ind_file)
        print >> fh, 'evecoutname: {0}'.format(evec_file)
        print >> fh, 'evaloutname: {0}'.format(eval_file)
        print >> fh, 'altnormstyle: NO'
        print >> fh, 'numoutevec: 2'

def do_smartpca(par_file):
    prog = 'smartpca'

    args = [ prog ]
    args.append('-p')
    args.append(par_file)

    stdoutdata, stderrdata = gd_util.run_program(prog, args)

    stats = []

    save_line = False
    for line in stdoutdata.split('\n'):
        if line.startswith(('## Average divergence', '## Anova statistics', '## Statistical significance')):
            stats.append('')
            save_line = True
        if line.strip() == '':
            save_line = False
        if save_line:
            stats.append(line)

    return '\n'.join(stats[1:])

def do_ploteig(evec_file, population_names):
    prog = 'gd_ploteig'

    args = [ prog ]
    args.append('-i')
    args.append(evec_file)
    args.append('-c')
    args.append('1:2')
    args.append('-p')
    args.append(':'.join(population_names))
    args.append('-x')

    gd_util.run_program(prog, args)

def do_eval2pct(eval_file, explained_file):
    prog = 'eval2pct'

    args = [ prog ]
    args.append(eval_file)

    with open(explained_file, 'w') as fh:
        gd_util.run_program(prog, args, stdout=fh)

def do_coords2admix(coords_file):
    prog = 'coords2admix'

    args = [ prog ]
    args.append(coords_file)

    with open('fake', 'w') as fh:
        gd_util.run_program(prog, args, stdout=fh)

    shutil.copy2('fake', coords_file)

ind_regex = re.compile('ind_([0-9]+)')

def fix_names(name_map, files):
    for file in files:
        tmp_filename = '%s.tmp' % file
        with open(tmp_filename, 'w') as ofh:
            with open(file) as fh:
                for line in fh:
                    line = line.rstrip('\r\n')
                    match = ind_regex.search(line)
                    if match:
                        idx = int(match.group(1))
                        old = 'ind_%s' % idx
                        new = name_map[idx].replace(' ', '_')
                        line = line.replace(old, new)
                    print >> ofh, line

        shutil.copy2(tmp_filename, file)
        os.unlink(tmp_filename)
        
################################################################################

if len(sys.argv) != 5:
    gd_util.die('Usage')

input, input_files_path, output, output_files_path = sys.argv[1:5]
gd_util.mkdir_p(output_files_path)

################################################################################

ped_file = os.path.join(input_files_path, 'admix.ped')
geno_file = os.path.join(output_files_path, 'admix.geno')
do_ped2geno(ped_file, geno_file)

################################################################################

map_file = os.path.join(input_files_path, 'admix.map')
snp_file = os.path.join(output_files_path, 'admix.snp')
do_map2snp(map_file, snp_file)

################################################################################

ind_file = os.path.join(output_files_path, 'admix.ind')
population_names, name_map = make_ind_file(ind_file, input)

################################################################################

par_file = os.path.join(output_files_path, 'par.admix')
evec_file = os.path.join(output_files_path, 'coordinates.txt')
eval_file = os.path.join(output_files_path, 'admix.eval')
make_par_file(par_file, geno_file, snp_file, ind_file, evec_file, eval_file)

################################################################################

smartpca_stats = do_smartpca(par_file)
fix_names(name_map, [ind_file, evec_file])

################################################################################

do_ploteig(evec_file, population_names)
plot_file = 'coordinates.txt.1:2.{0}.pdf'.format(':'.join(population_names))
output_plot_file = os.path.join(output_files_path, 'PCA.pdf')
shutil.copy2(plot_file, output_plot_file)
os.unlink(plot_file)

################################################################################

do_eval2pct(eval_file, os.path.join(output_files_path, 'explained.txt'))
os.unlink(eval_file)

################################################################################

do_coords2admix(evec_file)

################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('PCA Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='PCA.pdf', value='PCA.pdf', display_type=display_file)
out_evec = gd_composite.Parameter(name='coordinates.txt', value='coordinates.txt', display_type=display_file)
out_explained = gd_composite.Parameter(name='explained.txt', value='explained.txt', display_type=display_file)

evec_prefix = 'coordinates.txt.1:2.{0}'.format(':'.join(population_names))
ps_file = '{0}.ps'.format(evec_prefix)
xtxt_file = '{0}.xtxt'.format(evec_prefix)

os.unlink(os.path.join(output_files_path, ps_file))
os.unlink(os.path.join(output_files_path, xtxt_file))

info_page.add_output_parameter(out_pdf)
info_page.add_output_parameter(out_evec)
info_page.add_output_parameter(out_explained)

in_admix = gd_composite.Parameter(name='par.admix', value='par.admix', display_type=display_file)
in_geno = gd_composite.Parameter(name='admix.geno', value='admix.geno', display_type=display_file)
in_snp = gd_composite.Parameter(name='admix.snp', value='admix.snp', display_type=display_file)
in_ind = gd_composite.Parameter(name='admix.ind', value='admix.ind', display_type=display_file)

info_page.add_input_parameter(in_admix)
info_page.add_input_parameter(in_geno)
info_page.add_input_parameter(in_snp)
info_page.add_input_parameter(in_ind)

misc_stats = gd_composite.Parameter(description='Stats<p/><pre>\n{0}\n</pre>'.format(smartpca_stats), display_type=display_value)

info_page.add_misc(misc_stats)

with open (output, 'w') as ofh:
    print >> ofh, info_page.render()

sys.exit(0)

