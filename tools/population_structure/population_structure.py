#!/usr/bin/env python

import errno
import os
import shutil
import subprocess
import sys
from BeautifulSoup import BeautifulSoup
import gd_composite

################################################################################

def run_admixture(ped_file, populations):
    prog = 'admixture'

    args = []
    args.append(prog)
    args.append(input_ped_file)
    args.append(populations)

    #print "args:", ' '.join(args)
    ofh = open('/dev/null', 'w')
    p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=ofh, stderr=sys.stderr)
    rc = p.wait()
    ofh.close()

def run_r(input_file, output_file, populations):
    prog = 'R'

    args = []
    args.append(prog)
    args.append('--vanilla')
    args.append('--quiet')
    args.append('--args')
    args.append(input_file)
    args.append(output_file)
    args.append(populations)

    _realpath = os.path.realpath(__file__)
    _script_dir = os.path.dirname(_realpath)
    r_script_file = os.path.join(_script_dir, 'population_structure.r')

    ifh = open(r_script_file)
    ofh = open('/dev/null', 'w')
    p = subprocess.Popen(args, bufsize=-1, stdin=ifh, stdout=ofh, stderr=None)
    rc = p.wait()
    ifh.close()
    ofh.close()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

def get_populations(input):
    pops = []
    pop_names = {}

    with open(input) as fh:
        soup = BeautifulSoup(fh)
        misc = soup.find('div', {'id': 'gd_misc'})

        return 'Populations\n{0}'.format(misc('ul')[0])

################################################################################

if len(sys.argv) != 6:
    print >> sys.stderr, "Usage"
    sys.exit(1)

input_html_file, input_ped_file, output_file, extra_files_path, populations = sys.argv[1:6]
populations_html = get_populations(input_html_file)

run_admixture(input_ped_file, populations)

ped_base = os.path.basename(input_ped_file)
if ped_base.endswith('.ped'):
    ped_base = ped_base[:-4]

p_file = '%s.%s.P' % (ped_base, populations)
q_file = '%s.%s.Q' % (ped_base, populations)

mkdir_p(extra_files_path)
numeric_output_file = os.path.join(extra_files_path, 'numeric.txt')
shutil.copy2(q_file, numeric_output_file)
os.remove(p_file)
os.remove(q_file)

graphical_output_file = os.path.join(extra_files_path, 'graphical.pdf')
run_r(numeric_output_file, graphical_output_file, populations)

################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('Population structure Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='graphical.pdf', value='graphical.pdf', display_type=display_file)
out_txt = gd_composite.Parameter(name='numeric.txt', value='numeric.txt', display_type=display_file)

info_page.add_output_parameter(out_pdf)
info_page.add_output_parameter(out_txt)

in_pops = gd_composite.Parameter(description='Number of populations', value=populations, display_type=display_value)

info_page.add_input_parameter(in_pops)

misc_pops = gd_composite.Parameter(description=populations_html, display_type=display_value)

info_page.add_misc(misc_pops)


with open (output_file, 'w') as ofh:
    print >> ofh, info_page.render()


sys.exit(0)
