#!/usr/bin/env python

import os
import sys
import math

import matplotlib as mpl
mpl.use('PDF')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

################################################################################

def build_chrom_len_dict(dbkey, galaxy_data_index_dir):
    chrom_len_root = os.path.join(galaxy_data_index_dir, 'shared/ucsc/chrom')
    chrom_len_file = '{0}.len'.format(dbkey)
    chrom_len_path = os.path.join(chrom_len_root, chrom_len_file)

    chrom_len = {}

    try:
        with open(chrom_len_path) as fh:
            for line in fh:
                line = line.rstrip('\r\n')
                elems = line.split()
                if len(elems) == 2:
                    chrom = elems[0]
                    length = int(elems[1])
                    chrom_len[chrom] = length
    except:
        pass

    return chrom_len

def parse_input_file(input_file):
    chroms = []
    individuals = []
    data = {}
    chrom_len = {}
    used_states = []

    with open(input_file) as fh:
        for line in fh:
            line = line.strip()
            if line:
                elems = line.split()
                chrom = elems[0]
                p1, p2, state = map(int, elems[1:4])
                id = elems[4]

                if state not in used_states:
                    used_states.append(state)

                if chrom not in chroms:
                    chroms.append(chrom)

                if id not in individuals:
                    individuals.append(id)

                data.setdefault(chrom, {})
                data[chrom].setdefault(id, [])
                data[chrom][id].append((p1, p2, state))

                if p2 > chrom_len.setdefault(chrom, 0):
                    chrom_len[chrom] = p2

    return chroms, individuals, data, chrom_len, used_states

def check_chroms(chroms, chrom_len, dbkey):
    error = 0
    for chrom in chroms:
        if chrom not in chrom_len:
            print >> sys.stderr, "Can't find length for {0} chromosome {1}".format(dbkey, chrom)
            error = 1
    if error:
        sys.exit(1)

def check_data(data, chrom_len, dbkey):
    error = 0
    for chrom in data:
        chrom_beg = 0
        chrom_end = chrom_len[chrom]
        for individual in data[chrom]:
            for p1, p2, state in data[chrom][individual]:
                if p1 >= p2:
                    print >> sys.stderr, "Bad data line: begin >= end: {0} {1} {2} {3}".format(chrom, p1, p2, state, individual)
                    error = 1
                if p1 < chrom_beg or p2 > chrom_end:
                    print >> sys.stderr, "Bad data line: outside {0} boundaries[{1} - {2}]: {3} {4} {5} {6}".format(dbkey, chrom_beg, chrom_end, chrom, p1, p2, state, individual)
                    error = 1
    if error:
        sys.exit(1)

def make_rectangle(p1, p2, color, bottom=0.0, top=1.0):
    verts = [
        (p1, bottom),   # left, bottom
        (p1, top),      # left, top
        (p2, top),      # right, top
        (p2, bottom),   # right, bottom
        (0.0, 0.0)      # ignored
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY
    ]

    path = Path(verts, codes)
    return patches.PathPatch(path, facecolor=color, lw=0)

def make_split_rectangle(p1, p2, top_color, bottom_color):
    patch1 = make_rectangle(p1, p2, bottom_color, top=0.5)
    patch2 = make_rectangle(p1, p2, top_color, bottom=0.5)
    return [patch1, patch2]

def make_state_rectangle_2pop(p1, p2, state, chrom, individual):
    p1_color = 'r'
    p2_color = 'g'
    heterochromatin_color = '#c7c7c7'

    if state == 0:
        return [ make_rectangle(p1, p2, heterochromatin_color) ]
    elif state == 1:
        return [ make_rectangle(p1, p2, p1_color) ]
    elif state == 2:
        return [ make_rectangle(p1, p2, p2_color) ]
    elif state == 3:
        return make_split_rectangle(p1, p2, p1_color, p2_color)
    else:
        print >> sys.stderr, "Unknown state: {0}: {1} {2} {3} {4}".format(state, chrom, p1, p2, state, individual)
        sys.exit(1)

def make_state_rectangle_3pop(p1, p2, state, chrom, individual):
    p1_color = 'r'
    p2_color = 'g'
    p3_color = 'b'
    heterochromatin_color = '#c7c7c7'

    if state == 0:
        return [ make_rectangle(p1, p2, heterochromatin_color) ]
    if state == 1:
        return [ make_rectangle(p1, p2, p1_color) ]
    if state == 2:
        return [ make_rectangle(p1, p2, p2_color) ]
    if state == 3:
        return [ make_rectangle(p1, p2, p3_color) ]
    if state == 4:
        return make_split_rectangle(p1, p2, p1_color, p2_color)
    if state == 5:
        return make_split_rectangle(p1, p2, p1_color, p3_color)
    if state == 6:
        return make_split_rectangle(p1, p2, p2_color, p3_color)
    else:
        print >> sys.stderr, "Unknown state: {0}: {1} {2} {3} {4}".format(state, chrom, p1, p2, state, individual)
        sys.exit(1)

def nicenum(num, round=False):
    if num == 0:
        return 0.0

    exp = int(math.floor(math.log10(num)))
    f = num / math.pow(10, exp)

    if round:
        if f < 1.5:
            nf = 1.0
        elif f < 3.0:
            nf = 2.0
        elif f < 7.0:
            nf = 5.0
        else:
            nf = 10.0
    else:
        if f <= 1.0:
            nf = 1.0
        elif f <= 2.0:
            nf = 2.0
        elif f <= 5.0:
            nf = 5.0
        else:
            nf = 10.0

    return nf * pow(10, exp)

def tick_foo(beg, end, loose=False):
    ntick = 10

    range = nicenum(end - beg, round=False)
    d = nicenum(range/(ntick - 1), round=True)
    digits = int(math.floor(math.log10(d)))

    if loose:
        graph_min = math.floor(beg/d) * d
        graph_max = math.ceil(end/d) * d
    else:
        graph_min = beg
        graph_max = end

    nfrac = max([-1 * digits, 0])
    vals = []

    stop = graph_max
    if loose:
        stop = graph_max + (0.5 * d)

    x = graph_min
    while x <= stop:
        vals.append(int(x))
        x += d

    vals = vals[1:]

#    if not loose:
#        if vals[-1] < graph_max:
#            vals.append(int(graph_max))

    labels = []
    for val in vals:
        labels.append('{0}'.format(int(val/math.pow(10, digits))))

#   labels.append('{0:.1f}'.format(vals[-1]/math.pow(10, digits)))

    return vals, labels

################################################################################
################################################################################
################################################################################
################################################################################

def space_for_legend(plot_params):
    space = 0.0

    legend_states = plot_params['legend_states']
    if legend_states:
        ind_space = plot_params['ind_space']
        ind_height = plot_params['ind_height']
        space += len(legend_states) * (ind_space + ind_height) - ind_space

    return space

################################################################################

def space_for_chroms(plot_params, chroms, individuals, data):
    space_dict = {}

    chrom_height = plot_params['chrom_height']
    ind_space = plot_params['ind_space']
    ind_height = plot_params['ind_height']

    for chrom in chroms:
        space_dict[chrom] = chrom_height

        individual_count = 0
        for individual in individuals:
            if individual in data[chrom]:
                individual_count += 1

        space_dict[chrom] += individual_count * (ind_space + ind_height)

    return space_dict

################################################################################

def make_dpmix_plot(input_dbkey, input_file, output_file, galaxy_data_index_dir, state2name=None, populations=3):
    fs_chrom_len = build_chrom_len_dict(input_dbkey, galaxy_data_index_dir)
    chroms, individuals, data, chrom_len, used_states = parse_input_file(input_file)

    ## populate chrom_len
    for chrom in chrom_len.keys():
        if chrom in fs_chrom_len:
            chrom_len[chrom] = fs_chrom_len[chrom]

    #check_chroms(chroms, chrom_len, input_dbkey)
    check_data(data, chrom_len, input_dbkey)

    ## plot parameters
    plot_params = {
        'plot_dpi':        300,
        'page_width':     8.50,
        'page_height':   11.00,
        'top_margin':     0.10,
        'bottom_margin':  0.10,
        'chrom_space':    0.25,
        'chrom_height':   0.25,
        'ind_space':      0.10,
        'ind_height':     0.25,
        'legend_space':   0.10
    }

    ## in the legend, only print out states that are
    ##   1) in the data
    ##    - AND -
    ##   2) in the state2name map
    legend_states = []
    if state2name is not None:
        for state in used_states:
            if state in state2name:
                legend_states.append(state)

    plot_params['legend_states'] = legend_states

    ## choose the correct make_state_rectangle method
    if populations == 3:
        plot_params['rectangle_method'] = make_state_rectangle_3pop
    elif populations == 2:
        plot_params['rectangle_method'] = make_state_rectangle_2pop

    pdf_pages = PdfPages(output_file)

	## generate a list of chroms for each page

    needed_for_legend = space_for_legend(plot_params)
    needed_for_chroms = space_for_chroms(plot_params, chroms, individuals, data)

    chrom_space_per_page = plot_params['page_height']
    chrom_space_per_page -= plot_params['top_margin'] + plot_params['bottom_margin']
    chrom_space_per_page -= needed_for_legend + plot_params['legend_space']
    chrom_space_per_page -= plot_params['chrom_space']

    chroms_left = chroms[:]
    pages = []

    space_left = chrom_space_per_page
    chrom_list = []

    while chroms_left:
        chrom = chroms_left.pop(0)
        space_needed = needed_for_chroms[chrom] + plot_params['chrom_space']
        if (space_needed > chrom_space_per_page):
            print >> sys.stderr, 'Multipage chroms not yet supported'
            sys.exit(1)

		## sometimes 1.9 - 1.9 < 0 (-4.4408920985e-16)
		## so, we make sure it's not more than a millimeter over
        if space_left - space_needed > -0.04:
            chrom_list.append(chrom)
            space_left -= space_needed
        else:
            pages.append(chrom_list[:])
            chrom_list = []
            chroms_left.insert(0, chrom)
            space_left = chrom_space_per_page

    ############################################################################

    plot_dpi = plot_params['plot_dpi']
    page_width = plot_params['page_width']
    page_height = plot_params['page_height']
    top_margin = plot_params['top_margin']
    ind_space = plot_params['ind_space']
    ind_height = plot_params['ind_height']
    make_state_rectangle = plot_params['rectangle_method']
    legend_space = plot_params['legend_space']
    chrom_space = plot_params['chrom_space']
    chrom_height = plot_params['chrom_height']

    for page in pages:
        fig = plt.figure(figsize=(page_width, page_height), dpi=plot_dpi)
        bottom = 1.0 - (top_margin/page_height)

        # print legend
        if legend_states:
            top = True
            for state in sorted(legend_states):
                if top:
                    bottom -= ind_height/page_height
                    top = False
                else:
                    bottom -= (ind_space + ind_height)/page_height

                ax1 = fig.add_axes([0.0, bottom, 0.09, ind_height/page_height])
                plt.axis('off')
                ax1.set_xlim(0, 1)
                ax1.set_ylim(0, 1)
                for patch in make_state_rectangle(0, 1, state, 'legend', state2name[state]):
                    ax1.add_patch(patch)

                ax2 = fig.add_axes([0.10, bottom, 0.88, ind_height/page_height], frame_on=False)
                plt.axis('off')
                plt.text(0.0, 0.5, state2name[state], fontsize=10, ha='left', va='center')

            bottom -= legend_space/page_height

        # print chroms
        top = True
        for chrom in page:
            length = chrom_len[chrom]
            vals, labels = tick_foo(0, length)

            if top:
                bottom -= chrom_height/page_height
                top = False
            else:
                bottom -= (chrom_space + chrom_height)/page_height

            ax = fig.add_axes([0.0, bottom, 1.0, chrom_height/page_height])
            plt.axis('off')
            plt.text(0.5, 0.5, chrom, fontsize=14, ha='center')

            individual_count = 0
            for individual in individuals:
                if individual in data[chrom]:
                    individual_count += 1

            i = 0
            for individual in individuals:
                if individual in data[chrom]:
                    i += 1
                    bottom -= (ind_space + ind_height)/page_height

                    ax1 = fig.add_axes([0.0, bottom, 0.09, ind_height/page_height])
                    plt.axis('off')
                    plt.text(1.0, 0.5, individual, fontsize=10, ha='right', va='center')

                    ax2 = fig.add_axes([0.10, bottom, 0.88, ind_height/page_height], frame_on=False)
                    ax2.set_xlim(0, length)
                    ax2.set_ylim(0, 1)

                    if i != individual_count:
                        plt.axis('off')
                    else:
                        ax2.tick_params(top=False, left=False, right=False, labelleft=False)
                        ax2.set_xticks(vals)
                        ax2.set_xticklabels(labels)

                    for p1, p2, state in sorted(data[chrom][individual]):
                        for patch in make_state_rectangle(p1, p2, state, chrom, individual):
                            ax2.add_patch(patch)

                    # extend last state to end of chrom
                    if p2 < length:
                        for patch in make_state_rectangle(p2, length, state, chrom, individual):
                            ax2.add_patch(patch)


        pdf_pages.savefig(fig)
        plt.close(fig)

    pdf_pages.close()

################################################################################

if __name__ == '__main__':
    make_dpmix_plot('loxAfr3', 'output.dat', 'output2_files/picture.pdf', '/scratch/galaxy/home/oocyte/galaxy_oocyte/tool-data', state2name={0: 'heterochromatin', 1: 'reference', 2: 'asian'}, populations=2)
#    input_dbkey, input_file, output_file, galaxy_data_index_dir = sys.argv[1:5]
#    make_dpmix_plot(input_dbkey, input_file, output_file, galaxy_data_index_dir)
    sys.exit(0)

## notes
# 1) pass in a state to name mapping
# 2) only print out names for states which exist in the data, and are in the state to name mapping
