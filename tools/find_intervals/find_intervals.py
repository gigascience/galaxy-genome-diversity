#!/usr/bin/env python

import errno
import os
import subprocess
import sys

################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

def run_program(prog, args, stdout_file=None):
    #print "args:", ' '.join(args)
    p = subprocess.Popen(args, bufsize=-1, executable=prog, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate()
    rc = p.returncode

    if stdout_file is not None:
        with open(stdout_file, 'w') as ofh:
            print >> ofh, stdoutdata.rstrip('\r\n')

    if rc != 0:
        print >> sys.stderr, "FAILED: rc={0}: {1}".format(rc, ' '.join(args))
        print >> sys.stderr, stderrdata
        sys.exit(1)

################################################################################

if len(sys.argv) != 11:
    print "usage"
    sys.exit(1)

input, dbkey, output, output_files_path, chrom_col, pos_col, score_col, shuffles, cutoff, report_snps = sys.argv[1:11]

prog = 'sweep'

args = [ prog ]
args.append(input)
args.append(chrom_col)
args.append(pos_col)
args.append(score_col)
args.append(cutoff)
args.append(shuffles)
args.append(report_snps)

run_program(None, args, stdout_file=output)

if report_snps == "0":
    sys.exit(0)

################################################################################

mkdir_p(output_files_path)

bedgraph_filename = 'bedgraph.txt'
links_filename = os.path.join(output_files_path, 'links.txt')

data = []
links_data = []
    
with open(output) as fh:
    chrom = None
    for line in fh:
        line = line.rstrip('\r\n')
        if not line:
            continue
        if line[0] != ' ':
            # chrom line, add a link
            chrom, interval_begin, interval_end, interval_value = line.split('\t')
            links_data.append((chrom, int(interval_begin), int(interval_end)))
        else:
            # data line, add a bedgraph line
            begin, value = line.split()
            data.append((chrom, int(begin), value))

with open(bedgraph_filename, 'w') as ofh:
    print >> ofh, 'track type=bedGraph'
    for chrom, begin, value in sorted(data):
        print >> ofh, chrom, begin, begin+1, value

with open(links_filename, 'w') as ofh:
    for chrom, begin, end in sorted(links_data):
        print >> ofh, chrom, begin, end

################################################################################

chrom_sizes_filename = '{0}.chrom.sizes'.format(dbkey)

prog = 'fetchChromSizes'

args = [ prog ]
args.append(dbkey)

run_program(None, args, stdout_file=chrom_sizes_filename)

################################################################################

prog = 'bedGraphToBigWig'

args = [ prog ]
args.append(bedgraph_filename)
args.append(chrom_sizes_filename)
args.append(output)

run_program(None, args)

################################################################################

sys.exit(0)

