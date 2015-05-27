#!/usr/bin/env python

import base64
import json
import math
import re
import sys

identifier_regex = re.compile('[0-9A-Z_a-z]+$')

def unwrap_column_names(string):
    column_names = []
    string = unwrap_string(string)
    for line in string.split('\n'):
        line = line.strip()
        if is_identifier(line):
            column_names.append(line)
        else:
            die('invalid column name format: {}'.format(line))
    return column_names

def unwrap_string(string):
    try:
        decoded = base64.b64decode(string)
    except:
        die('invalid base64 string: {}'.format(string))
    return decoded

def is_identifier(string):
    match = identifier_regex.match(string)
    if match:
        return True
    else:
        return False

def read_individual_names(filename):
    tokens = []
    names = []
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            elems = line.split()

            columns = len(elems)
            if columns == 0:
                continue

            first_token = elems[0]

            if columns == 1:
                name = first_token
            else:
                keywords = ' '.join(elems[1:])
                name = ' '.join([first_token, keywords])

            if first_token not in tokens:
                tokens.append(first_token)
                names.append(name)
            else:
                die('duplicate first column entry in Names dataset: {}'.format(first_token))
    return names

def fold_line(line, maxlen, prefix):
    prefix_len = len(prefix)

    lines = []

    while len(line) > maxlen:
        split_points = []
        state = 0
        for i in range(maxlen - prefix_len):
            c = line[i]
            if state == 0:
                if c == '"':
                    state = 1
                elif c in [ '{', ':', ',', '}', '[', ']' ]:
                    split_points.append(i)
            elif state == 1:
                if c == '"':
                    state = 0
                elif c == '\\':
                    state = 2
            elif state == 2:
                state = 1
        idx = split_points[-1]
        lines.append('{0}{1}'.format(prefix, line[:idx+1]))
        line = line[idx+1:]

    lines.append('{0}{1}'.format(prefix, line))

    return lines

def die(message):
    print >> sys.stderr, message
    sys.exit(1)

################################################################################

type_to_columns = {
    'gd_snp':4,
    'gd_genotype':1
}

if len(sys.argv) != 12:
    print >> sys.stderr, 'Usage'
    sys.exit(1)

input, scaffold_col, pos_col, ref_col, rPos_col, preamble_arg, names, species_arg, dbkey, output_type, output = sys.argv[1:12]

preamble_column_names = unwrap_column_names(preamble_arg)
first_individual_column = len(preamble_column_names) + 1

individual_names = read_individual_names(names)

species = unwrap_string(species_arg)
if not is_identifier(species):
    die('invalid species format: {}'.format(species))

if not output_type in type_to_columns:
    die('unknown output type: {}'.format(output_type))
columns_per_individual = type_to_columns[output_type]

jdict = {}

column_names = preamble_column_names[:]
for i in range(1, len(individual_names) + 1):
    if output_type == 'gd_snp':
        column_names.append('{}A'.format(i))
        column_names.append('{}B'.format(i))
        column_names.append('{}G'.format(i))
        column_names.append('{}Q'.format(i))
    elif output_type == 'gd_genotype':
        column_names.append('{}G'.format(i))
    else:
        die('unknown output type: {}'.format(output_type))

jdict['column_names'] = column_names

individuals = []

for pos, individual in enumerate(individual_names):
    col = first_individual_column + pos * columns_per_individual
    individuals.append([individual, col])

jdict['individuals'] = individuals

jdict['scaffold'] = int(scaffold_col)
jdict['pos'] = int(pos_col)
jdict['ref'] = int(ref_col)
jdict['rPos'] = int(rPos_col)

jdict['species'] = species
jdict['dbkey'] = dbkey

json_string = json.dumps(jdict, separators=(',',':'), sort_keys=True)

min_cols = len(column_names)
pos_col = int(pos_col) - 1
rPos_col = int(rPos_col) - 1

def is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

with open(output, 'w') as ofh:
    lines = fold_line(json_string, 200, '#')
    for line in lines:
        print >> ofh, line

    with open(input) as fh:
        line_number = 0
        for line in fh:
            line_number += 1
            if line[0] == '#':
                continue
            line = line.rstrip('\r\n')
            elems = line.split('\t')
            if len(elems) < min_cols:
                die('Too few columns on line {0} of input file.  Expecting {1}, saw {2}.'.format(line_number, min_cols, len(elems)))
            if not is_int(elems[pos_col]):
                die('bad pos on line {0} column {1} of input file: {2}'.format(line_number, pos_col+1, elems[pos_col]))
            if not is_int(elems[rPos_col]):
                die('bad rPos on line {0} column {1} of input file: {2}'.format(line_number, rPos_col+1, elems[rPos_col]))
            print >> ofh, line
