#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

def parse_string(str_arg, ind_token2col):
    columns = []

    string = gd_util.unwrap_string(str_arg)
    tokens = find_tokens(string, ind_token2col)

    for token in tokens:
        col = ind_token2col[token]
        if col not in columns:
            columns.append(col)

    return columns

def find_tokens(string, tokens):
    rv = []
    for token in tokens:
        if token in string:
            if token not in rv:
                rv.append(token)
    return rv

################################################################################

if len(sys.argv) != 6:
    gd_util.die('Usage')

input, output, ind_arg, cb_arg, str_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p_cb = Population()
p_cb.from_wrapped_dict(cb_arg)

if not p_total.is_superset(p_cb):
    gd_util.die('There is a checked individual that does not appear in the SNP table')

################################################################################

ind_col2name = {}
ind_token2col = {}
for col in p_total.column_list():
    individual = p_total.individual_with_column(col)
    name = individual.name
    ind_col2name[col] = name
    first_token = name.split()[0]
    if first_token not in ind_token2col:
        ind_token2col[first_token] = col
    else:
        gd_util.die('duplicate first token: {0}'.format(first_token))

out_cols = p_cb.column_list()
str_cols = parse_string(str_arg, ind_token2col)

with open(output, 'w') as fh:
    for col in sorted(ind_col2name.keys()):
        if col in out_cols or col in str_cols:
            print >> fh, '\t'.join([str(x) for x in [col, ind_col2name[col], '']])

sys.exit(0)

