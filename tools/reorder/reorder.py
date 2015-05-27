#!/usr/bin/env python

import sys

def parse_rangelist(string):
    rv = []

    tokens = strip_split(string, ',')
    for token in tokens:
        int_list = parse_token(token)
        for int_val in int_list:
            int_val -= 1
            if int_val not in rv:
                rv.append(int_val)

    return rv

def parse_token(token):
    values = strip_split(token, '-')
    num_values = len(values)

    if num_values not in [1, 2]:
        print >> sys.stderr, 'Error: "%s" is not a valid range' % token
        sys.exit(1)

    int_list = []
    for value in values:
        if value:
            int_val = as_int(value)

            if int_val < 1:
                print >> sys.stderr, 'Error: "%s" is not >= 1' % value
                sys.exit(1)

            int_list.append(int_val)
        else:
            print >> sys.stderr, 'Error: "%s" is not a valid range' % token
            sys.exit(1)

    if num_values == 1:
        return int_list

    a, b = int_list

    if a <= b:
        return range(a, b+1)
    else:
        return range(a, b-1, -1)

def strip_split(string, delim):
    return [elem.strip() for elem in string.split(delim)]

def as_int(string):
    try:
        val = int(string)
    except:
        print >> sys.stderr, 'Error: "%s" does not appear to be an integer' % string
        sys.exit(1)
    return val

def get_lines(filename):
    rv = []

    fh = open(filename)
    for line in fh:
        line = line.rstrip('\r\n')
        rv.append(line)
    fh.close()

    return rv

def reorder(old_lines, new_order, filename):
    max_index = len(old_lines) - 1

    fh = open(filename, 'w')

    for index in new_order:
        if index <= max_index:
            print >> fh, old_lines[index]
            old_lines[index] = None

    for line in old_lines:
        if line is not None:
            print >> fh, line

    fh.close()

if len(sys.argv) != 4:
    print >> sys.stderr, "Usage"
    sys.exit(1)

input, output, order_string = sys.argv[1:]

new_order = parse_rangelist(order_string)
old_lines = get_lines(input)
reorder(old_lines, new_order, output)

sys.exit(0)
