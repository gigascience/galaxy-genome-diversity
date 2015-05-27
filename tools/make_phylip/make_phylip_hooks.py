#!/usr/bin/env python

def exec_before_job(app, inp_data=None, out_data=None, tool=None, param_dict=None):
    pass

def exec_after_process(app, inp_data=None, out_data=None, param_dict=None, tool=None, stdout=None, stderr=None):
    output_name = 'output1'

    ## check for output
    try:
        first_output = out_data[output_name]
    except:
        return

    ## check for collected datasets
    try:
        collected_dict = param_dict['__collected_datasets__']['primary'][output_name]
    except:
        return

    if len(collected_dict.keys()) == 0:
        return

    ## check for fasta file
    try:
        fasta_file = inp_data['fasta_input']
    except:
        return

    ## find missing fasta header
    first_output_name = None
    with open(fasta_file.get_file_name()) as fh:
        for line in fh:
            if line[0] != '>':
                continue
            name = line[1:]
            name = name.strip()
            name = name.split()[0]
            name = name.replace('_', '-')
            if name not in collected_dict:
                first_output_name = name
                break

    ## fix name
    if first_output_name is not None:
        first_output.name = '%s (%s)' % (first_output.name, first_output_name)

