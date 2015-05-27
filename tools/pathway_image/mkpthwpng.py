#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       mkpthwpng.py
#       
#       Copyright 2011 Oscar Bedoya-Reina <oscar@niska.bx.psu.edu>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

import argparse
import mechanize
import os
import sys

# this return an image made up from a list of genes and pathway code
def rtnHTMLformat(tmpddGenrcgenPresent, sppPrefx, pthwcod, ouPthwpng):
    inpx = '\n'.join(tmpddGenrcgenPresent)  # inpx="ALDH2 color \nALDH3A1	color"
    request = mechanize.Request(
        "http://www.genome.jp/kegg/tool/map_pathway2.html")
    response = mechanize.urlopen(request)
    forms = mechanize.ParseResponse(response, backwards_compat=False)
    form = forms[0]
    form["unclassified"] = inpx
    form["org"] = sppPrefx
    request2 = form.click()
    response2 = mechanize.urlopen(request2)
    a = str(response2.read()).split('href="/kegg-bin/show_pathway?')[1]
    code = a.split('/')[0]  # response2.read()
    request = mechanize.Request(
        "http://www.genome.jp/kegg-bin/show_pathway?%s/%s.args" % (code, pthwcod))  # request=mechanize.Request("http://www.genome.jp/kegg-bin/show_pathway?%s/%s.args"%('13171478854246','hsa00410'))
    response = mechanize.urlopen(request)
    forms = mechanize.ParseResponse(response, backwards_compat=False)
    form = forms[1]
    status = ' NOT '
    try:
        imgf = str(forms[1]).split('/mark_pathway')[1].split('/')[0]
        os.system("wget --quiet http://www.genome.jp/tmp/mark_pathway%s/%s.png -O %s" % (imgf, pthwcod, ouPthwpng))
        status = ' '
    except:
        pass
    return 'A pathway image was%ssuccefully produced...' % status


def main():
    parser = argparse.ArgumentParser(description='Obtain KEGG images from a list of genes.')
    parser.add_argument('--input', metavar='input TXT file', type=str,
                        help='the input file with the table in txt format')
    parser.add_argument('--output', metavar='output PNG image', type=str,
                        help='the output image file in png format')
    parser.add_argument('--KEGGpath',
                        metavar='KEGG pathway code (i.e. cfa00230)', type=str,
                        help='the code of the pathway of interest')
    parser.add_argument('--posKEGGclmn', metavar='column number', type=int,
                        help='the column with the KEGG pathway code/name')
    parser.add_argument('--KEGGgeneposcolmn', metavar='column number', type=int,
                        help='column with the KEGG gene code')

    # ~Open arguments
    class C(object):
        pass

    fulargs = C()
    parser.parse_args(sys.argv[1:], namespace=fulargs)
    # test input vars
    inputf, outputf, KEGGpathw, posKEGGclmn, Kgeneposcolmn = fulargs.input, fulargs.output, fulargs.KEGGpath, fulargs.posKEGGclmn, fulargs.KEGGgeneposcolmn
    # make posKEGGclmn, Kgeneposcolmn 0-based
    sppPrefx = KEGGpathw[:3]
    posKEGGclmn -= 1
    Kgeneposcolmn -= 1
    # make a dictionary of valid genes
    dKEGGcPthws = dict([(x.split('\t')[Kgeneposcolmn], set(
        [y.split('=')[0] for y in x.split('\t')[posKEGGclmn].split('.')])) for x
                        in open(inputf).read().splitlines()[1:] if x.strip()])
    for mt1gene in [x for x in dKEGGcPthws.keys() if x.find('.') > -1]:  # to correct names with more than one gene
        pthwsAssotd = dKEGGcPthws.pop(mt1gene)
        for eachg in mt1gene.split('.'):
            dKEGGcPthws[eachg] = pthwsAssotd
    tmpddGenrcgenPresent = set()
    sKEGGc = dKEGGcPthws.keys()
    lsKEGGc = len(sKEGGc)
    ctPthw = 0
    while ctPthw < lsKEGGc:  # to save memory
        eachK = sKEGGc.pop()
        alPthws = dKEGGcPthws[eachK]
        if KEGGpathw in alPthws:
            tmpddGenrcgenPresent.add('\t'.join([eachK, 'red']))
        ctPthw += 1
    # run the program
    rtnHTMLformat(tmpddGenrcgenPresent, sppPrefx, KEGGpathw, outputf)
    return 0


if __name__ == '__main__':
    main()
