#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Cluster_GOKEGG.py
#       
#       Copyright 2013 Oscar Reina <oscar@niska.bx.psu.edu>
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
import os
from networkx import connected_components, Graph, clustering
from numpy import percentile
from decimal import Decimal, getcontext
from itertools import permutations, combinations
import sys


def rtrnClustrsOnCltrCoff(dNodesWeightMin, threshold, perctile=True):
    """
    From a file with three columns: nodeA, nodeB and a score, it returns
    the strong and weak connected components produced when the edges
    """

    # ~
    Gmin = Graph()
    for nodeA, nodeB in dNodesWeightMin:
        wMin = dNodesWeightMin[nodeA, nodeB]
        Gmin.add_edge(nodeA, nodeB, weight=wMin)
    # ~
    clstrCoffcMin = clustering(Gmin, weight='weight')
    # ~
    if perctile:
        umbralMin = percentile(clstrCoffcMin.values(), threshold)
    else:
        umbralMin = threshold
    # ~
    GminNdsRmv = [x for x in clstrCoffcMin if clstrCoffcMin[x] < umbralMin]
    # ~
    Gmin.remove_nodes_from(GminNdsRmv)
    # ~
    dTermCmptNumbWkMin = rtrndata(Gmin)
    # ~
    salelClustr = []
    srtdterms = sorted(dTermCmptNumbWkMin.keys())
    for echTerm in srtdterms:
        try:
            MinT = dTermCmptNumbWkMin[echTerm]
        except:
            MinT = '-'
        salelClustr.append('\t'.join([echTerm, MinT]))
    # ~
    return salelClustr


def rtrndata(G):
    """
    return a list of terms and its clustering, as well as clusters from
    a networkx formatted file.
    """
    # ~
    cntCompnts = 0
    dTermCmptNumbWk = {}
    for echCompnt in connected_components(G):
        cntCompnts += 1
        # print '.'.join(echCompnt)
        for echTerm in echCompnt:
            dTermCmptNumbWk[echTerm] = str(cntCompnts)
    # ~
    return dTermCmptNumbWk


def rtrnCATcENSEMBLc(inCATfile, classClmns, ENSEMBLTcolmn, nonHdr=True):
    """
    return a dictionary of all the categories in an input file with
    a set of genes. Takes as input a file with categories an genes.
    """

    dCAT = {}
    dENSEMBLTCAT = {}

    for eachl in open(inCATfile, 'r'):
        if nonHdr and eachl.strip():
            ENSEMBLT = eachl.splitlines()[0].split('\t')[ENSEMBLTcolmn]
            sCAT = set()
            for CATcolmn in classClmns:
                sCAT.update(set([x for x in eachl.splitlines()[0].split('\t')[
                    CATcolmn].split('.')]))
            sCAT = sCAT.difference(set(['', 'U', 'N']))
            if len(sCAT) > 0:
                dENSEMBLTCAT[ENSEMBLT] = sCAT
            for CAT in sCAT:
                try:
                    dCAT[CAT].add(ENSEMBLT)
                except:
                    dCAT[CAT] = set([ENSEMBLT])
        nonHdr = True
    # ~
    dCAT = dict([(x, len(dCAT[x])) for x in dCAT.keys()])
    # ~
    return dCAT, dENSEMBLTCAT


def calcDistance(sCAT1, sCAT2):
    """
    takes as input two set of genesin different categories and returns
    a value 1-percentage of gene shared cat1->cat2, and cat2->cat1.
    """

    getcontext().prec = 5
    lgensS1 = Decimal(len(sCAT1))
    lgensS2 = Decimal(len(sCAT2))
    shrdGns = sCAT1.intersection(sCAT2)
    lenshrdGns = len(shrdGns)
    # ~
    dC1C2 = 1 - (lenshrdGns / lgensS1)
    dC2C1 = 1 - (lenshrdGns / lgensS2)
    # ~
    return dC1C2, dC2C1


def rtnPrwsdtncs(dCAT, dENSEMBLTCAT):
    """
    return a mcl formated pairwise distances from a list of categories
    """

    # ~
    getcontext().prec = 5
    dCATdst = {}
    lENSEMBL = dENSEMBLTCAT.keys()
    l = len(lENSEMBL)
    c = 0
    for ENSEMBL in lENSEMBL:
        c += 1
        lCAT = dENSEMBLTCAT.pop(ENSEMBL)
        for CAT1, CAT2 in combinations(lCAT, 2):
            try:
                dCATdst[CAT1, CAT2] += 1
            except:
                dCATdst[CAT1, CAT2] = 1
            try:
                dCATdst[CAT2, CAT1] += 1
            except:
                dCATdst[CAT2, CAT1] = 1
    # ~
    dNodesWeightMin = {}
    l = len(dCATdst)
    for CAT1, CAT2 in dCATdst.keys():
        shrdGns = dCATdst.pop((CAT1, CAT2))
        dC1C2 = float(shrdGns)
        nodeA, nodeB = sorted([CAT1, CAT2])
        try:
            cscor = dNodesWeightMin[nodeA, nodeB]
            if cscor >= dC1C2:
                dNodesWeightMin[nodeA, nodeB] = dC1C2
        except:
            dNodesWeightMin[nodeA, nodeB] = dC1C2
    #
    return dNodesWeightMin


def parse_class_columns(val, max_col):
    int_list = []

    for elem in [x.strip() for x in val.split(',')]:
        if elem[0].lower() != 'c':
            print >> sys.stderr, "bad column format:", elem
            sys.exit(1)

    int_val = as_int(elem[1:])

    if int_val is None:
        print >> sys.stderr, "bad column format:", elem
        sys.exit(1)
    elif not 1 <= int_val <= max_col:
        print >> sys.stderr, "column out of range:", elem
        sys.exit(1)

    int_list.append(int_val - 1)

    return int_list


def as_int(val):
    try:
        return int(val)
    except ValueError:
        return None
    else:
        raise


def main():
    """
    """

    # ~ bpython cluster_onConnctdComps.py --input=../conctFinal_CEU.tsv --outfile=../borrar.txt --threshold=90 --ENSEMBLTcolmn=1 --classClmns='20 22'
    parser = argparse.ArgumentParser(
        description='Returns the count of genes in ...')
    parser.add_argument('--input', metavar='input TXT file', type=str,
                        help='the input file with the table in txt format.',
                        required=True)
    parser.add_argument('--input_columns', metavar='input INT value', type=int,
                        help='the number of columns in the input file.',
                        required=True)
    parser.add_argument('--outfile', metavar='input TXT file', type=str,
                        help='the output file with the connected components.',
                        required=True)
    parser.add_argument('--threshold', metavar='input FLOAT value', type=float,
                        help='the threshold to disconnect the nodes.',
                        required=True)
    parser.add_argument('--ENSEMBLTcolmn', metavar='input INT file', type=int,
                        help='the column with the ENSEMBLE code in the input.',
                        required=True)
    parser.add_argument('--classClmns', metavar='input STR value', type=str,
                        help='the list of columns with the gene categories separated by space.',
                        required=True)
    args = parser.parse_args()
    infile = args.input
    threshold = args.threshold
    outfile = args.outfile
    ENSEMBLTcolmn = args.ENSEMBLTcolmn
    classClmns = parse_class_columns(args.classClmns, args.input_columns)
    # ~
    dCAT, dENSEMBLTCAT = rtrnCATcENSEMBLc(infile, classClmns, ENSEMBLTcolmn)
    dNodesWeightMin = rtnPrwsdtncs(dCAT, dENSEMBLTCAT)
    salelClustr = rtrnClustrsOnCltrCoff(dNodesWeightMin, threshold)
    # ~
    with open(outfile, 'w') as salef:
        print >> salef, '\n'.join(salelClustr)
    # ~
    # ~


if __name__ == '__main__':
    main()
