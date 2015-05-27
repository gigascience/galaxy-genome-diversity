#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       GOFisher.py
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
import sys
from fisher import pvalue as fisher
from decimal import Decimal, getcontext
from math import lgamma, exp, factorial


def binProb(SAPs_GO, NoSAPs_GO, SAPs_all, NoSAPs_all, CntGO_All, totalSAPs,
            pGO):
    """
    Returns binomial probability.
    """

    def comb(CntGO_All, k):
        return factorial(CntGO_All) / Decimal(
            str(factorial(k) * factorial(CntGO_All - k)))

    probLow = 0
    for k in range(0, SAPs_GO + 1):
        cp = Decimal(str(comb(CntGO_All, k)))
        bp = Decimal(str(pGO ** k))
        dp = Decimal(str(1.0 - pGO)) ** Decimal(str(CntGO_All - k))
        probLow += cp * bp * dp
    # ~
    probHigh = 0
    for k in range(int(SAPs_GO), CntGO_All + 1):
        cp = Decimal(str(comb(CntGO_All, k)))
        bp = Decimal(str(pGO ** k))
        dp = Decimal(str(1.0 - pGO)) ** Decimal(str(CntGO_All - k))
        probHigh += cp * bp * dp
    return probLow, probHigh


def gauss_hypergeom(X, CntGO_All, SAPs_all, totalSAPs):
    CntGO_All, SAPs_all, totalSAPs
    """
    Returns the probability of drawing X successes of SAPs_all marked items
    in CntGO_All draws from a bin of totalSAPs total items
    """

    def logchoose(ni, ki):
        try:
            lgn1 = lgamma(ni + 1)
            lgk1 = lgamma(ki + 1)
            lgnk1 = lgamma(ni - ki + 1)
        except ValueError:
            raise ValueError
        return lgn1 - (lgnk1 + lgk1)

    # ~
    r1 = logchoose(SAPs_all, X)
    try:
        r2 = logchoose(totalSAPs - SAPs_all, CntGO_All - X)
    except ValueError:
        return 0
    r3 = logchoose(totalSAPs, CntGO_All)
    return exp(r1 + r2 - r3)


def hypergeo_sf(SAPs_GO, NoSAPs_GO, SAPs_all, NoSAPs_all, CntGO_All, totalSAPs,
                pGO):
    """
    Runs Hypergeometric probability test
    """
    s = 0
    t = 0
    for i in range(SAPs_GO, min(SAPs_all, CntGO_All) + 1):
        s += max(gauss_hypergeom(i, CntGO_All, SAPs_all, totalSAPs), 0.0)
    for i in range(0, SAPs_GO + 1):
        t += max(gauss_hypergeom(i, CntGO_All, SAPs_all, totalSAPs), 0.0)
    return min(max(t, 0.0), 1), min(max(s, 0.0), 1)


def fisherexct(SAPs_GO, NoSAPs_GO, SAPs_all, NoSAPs_all, CntGO_All, totalSAPs,
               pGO):
    """
    Runs Fisher's exact test
    """
    ftest = fisher(SAPs_GO, NoSAPs_GO, SAPs_all, NoSAPs_all)
    probLow, probHigh = ftest.left_tail, ftest.right_tail
    return probLow, probHigh


def rtrnGOcENSEMBLc(inExtnddfile, columnENSEMBLTExtndd, columnGOExtndd):
    """
    """
    dGOTENSEMBLT = {}
    for eachl in open(inExtnddfile, 'r'):
        if eachl.strip():
            ENSEMBLT = eachl.splitlines()[0].split('\t')[columnENSEMBLTExtndd]
            GOTs = set(
                eachl.splitlines()[0].split('\t')[columnGOExtndd].split('.'))
            GOTs = GOTs.difference(set(['', 'U', 'N']))
            for GOT in GOTs:
                try:
                    dGOTENSEMBLT[GOT].add(ENSEMBLT)
                except:
                    dGOTENSEMBLT[GOT] = set([ENSEMBLT])
    ENSEMBLTGinGO = set.union(*dGOTENSEMBLT.values())
    return dGOTENSEMBLT, ENSEMBLTGinGO


def rtrnENSEMBLcSAPs(inSAPsfile, columnENSEMBLT, ENSEMBLTGinGO):
    """
    returns a set of the ENSEMBLT codes present in the input list and
    in the GO file
    """
    sENSEMBLTSAPsinGO = set()
    for eachl in open(inSAPsfile, 'r'):
        ENSEMBLT = eachl.splitlines()[0].split('\t')[columnENSEMBLT]
        if ENSEMBLT in ENSEMBLTGinGO:
            sENSEMBLTSAPsinGO.add(ENSEMBLT)
    return sENSEMBLTSAPsinGO


def rtrnCounts(dGOTENSEMBLT, ENSEMBLTGinGO, sENSEMBLTSAPsinGO, statsTest):
    """
    returns a list of the ENSEMBLT codes present in the input list and
    in the GO file. The terms in this list are: 'Go Term','# Genes in
    the GO Term','# Genes in the list and in the GO Term','Enrichement
    of the GO Term for genes in the input list','Genes in the input list
    present in the GO term'
    """
    totalSAPs = len(ENSEMBLTGinGO)
    SAPs_all = len(sENSEMBLTSAPsinGO)
    NoSAPs_all = totalSAPs - SAPs_all
    pGO = SAPs_all / float(totalSAPs)
    # ~
    lp = len(dGOTENSEMBLT)
    cnt = 0
    # ~
    if statsTest == 'fisher':
        ptest = fisherexct
    elif statsTest == 'hypergeometric':
        ptest = hypergeo_sf
    elif statsTest == 'binomial':
        ptest = binProb
    # ~
    ltfreqs = []
    for echGOT in dGOTENSEMBLT:
        cnt += 1
        CntGO_All = len(dGOTENSEMBLT[echGOT])
        SAPs_GO = len(dGOTENSEMBLT[echGOT].intersection(sENSEMBLTSAPsinGO))
        NoSAPs_GO = CntGO_All - SAPs_GO
        probLow, probHigh = ptest(SAPs_GO, NoSAPs_GO, SAPs_all, NoSAPs_all,
                                  CntGO_All, totalSAPs, pGO)
        ltfreqs.append(
            [(SAPs_GO / Decimal(CntGO_All)), SAPs_GO, probLow, probHigh,
             echGOT])
    # ~
    ltfreqs.sort()
    ltfreqs.reverse()
    outl = []
    cper, crank = Decimal('2'), 0
    # ~
    getcontext().prec = 2  # set 2 decimal places
    for perc, cnt_go, pvalLow, pvalHigh, goTerm in ltfreqs:
        if perc < cper:
            crank += 1
            cper = perc
        outl.append('\t'.join(
            [str(cnt_go), str(Decimal(perc) * Decimal('1.0')), str(crank),
             str(Decimal(pvalLow) * Decimal('1.0')),
             str(Decimal(pvalHigh) * Decimal('1.0')), goTerm]))
    # ~
    return outl


def main():
    # ~
    parser = argparse.ArgumentParser(
        description='Returns the count of genes in GO categories and their statistical overrrepresentation, from a list of genes and an extended file (i.e. plane text with ENSEMBLT and GO terms).')
    parser.add_argument('--input', metavar='input TXT file', type=str,
                        help='the input file with the table in txt format.',
                        required=True)
    parser.add_argument('--inExtnddfile', metavar='input TXT file', type=str,
                        help='the input file with the extended table in txt format.',
                        required=True)
    parser.add_argument('--output', metavar='output TXT file', type=str,
                        help='the output file with the table in txt format.',
                        required=True)
    parser.add_argument('--columnENSEMBLT', metavar='column number', type=int,
                        help='column with the ENSEMBL transcript code in the input file.',
                        required=True)
    parser.add_argument('--columnENSEMBLTExtndd', metavar='column number',
                        type=int,
                        help='column with the ENSEMBL transcript code in the extended file.',
                        required=True)
    parser.add_argument('--columnGOExtndd', metavar='column number', type=int,
                        help='column with the GO terms in the extended file.',
                        required=True)
    parser.add_argument('--statsTest', metavar='input TXT file', type=str,
                        help='statistical test to compare GO terms (i.e. fisher, hypergeometric, binomial).',
                        required=True)

    args = parser.parse_args()

    inSAPsfile = args.input
    inExtnddfile = args.inExtnddfile
    saleGOPCount = args.output
    columnENSEMBLT = args.columnENSEMBLT
    columnENSEMBLTExtndd = args.columnENSEMBLTExtndd
    columnGOExtndd = args.columnGOExtndd
    statsTest = args.statsTest

    # ~
    dGOTENSEMBLT, ENSEMBLTGinGO = rtrnGOcENSEMBLc(inExtnddfile,
                                                  columnENSEMBLTExtndd,
                                                  columnGOExtndd)
    sENSEMBLTSAPsinGO = rtrnENSEMBLcSAPs(inSAPsfile, columnENSEMBLT,
                                         ENSEMBLTGinGO)
    outl = rtrnCounts(dGOTENSEMBLT, ENSEMBLTGinGO, sENSEMBLTSAPsinGO, statsTest)
    # ~
    saleGOPCount = open(saleGOPCount, 'w')
    saleGOPCount.write('\n'.join(outl))
    saleGOPCount.close()
    # ~
    return 0


if __name__ == '__main__':
    main()
