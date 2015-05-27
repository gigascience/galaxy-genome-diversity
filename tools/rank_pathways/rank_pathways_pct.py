#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       KEGGFisher.py
#       
#       Copyright 2013 Oscar Reina <oscar@niska.bx.psu.edu>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the pathways of the GNU General Public License as published by
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


def binProb(SAPs_KEGG, NoSAPs_KEGG, SAPs_all, NoSAPs_all, CntKEGG_All,
            totalSAPs, pKEGG):
    """
    Returns binomial probability.
    """

    def comb(CntKEGG_All, k):
        return factorial(CntKEGG_All) / Decimal(
            str(factorial(k) * factorial(CntKEGG_All - k)))

    probLow = 0
    for k in range(0, SAPs_KEGG + 1):
        cp = Decimal(str(comb(CntKEGG_All, k)))
        bp = Decimal(str(pKEGG ** k))
        dp = Decimal(str(1.0 - pKEGG)) ** Decimal(str(CntKEGG_All - k))
        probLow += cp * bp * dp
    # ~
    probHigh = 0
    for k in range(int(SAPs_KEGG), CntKEGG_All + 1):
        cp = Decimal(str(comb(CntKEGG_All, k)))
        bp = Decimal(str(pKEGG ** k))
        dp = Decimal(str(1.0 - pKEGG)) ** Decimal(str(CntKEGG_All - k))
        probHigh += cp * bp * dp
    return probLow, probHigh


def gauss_hypergeom(X, CntKEGG_All, SAPs_all, totalSAPs):
    CntKEGG_All, SAPs_all, totalSAPs
    """
    Returns the probability of drawing X successes of SAPs_all marked items
    in CntKEGG_All draws from a bin of totalSAPs total items
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
        r2 = logchoose(totalSAPs - SAPs_all, CntKEGG_All - X)
    except ValueError:
        return 0
    r3 = logchoose(totalSAPs, CntKEGG_All)
    return exp(r1 + r2 - r3)


def hypergeo_sf(SAPs_KEGG, NoSAPs_KEGG, SAPs_all, NoSAPs_all, CntKEGG_All,
                totalSAPs, pKEGG):
    """
    Runs Hypergeometric probability test
    """
    s = 0
    t = 0
    for i in range(SAPs_KEGG, min(SAPs_all, CntKEGG_All) + 1):
        s += max(gauss_hypergeom(i, CntKEGG_All, SAPs_all, totalSAPs), 0.0)
    for i in range(0, SAPs_KEGG + 1):
        t += max(gauss_hypergeom(i, CntKEGG_All, SAPs_all, totalSAPs), 0.0)
    return min(max(t, 0.0), 1), min(max(s, 0.0), 1)


def fisherexct(SAPs_KEGG, NoSAPs_KEGG, SAPs_all, NoSAPs_all, CntKEGG_All,
               totalSAPs, pKEGG):
    """
    Runs Fisher's exact test
    """
    ftest = fisher(SAPs_KEGG, NoSAPs_KEGG, SAPs_all, NoSAPs_all)
    probLow, probHigh = ftest.left_tail, ftest.right_tail
    return probLow, probHigh


def rtrnKEGGcENSEMBLc(inBckgrndfile, columnENSEMBLTBckgrnd, columnKEGGBckgrnd):
    """
    """
    dKEGGTENSEMBLT = {}
    for eachl in open(inBckgrndfile, 'r'):
        if eachl.strip():
            ENSEMBLT = eachl.splitlines()[0].split('\t')[columnENSEMBLTBckgrnd]
            KEGGTs = set(
                eachl.splitlines()[0].split('\t')[columnKEGGBckgrnd].split('.'))
            KEGGTs = KEGGTs.difference(set(['', 'U', 'N']))
            for KEGGT in KEGGTs:
                try:
                    dKEGGTENSEMBLT[KEGGT].add(ENSEMBLT)
                except:
                    dKEGGTENSEMBLT[KEGGT] = set([ENSEMBLT])
    ENSEMBLTGinKEGG = set.union(*dKEGGTENSEMBLT.values())
    return dKEGGTENSEMBLT, ENSEMBLTGinKEGG


def rtrnENSEMBLcSAPs(inSAPsfile, columnENSEMBLT, ENSEMBLTGinKEGG):
    """
    returns a set of the ENSEMBLT codes present in the input list and
    in the KEGG file
    """
    sENSEMBLTSAPsinKEGG = set()
    for eachl in open(inSAPsfile, 'r'):
        ENSEMBLT = eachl.splitlines()[0].split('\t')[columnENSEMBLT]
        if ENSEMBLT in ENSEMBLTGinKEGG:
            sENSEMBLTSAPsinKEGG.add(ENSEMBLT)
    return sENSEMBLTSAPsinKEGG


def rtrnCounts(dKEGGTENSEMBLT, ENSEMBLTGinKEGG, sENSEMBLTSAPsinKEGG, statsTest):
    """
    returns a list of the ENSEMBLT codes present in the input list and
    in the KEGG file. The pathways in this list are: 'Go Term','# Genes in
    the KEGG Term','# Genes in the list and in the KEGG Term','Enrichement
    of the KEGG Term for genes in the input list','Genes in the input list
    present in the KEGG term'
    """
    totalSAPs = len(ENSEMBLTGinKEGG)
    SAPs_all = len(sENSEMBLTSAPsinKEGG)
    NoSAPs_all = totalSAPs - SAPs_all
    pKEGG = SAPs_all / float(totalSAPs)
    # ~
    lp = len(dKEGGTENSEMBLT)
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
    for echKEGGT in dKEGGTENSEMBLT:
        cnt += 1
        CntKEGG_All = len(dKEGGTENSEMBLT[echKEGGT])
        SAPs_KEGG = len(
            dKEGGTENSEMBLT[echKEGGT].intersection(sENSEMBLTSAPsinKEGG))
        NoSAPs_KEGG = CntKEGG_All - SAPs_KEGG
        probLow, probHigh = ptest(SAPs_KEGG, NoSAPs_KEGG, SAPs_all, NoSAPs_all,
                                  CntKEGG_All, totalSAPs, pKEGG)
        ltfreqs.append(
            [(SAPs_KEGG / Decimal(CntKEGG_All)), SAPs_KEGG, probLow, probHigh,
             echKEGGT])
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
        description='Returns the count of genes in KEGG categories and their statistical overrrepresentation, from a list of genes and an background file (i.e. plane text with ENSEMBLT and KEGG pathways).')
    parser.add_argument('--input', metavar='input TXT file', type=str,
                        help='the input file with the table in txt format.',
                        required=True)
    parser.add_argument('--inBckgrndfile', metavar='input TXT file', type=str,
                        help='the input file with the background table in txt format.',
                        required=True)
    parser.add_argument('--output', metavar='output TXT file', type=str,
                        help='the output file with the table in txt format.',
                        required=True)
    parser.add_argument('--columnENSEMBLT', metavar='column number', type=int,
                        help='column with the ENSEMBL transcript code in the input file.',
                        required=True)
    parser.add_argument('--columnENSEMBLTBckgrnd', metavar='column number',
                        type=int,
                        help='column with the ENSEMBL transcript code in the background file.',
                        required=True)
    parser.add_argument('--columnKEGGBckgrnd', metavar='column number',
                        type=int,
                        help='column with the KEGG pathways in the background file.',
                        required=True)
    parser.add_argument('--statsTest', metavar='input TXT file', type=str,
                        help='statistical test to compare KEGG pathways (i.e. fisher, hypergeometric, binomial).',
                        required=True)

    args = parser.parse_args()

    inSAPsfile = args.input
    inBckgrndfile = args.inBckgrndfile
    saleKEGGPCount = args.output
    columnENSEMBLT = args.columnENSEMBLT
    columnENSEMBLTBckgrnd = args.columnENSEMBLTBckgrnd
    columnKEGGBckgrnd = args.columnKEGGBckgrnd
    statsTest = args.statsTest
    columnENSEMBLT -= 1
    columnENSEMBLTBckgrnd -= 1
    columnKEGGBckgrnd = -1
    # ~
    dKEGGTENSEMBLT, ENSEMBLTGinKEGG = rtrnKEGGcENSEMBLc(inBckgrndfile,
                                                        columnENSEMBLTBckgrnd,
                                                        columnKEGGBckgrnd)
    sENSEMBLTSAPsinKEGG = rtrnENSEMBLcSAPs(inSAPsfile, columnENSEMBLT,
                                           ENSEMBLTGinKEGG)
    outl = rtrnCounts(dKEGGTENSEMBLT, ENSEMBLTGinKEGG, sENSEMBLTSAPsinKEGG,
                      statsTest)
    # ~
    saleKEGGPCount = open(saleKEGGPCount, 'w')
    saleKEGGPCount.write('\n'.join(outl))
    saleKEGGPCount.close()
    # ~
    return 0


if __name__ == '__main__':
    main()
