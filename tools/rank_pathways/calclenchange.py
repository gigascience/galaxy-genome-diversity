#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       calclenchange.py
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

import argparse, mechanize, os, sys
from decimal import Decimal, getcontext
from xml.etree.ElementTree import ElementTree, tostring
import networkx as nx
from copy import copy

# method to rank the the pthways by mut. freq.
def rankdN(ltfreqs):
    ordvals = sorted(ltfreqs)  # sort and reverse freqs.
    # ~
    outrnk = []
    tmpChng0, tmpOri, tmpMut, tmpPthw = ordvals.pop()  # the highest possible value
    if tmpOri == 'C':
        if tmpMut != 'C':
            tmpChng0 = 'C-%s' % tmpMut
        else:
            tmpChng0 = Decimal('0')
    crank = 1
    outrnk.append(
        [str(tmpChng0), str(tmpOri), str(tmpMut), str(crank), tmpPthw])
    totalnvals = len(ordvals)
    cnt = 0
    while totalnvals > cnt:
        cnt += 1
        tmpChng, tmpOri, tmpMut, tmpPthw = ordvals.pop()
        if tmpOri == 'C':
            if tmpMut != 'C':
                tmpChng = 'C-%s' % tmpMut
            else:
                tmpChng = Decimal('0')
        if tmpChng != tmpChng0:
            crank = len(outrnk) + 1
            tmpChng0 = tmpChng
        outrnk.append(
            [str(tmpChng), str(tmpOri), str(tmpMut), str(crank), tmpPthw])
    return outrnk


# method to rank the the pthways by mut. freq.
def rankdAvr(ltfreqs):
    ordvals = sorted(ltfreqs)  # sort and reverse freqs.
    # ~
    outrnk = {}
    tmpChng0, tmpOri, tmpMut, tmpPthw = ordvals.pop()  # the highest possible value
    if tmpOri == 'I':
        if tmpMut != 'I':
            tmpChng0 = 'I-%s' % tmpMut
        else:
            tmpChng0 = Decimal('0')
    crank = 1
    outrnk[tmpPthw] = '\t'.join(
        [str(tmpChng0), str(tmpOri), str(tmpMut), str(crank)])
    totalnvals = len(ordvals)
    cnt = 0
    while totalnvals > cnt:
        cnt += 1
        tmpChng, tmpOri, tmpMut, tmpPthw = ordvals.pop()
        if tmpOri == 'I':
            if tmpMut != 'I':
                tmpChng = 'I-%s' % tmpMut
            else:
                tmpChng = Decimal('0')
        if tmpChng != tmpChng0:
            crank = len(outrnk) + 1
            tmpChng0 = tmpChng
        outrnk[tmpPthw] = '\t'.join(
            [str(tmpChng), str(tmpOri), str(tmpMut), str(crank)])
    return outrnk


# this method takes as input a list of pairs of edges(beginNod,endNod) and returns a list of nodes with indegree 0 and outdegree 0
def returnstartanendnodes(edges):
    listID0st = set()  # starts
    listOD0en = set()  # end
    for beginNod, endNod in edges:  # O(n)
        listID0st.add(beginNod)
        listOD0en.add(endNod)
    startNdsID0 = listID0st.difference(listOD0en)
    endNdsOD0 = listOD0en.difference(listID0st)
    return startNdsID0, endNdsOD0


# ~ Method to return nodes and edges
def returnNodesNEdgesfKXML(fpthwKGXML):
    # ~
    tree = ElementTree()
    ptree = tree.parse(fpthwKGXML)
    # ~
    title = ptree.get('title')
    prots = ptree.findall('entry')
    reactns = ptree.findall('reaction')
    # ~
    edges, ndstmp = set(), set()
    nreactns = len(reactns)
    cr = 0  # count reacts
    while nreactns > cr:
        cr += 1
        reactn = reactns.pop()
        mainid = reactn.get('id')
        ndstmp.add(mainid)  # add node
        reacttyp = reactn.get('type')
        sbstrts = reactn.findall('substrate')
        while len(sbstrts) > 0:
            csbstrt = sbstrts.pop()
            csbtsid = csbstrt.get('id')
            ndstmp.add(csbtsid)  # add node
            if reacttyp == 'irreversible':
                edges.add((csbtsid, mainid))  # add edges
            elif reacttyp == 'reversible':
                edges.add((mainid, csbtsid))  # add edges
                edges.add((csbtsid, mainid))  # add edges
        # ~
        prdcts = reactn.findall('product')
        while len(prdcts) > 0:
            prdct = prdcts.pop()
            prodctid = prdct.get('id')
            ndstmp.add(prodctid)  # add node
            if reacttyp == 'irreversible':
                edges.add((mainid, prodctid))  # add edges
            elif reacttyp == 'reversible':
                edges.add((mainid, prodctid))  # add edges
                edges.add((prodctid, mainid))  # add edges
    # ~ Nodes
    nprots = len(prots)
    cp = 0  # count prots
    dnodes = {}
    while nprots > cp:
        cp += 1
        prot = prots.pop()
        tmpProtnm = prot.get('id')
        if tmpProtnm in ndstmp:
            dnodes[prot.get('id')] = set(
                prot.get('name').split())  # each genename for each Id
    return dnodes, edges, title


# ~ make calculation on pathways
def rtrnAvrgLen(edges, strNds, endNds):
    wG = nx.DiGraph()  # reference graph
    wG.add_edges_from(edges)
    dPairsSrcSnks = nx.all_pairs_shortest_path_length(
        wG)  # dictionary between sources and sink and length
    nstartNdsID0 = len(strNds)
    cstrtNds = 0
    nPaths = 0
    lPathLen = []
    while nstartNdsID0 > cstrtNds:
        cStartNd = strNds.pop()  # current start node
        dEndNdsLen = dPairsSrcSnks.pop(cStartNd)
        for cendNd in dEndNdsLen:
            if cendNd in endNds:
                lPathLen.append(dEndNdsLen[cendNd])
                nPaths += 1
        cstrtNds += 1
    AvrgPthLen = 0
    if nPaths != 0:
        AvrgPthLen = Decimal(sum(lPathLen)) / Decimal(str(nPaths))
    return nPaths, AvrgPthLen


def main():
    parser = argparse.ArgumentParser(
        description='Rank pathways based on the change in length and number of paths connecting sources and sinks.')
    parser.add_argument('--loc_file', metavar='correlational database',
                        type=str, help='correlational database')
    parser.add_argument('--species', metavar='species name', type=str,
                        help='the species of interest in loc_file')
    parser.add_argument('--output', metavar='output TXT file', type=str,
                        help='the output file with the table in txt format. Column 1 is the diference between column 2 and column 3, Column 2 is the pathway average length (between sources and sinks) including the genes in the input list, Column 3 is the pathway average length EXCLUDING the genes in the input list, Column 4 is the rank based on column 1. Column 5 is the diference between column 6 and column 7, Column 6 is the number of paths between sources and sinks, including the genes in the input list, Column 7 is the number of paths between sources and sinks EXCLUDING the genes in the input list, Column 8 is the rank based on column 5. Column 9 I the pathway name')
    parser.add_argument('--posKEGGclmn', metavar='column number', type=int,
                        help='the column with the KEGG pathway code/name')
    parser.add_argument('--KEGGgeneposcolmn', metavar='column number', type=int,
                        help='column with the KEGG gene code')
    parser.add_argument('--input', metavar='input TXT file', type=str,
                        help='the input file with the table in txt format')
    # ~
    # ~Open arguments
    class C(object):
        pass

    fulargs = C()
    parser.parse_args(sys.argv[1:], namespace=fulargs)
    # test input vars
    inputf, loc_file, species, output, posKEGGclmn, Kgeneposcolmn = fulargs.input, fulargs.loc_file, fulargs.species, fulargs.output, fulargs.posKEGGclmn, fulargs.KEGGgeneposcolmn
    posKEGGclmn -= 1  # correct pos
    Kgeneposcolmn -= 1
    # ~ Get the extra variables
    crDB = [x.split() for x in open(loc_file).read().splitlines() if
            x.split()[0] == species][0]
    sppPrefx, dinput = crDB[1], crDB[2]
    # ~ set decimal positions
    getcontext().prec = 3
    # make a dictionary of valid genes
    dKEGGcPthws = dict([(x.split('\t')[Kgeneposcolmn], set(
        [y.split('=')[0] for y in x.split('\t')[posKEGGclmn].split('.')])) for x
                        in open(inputf).read().splitlines()[1:] if x.strip()])
    sdGenes = set([x for x in dKEGGcPthws.keys() if x.find('.') > -1])
    while True:  # to crrect names with more than one gene
        try:
            mgenes = sdGenes.pop()
            pthwsAssotd = dKEGGcPthws.pop(mgenes)
            mgenes = mgenes.split('.')
            for eachg in mgenes:
                dKEGGcPthws[eachg] = pthwsAssotd
        except:
            break
    # ~
    lPthwsF = [x for x in os.listdir(dinput) if x.find('.xml') > -1 if
               x not in ['cfa04070.xml']]
    nPthws = len(lPthwsF)
    cPthw = 0
    lPthwPthN = []  # the output list for number of paths
    lPthwPthAvr = []  # the output list for the length of paths
    # ~
    while cPthw < nPthws:
        cPthw += 1
        KEGGpathw = lPthwsF.pop()
        comdKEGGpathw = KEGGpathw.split('.')[0]
        tmpddGenrcgenPresent = set()
        sKEGGc = dKEGGcPthws.keys()
        lsKEGGc = len(sKEGGc)
        ctPthw = 0
        while ctPthw < lsKEGGc:  # to save memory
            eachK = sKEGGc.pop()
            alPthws = dKEGGcPthws[eachK]
            if comdKEGGpathw in alPthws:
                tmpddGenrcgenPresent.add(':'.join([sppPrefx, eachK]))
            ctPthw += 1
        # ~ Make graph calculations
        dnodes, edges, title = returnNodesNEdgesfKXML(
            open(os.path.join(dinput, KEGGpathw)))
        startNdsID0, endNdsOD0 = returnstartanendnodes(edges)
        startNdsOri = copy(startNdsID0)
        # ~
        nPaths = 'C'  # stands for circuit
        AvrgPthLen = 'I'  # stand for infinite
        if len(startNdsID0) > 0 and len(endNdsOD0) > 0:
            nPaths, AvrgPthLen = rtrnAvrgLen(edges, startNdsID0, endNdsOD0)
        # ~ work with the genes in the list
        genestodel = set()
        lnodes = len(dnodes)
        sNds = set(dnodes)
        ctPthw = 0
        while ctPthw < lnodes:
            ctPthw += 1
            cNod = sNds.pop()
            sgenes = dnodes.pop(cNod)
            if len(sgenes.intersection(tmpddGenrcgenPresent)) == len(sgenes):
                genestodel.add(cNod)
        # ~ del nodes from graph edges
        wnPaths, wAvrgPthLen = copy(nPaths), copy(AvrgPthLen)
        if len(genestodel) > 0:
            wedges = set(
                [x for x in edges if len(set(x).intersection(genestodel)) == 0])
            wstartNds, wendNds = returnstartanendnodes(wedges)
            if nPaths != 'C':
                wstartNds = [x for x in wstartNds if x in startNdsOri]
                wendNds = [x for x in wendNds if x in endNdsOD0]
            if len(wstartNds) > 0 and len(wendNds) > 0:
                wnPaths, wAvrgPthLen = rtrnAvrgLen(wedges, wstartNds, wendNds)
        # ~ Calculate the differences
        orNP, mutNP, oriLen, mutLen = nPaths, wnPaths, AvrgPthLen, wAvrgPthLen
        if nPaths == 'C':
            orNP = Decimal('1000')
            oriLen = Decimal('1000')
        if wnPaths == 'C':
            mutNP = Decimal('1000')
            mutLen = Decimal('1000')
        lPthwPthN.append([orNP - mutNP, nPaths, wnPaths, '='.join(
            [comdKEGGpathw, title])])  # print nPaths,AvrgPthLen
        lPthwPthAvr.append([oriLen - mutLen, AvrgPthLen, wAvrgPthLen, '='.join(
            [comdKEGGpathw, title])])  # print nPaths,AvrgPthLen
    doutrnkPthN = rankdN(lPthwPthN)
    doutrnkPthAvr = rankdAvr(lPthwPthAvr)
    # ~
    sall = ['\t'.join([doutrnkPthAvr[x[4]], '\t'.join(x)]) for x in doutrnkPthN]
    salef = open(output, 'w')
    salef.write('\n'.join(sall))
    salef.close()
    return 0


if __name__ == '__main__':
    main()
