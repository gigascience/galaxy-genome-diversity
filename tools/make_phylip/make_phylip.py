#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       mkFastas.py
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
import glob
import os
import shutil
from Population import Population

def revseq(seq):
    seq=list(seq)
    seq.reverse()
    seq=''.join(seq)
    return seq
    
def revComp(allPop):
    dAllCompAll={'A':'T','T':'A','C':'G','G':'C','N':'N','M':'K','K':'M','R':'Y','Y':'R','W':'W','S':'S'}
    allPopsComp=dAllCompAll[allPop]
    return allPopsComp

def rtrnCons(ntA,ntB):
    srtdPairs=''.join(sorted([ntA,ntB]))
    dpairsCons={'AC':'M', 'AG':'R', 'AT':'W', 'CG':'S', 'CT':'Y', 'GT':'K', 'AN':'A', 'CN':'C', 'GN':'G', 'NT':'T'}
    cons=dpairsCons[srtdPairs]
    return cons

def rtrnFxdChrPos(inSNPf,dPopsinSNPfPos,pxchrx,pxpos,pxntA,pxntB,fulldChrdPosdPopsAlllsInit=False,cvrgTreshold=False,indvlsPrctTrshld=False):
    """
    """
    dChrdPosdPopsAlllsInit={}
    seqref=[]
    for eachl in open(inSNPf,'r'):
        if eachl.strip() and eachl[0]!='#':
            fllInfoSplt=eachl.splitlines()[0].split('\t')
            chrx=fllInfoSplt[pxchrx]
            pos=int(fllInfoSplt[pxpos])
            ntA=fllInfoSplt[pxntA]
            ntB=fllInfoSplt[pxntB]
            seqref.append([pos,ntA])
            dPopsAllls={}
            if fulldChrdPosdPopsAlllsInit:
                #~ 
                cntIndv=0
                #
                try:
                    fulldPopsAllls=fulldChrdPosdPopsAlllsInit[chrx][pos]
                except:
                    fulldPopsAllls=dict([(echPop,ntA) for echPop in dPopsinSNPfPos])            
                #
                for eachPop in dPopsinSNPfPos:
                    clmnCvrg=dPopsinSNPfPos[eachPop]
                    if clmnCvrg:
                        eachPopCvrg=int(fllInfoSplt[clmnCvrg])
                    else:
                        #~ eachPopCvrg=0
                        eachPopCvrg=cvrgTreshold
                    if eachPopCvrg>=cvrgTreshold:
                        dPopsAllls[eachPop]=fulldPopsAllls[eachPop]
                        cntIndv+=1
                    else:
                        dPopsAllls[eachPop]='N'        
                #~ 
                if indvlsPrctTrshld>(cntIndv/float(len(dPopsinSNPfPos))):
                    dPopsAllls=dict([(echPop,'N') for echPop in dPopsinSNPfPos])
            #~ 
            else:
                for eachPop in dPopsinSNPfPos:
                    if dPopsinSNPfPos[eachPop]:
                        eachPopAll=int(fllInfoSplt[dPopsinSNPfPos[eachPop]])
                        if eachPopAll==0:
                            dPopsAllls[eachPop]=ntB
                        elif eachPopAll==2:
                            dPopsAllls[eachPop]=ntA
                        elif eachPopAll==1:
                            dPopsAllls[eachPop]=rtrnCons(ntA,ntB)
                        else:
                            dPopsAllls[eachPop]='N'
                    else:
                        dPopsAllls[eachPop]=ntA
            try:
                dChrdPosdPopsAlllsInit[chrx][pos]=dPopsAllls
            except:
                dChrdPosdPopsAlllsInit[chrx]={pos:dPopsAllls}
    #~ 
    seqref.sort()
    startExs=[seqref[0][0]]
    endExs=[seqref[-1][0]+1]
    seqref=''.join(x[1] for x in seqref)
    #~ 
    return dChrdPosdPopsAlllsInit,seqref,chrx,startExs,endExs


def rtrndENSEMBLTseq(inCDSfile,inUCSCfile,fchrClmn,txStartClmn,txEndClmn,strandClmn,geneNameClmn,startExsClmn,endExsClmn,cdsStartClmn,cdsEndClmn):
    """
    """
    dENSEMBLTchrxStEndEx={}
    dChrdStrtEndENSEMBLT={}
    for eachl in open(inUCSCfile,'r'):
        if eachl.strip():
            rvrse=False
            allVls=eachl.split('\t')
            txStart=allVls[txStartClmn]
            txEnd=allVls[txEndClmn]
            ENSEMBLT=allVls[geneNameClmn]
            strand=allVls[strandClmn]
            chrx=allVls[fchrClmn]
            if cdsStartClmn and cdsEndClmn:
                cdsStart=allVls[cdsStartClmn]
                cdsEnd=allVls[cdsEndClmn]
            if startExsClmn and endExsClmn:
                startExs=allVls[startExsClmn]
                endExs=allVls[endExsClmn]
            if strand=='-':
                rvrse=True
            try:
                dChrdStrtEndENSEMBLT[chrx][int(txStart),int(txEnd)]=ENSEMBLT
            except:
                try:
                    dChrdStrtEndENSEMBLT[chrx]={(int(txStart),int(txEnd)):ENSEMBLT}
                except:
                    dChrdStrtEndENSEMBLT={chrx:{(int(txStart),int(txEnd)):ENSEMBLT}}
            #~ 
            if cdsStartClmn and cdsEndClmn and startExsClmn and endExsClmn:        
                startExs,endExs=rtrnExnStarEndCorrc(startExs,endExs,cdsStart,cdsEnd)
            else:
                startExs,endExs=[int(txStart)],[int(txEnd)]
            dENSEMBLTchrxStEndEx[ENSEMBLT]=(chrx,startExs,endExs,rvrse)
    #~ 
    dENSEMBLTseq={}
    ENSEMBLTseqs=[(x.splitlines()[0],''.join(x.splitlines()[1:])) for x in open(inCDSfile).read().split('>') if x.strip()]
    for ENSEMBLT,seq in ENSEMBLTseqs:
        dENSEMBLTseq[ENSEMBLT]=seq
    #~ 
    dENSEMBLTseqChrStEnEx={}
    for ENSEMBLT in dENSEMBLTchrxStEndEx:
        chrx,startExs,endExs,rvrse=dENSEMBLTchrxStEndEx[ENSEMBLT]
        addEseqChrStEnEx=True
        try:
            seq=dENSEMBLTseq[ENSEMBLT]
            if rvrse:
                seq=revseq(seq)
        except:
            addEseqChrStEnEx=False
        if addEseqChrStEnEx:
            dENSEMBLTseqChrStEnEx[ENSEMBLT]=(seq,chrx,startExs,endExs,rvrse)
    return dENSEMBLTseqChrStEnEx,dChrdStrtEndENSEMBLT
    

def rtrnFxdChrPosinCodReg(dChrdStrtEndENSEMBLT,dChrdPosdPopsAlllsInit):
    """
    """
    dENSEMBLTChrPosdAlls={}
    dChrPosdPopsAllls={}
    todel=set(dChrdPosdPopsAlllsInit.keys()).difference(set(dChrdStrtEndENSEMBLT.keys()))
    for x in todel:
        x=dChrdPosdPopsAlllsInit.pop(x)
    #---
    while len(dChrdPosdPopsAlllsInit)>0:
        chrx=dChrdPosdPopsAlllsInit.keys()[0]
        dStrtEndENSEMBLT=dChrdStrtEndENSEMBLT.pop(chrx)
        dPosdPopsAllls=dChrdPosdPopsAlllsInit.pop(chrx)
        #~ 
        srtdStrtEndENSEMBLT=sorted(dStrtEndENSEMBLT.keys())
        srtdPosdPopsAllls=sorted(dPosdPopsAllls.keys())
        #~
        pos=srtdPosdPopsAllls.pop(0)
        strt,end=srtdStrtEndENSEMBLT.pop(0)
        ENSEMBLT=dStrtEndENSEMBLT[strt,end]
        dPopsAllls=dPosdPopsAllls[pos]
        keePloop=True
        #~ 
        while keePloop:
            if strt<=pos<=end:
                for tmpstrt,tmpend in [(strt,end)]+srtdStrtEndENSEMBLT:
                    if tmpstrt<=pos<=tmpend:
                        dPopsAllls=dPosdPopsAllls[pos]
                        dChrPosdPopsAllls[chrx,pos]=dPopsAllls
                        try:
                            dENSEMBLTChrPosdAlls[ENSEMBLT][chrx,pos]=dPopsAllls
                        except:
                            dENSEMBLTChrPosdAlls[ENSEMBLT]={(chrx,pos):dPopsAllls}
                    else:
                        continue
                if len(srtdPosdPopsAllls)>0:
                    pos=srtdPosdPopsAllls.pop(0)
                    dPopsAllls=dPosdPopsAllls[pos]
                else:
                    keePloop=False
            #~ 
            elif pos<=strt:
                if len(srtdPosdPopsAllls)>0:
                    pos=srtdPosdPopsAllls.pop(0)
                    dPopsAllls=dPosdPopsAllls[pos]
                else:
                    keePloop=False
            else:
                if len(srtdStrtEndENSEMBLT)>0:
                    strt,end=srtdStrtEndENSEMBLT.pop(0)
                    ENSEMBLT=dStrtEndENSEMBLT[strt,end]
                else:
                    keePloop=False
    return dENSEMBLTChrPosdAlls,dChrPosdPopsAllls

def rtrnExnStarEndCorrc(startExs,endExs,cdsStart,cdsEnd):
    """
    """
    cdsStart,cdsEnd=int(cdsStart),int(cdsEnd)
    crrctdstartExs=set([int(x) for x in startExs.split(',') if x.strip()])
    crrctdendExs=set([int(x) for x in endExs.split(',') if x.strip()])
    crrctdstartExs.add(cdsStart)
    crrctdendExs.add(cdsEnd)
    sStartDel=set()
    sEndDel=set()
    #~ 
    for echvl in crrctdstartExs:
        if echvl<cdsStart or echvl>cdsEnd:
            sStartDel.add(echvl)
    #~ 
    for echvl in crrctdendExs:
        if echvl<cdsStart or echvl>cdsEnd:
            sEndDel.add(echvl)
    #~ 
    return sorted(crrctdstartExs.difference(sStartDel)),sorted(crrctdendExs.difference(sEndDel))

def rtrndPopsFasta(seq,chrx,startExs,endExs,rvrse,dChrPosdPopsAllls,ENSEMBLT):
    """
    """
    exnIntrvl=zip(startExs,endExs)
    CDSinitPos=exnIntrvl[0][0]
    dexnIntrvlSeq={}    
    for exStart,exEnd in exnIntrvl:
        lenEx=exEnd-exStart
        dexnIntrvlSeq[exStart,exEnd]=seq[:lenEx]
        seq=seq[lenEx:]
        
    ldexnIntrvlSeq=len(dexnIntrvlSeq)
    #~ 
    dPopsFasta={}
    #~ 
    strePos=set()
    dStrePosAbsPos={}
    tmpAcmltdPos=0
    #~ 
    exStart,exEnd=sorted(dexnIntrvlSeq.keys())[0]
    seq=dexnIntrvlSeq.pop((exStart,exEnd))
    chrx,pos=sorted(dChrPosdPopsAllls.keys())[0]
    dPopsAllls=dChrPosdPopsAllls.pop((chrx,pos))
    tmpdPopsFasta=dict([(x,list(seq)) for x in dPopsAllls])
    cntExns=0
    while True:
        if  exStart<=pos<=exEnd-1:
            relPos=tmpAcmltdPos+pos-exStart
            strePos.add(relPos)
            dStrePosAbsPos[relPos]=pos
            for echPop in tmpdPopsFasta:
                allPop=dPopsAllls[echPop]
                if rvrse:
                    allPop=revComp(allPop)
                tmpdPopsFasta[echPop][pos-exStart]=allPop
            if len(dChrPosdPopsAllls)>0:
                chrx,pos=sorted(dChrPosdPopsAllls.keys())[0]
                dPopsAllls=dChrPosdPopsAllls.pop((chrx,pos))
            else:
                pos=endExs[-1]+100#max pos of exns
        elif pos<exStart:
            if len(dChrPosdPopsAllls)>0:
                chrx,pos=sorted(dChrPosdPopsAllls.keys())[0]
                dPopsAllls=dChrPosdPopsAllls.pop((chrx,pos))
            else:
                pos=endExs[-1]+100#max pos of exns
        elif pos>exEnd-1:# or len(dChrPosdPopsAllls)==0:
            for echPop in tmpdPopsFasta:
                try:
                    dPopsFasta[echPop]+=''.join(tmpdPopsFasta[echPop])
                except:
                    dPopsFasta[echPop]=''.join(tmpdPopsFasta[echPop])
            cntExns+=1
            tmpAcmltdPos+=len(seq)
            if len(dexnIntrvlSeq)>0:
                exStart,exEnd=sorted(dexnIntrvlSeq.keys())[0]
                seq=dexnIntrvlSeq.pop((exStart,exEnd))
                tmpdPopsFasta=dict([(x,list(seq)) for x in dPopsAllls])
            else:
                break
    if ldexnIntrvlSeq!=cntExns:
        for echPop in tmpdPopsFasta:
            dPopsFasta[echPop]+=''.join(tmpdPopsFasta[echPop])
    #~ 
    lchrStartexEndpos=[]
    if rvrse:
        dPopsFasta=dict([(echPop,revseq(dPopsFasta[echPop])) for echPop in dPopsFasta])#[echPop]+=''.join(tmpdPopsFasta[echPop])
        for ePos in strePos:
            lchrStartexEndpos.append('\t'.join([ENSEMBLT,chrx,str(tmpAcmltdPos-ePos-1),str(dStrePosAbsPos[ePos])]))
    else:
        for ePos in strePos:
            lchrStartexEndpos.append('\t'.join([ENSEMBLT,chrx,str(ePos),str(dStrePosAbsPos[ePos])]))
    #~ 
    return dPopsFasta,lchrStartexEndpos

def rtrnSeqVars(dENSEMBLTseqChrStEnEx,dENSEMBLTChrPosdAlls):
    """
    """
    dENSEMBLTPopsFasta={}
    lchrStartexEndposAll=[]
    #~ 
    sENSEMBLTcmmn=set(dENSEMBLTChrPosdAlls.keys()).intersection(set(dENSEMBLTseqChrStEnEx.keys()))#sENSEMBLTcmmn between UCSC and ENSEMBLE
    #~ 
    for ENSEMBLT in sENSEMBLTcmmn:
        seq,chrx,startExs,endExs,rvrse=dENSEMBLTseqChrStEnEx[ENSEMBLT]
        dChrPosdPopsAllls=dENSEMBLTChrPosdAlls[ENSEMBLT]
        if len(startExs)>0 and len(endExs)>0:
            dPopsFasta,lchrStartexEndpos=rtrndPopsFasta(seq,chrx,startExs,endExs,rvrse,dChrPosdPopsAllls,ENSEMBLT)
            lchrStartexEndposAll.extend(lchrStartexEndpos)
            if dPopsFasta:#to correct a bug of the input table, in cases in which endExons<startExn (!). See ENSCAFT00000000145 (MC4R) in canFam2 for example.
                dENSEMBLTPopsFasta[ENSEMBLT]=dPopsFasta
    return dENSEMBLTPopsFasta,lchrStartexEndposAll



def rtrnPhy(dPopsFasta,ENSEMBLT):
    """
    """
    dPopsFormPhy={}
    for eachPop in dPopsFasta:
        hader='%s'%eachPop
        #~ hader='>%s'%eachPop
        seq=dPopsFasta[eachPop]
        formtd='\t'.join([hader,seq])
        #~ formtd='\n'.join([hader,seq])
        dPopsFormPhy[eachPop]=formtd
    #~ 
    return dPopsFormPhy,len(seq)

def wrapSeqsFasta(dENSEMBLTPopsFasta,sPopsIntrst):
    """
    """
    ENSEMBLTKaKs=[]
    nonHeader=True
    cnt=0
    lENSEMBLT=len(dENSEMBLTPopsFasta)
    #~     
    for ENSEMBLT in sorted(dENSEMBLTPopsFasta.keys()):
        cnt+=1
        dPopsFasta=dENSEMBLTPopsFasta[ENSEMBLT]
        dPopsFormPhy,lenseq=rtrnPhy(dPopsFasta,ENSEMBLT)
        #~ 
        seqPMLformat=['%s %s'%(len(dPopsFormPhy),lenseq)]#generate new PHYML sequence
        #~ seqPMLformat=[]#generate new PHYML sequence
        for namex in sorted(sPopsIntrst):
            seqPMLformat.append(dPopsFormPhy[namex])
        #~ 
        outFastaf=open('%s.phy'%ENSEMBLT,'w')
        outFastaf.write('\n'.join(seqPMLformat))
        outFastaf.close()
        #~ 
    return 0

def pos_dict(gd_indivs_file, input_type):
    rv = {}

    p = Population()
    p.from_population_file(gd_indivs_file)

    for tag in p.tag_list():
        column, name = tag.split(':')
        column = int(column) - 1

        if input_type == 'gd_genotype':
            column -= 2

        rv[name] = column

    return rv

def main():
    #~ 
    #~bpython mkPhyl.py --input=colugo_mt_Galaxy_genotypes.txt --chrClmn=0 --posClmn=1 --refClmn=2 --altrClmn=3 --output=out.d --gd_indivs=genotypes.gd_indivs --inputCover=colugo_mt_Galaxy_coverage.txt --gd_indivs_cover=coverage.gd_indivs --cvrgTreshold=0 --chrClmnCvrg=0 --posClmnCvrg=1 --refClmnCvrg=2 --altrClmnCvrg=3 --indvlsPrctTrshld=0
    parser = argparse.ArgumentParser(description='Returns the count of genes in KEGG categories and their statistical overrrepresentation, from a list of genes and an background file (i.e. plane text with ENSEMBLT and KEGG pathways).')
    parser.add_argument('--input',metavar='input gd_snp file',type=str,help='the input file with the table in gd_snp/gd_genotype format.',required=True)    
    parser.add_argument('--input_type',metavar='input type',type=str,help='the input file type (gd_snp or gd_genotype)',required=True)    
    parser.add_argument('--chrClmn',metavar='int',type=int,help='the column with the chromosome.',required=True)
    parser.add_argument('--posClmn',metavar='int',type=int,help='the column with the SNPs position.',required=True)
    parser.add_argument('--refClmn',metavar='int',type=int,help='the column with the reference nucleotide.',required=True)
    parser.add_argument('--altrClmn',metavar='int',type=int,help='the column with the derived nucleotide.',required=True)
    parser.add_argument('--output',metavar='output',type=str,help='the output',required=True)
    parser.add_argument('--output_id',metavar='int',type=int,help='the output id',required=True)
    parser.add_argument('--gd_indivs',metavar='input gd_indivs file',type=str,help='the input reference species columns in the input file.',required=True)
    #~ 
    parser.add_argument('--inputCover',metavar='input gd_snp cover file',type=str,help='the input file with the table in gd_snp/gd_genotype cover format.',required=False,default=False)
    parser.add_argument('--inputCover_type',metavar='input cover type',type=str,help='the cover input file type (gd_snp or gd_genotype)',required=False,default=False)
    parser.add_argument('--gd_indivs_cover',metavar='input gd_indivs file',type=str,help='the input reference species columns in the input cover file.',required=False,default=False)
    parser.add_argument('--cvrgTreshold',metavar='input coverage threshold',type=int,help='the coverage threshold above which nucleotides are included, else "N".',required=False,default=False)
    parser.add_argument('--chrClmnCvrg',metavar='int',type=int,help='the column with the chromosome in the input coverage file.',required=False,default=False)
    parser.add_argument('--posClmnCvrg',metavar='int',type=int,help='the column with the SNPs position in the input coverage file.',required=False,default=False)
    parser.add_argument('--refClmnCvrg',metavar='int',type=int,help='the column with the reference nucleotide in the input coverage file.',required=False,default=False)
    parser.add_argument('--altrClmnCvrg',metavar='int',type=int,help='the column with the derived nucleotide in the input coverage file.',required=False,default=False)
    parser.add_argument('--indvlsPrctTrshld',metavar='int',type=float,help='the percentage of individual above which nucleotides are included, else "N".',required=False,default=False)
    #~ 
    parser.add_argument('--sequence',metavar='input fasta file',type=str,help='the input file with the sequence whose SNPs are in the input.',required=False,default=False)
    parser.add_argument('--gene_info',metavar='input interval file',type=str,help='the input interval file with the the information on the genes.',required=False,default=False)
    parser.add_argument('--fchrClmn',metavar='int',type=int,help='the column with the chromosome in the gene_info file.',required=False,default=False)
    parser.add_argument('--txStartClmn',metavar='int',type=int,help='the column with the transcript start column in the gene_info file.',required=False,default=False)
    parser.add_argument('--txEndClmn',metavar='int',type=int,help='the column with the transcript end column in the gene_info file.',required=False,default=False)
    parser.add_argument('--strandClmn',metavar='int',type=int,help='the column with the strand column in the gene_info file.',required=False,default=False)
    parser.add_argument('--geneNameClmn',metavar='int',type=int,help='the column with the gene name column in the gene_info file.',required=False,default=False)
    parser.add_argument('--cdsStartClmn',metavar='int',type=int,help='the column with the coding start column in the gene_info file.',required=False,default=False)
    parser.add_argument('--cdsEndClmn',metavar='int',type=int,help='the column with the coding end column in the gene_info file.',required=False,default=False)
    parser.add_argument('--startExsClmn',metavar='int',type=int,help='the column with the exon start positions column in the gene_info file.',required=False,default=False)
    parser.add_argument('--endExsClmn',metavar='int',type=int,help='the column with the exon end positions column in the gene_info file.',required=False,default=False)
 
    args = parser.parse_args()

    inSNPf = args.input
    inSNPf_type = args.input_type
    outfile = args.output
    outfile_id = args.output_id
    gd_indivs = args.gd_indivs    
    pxchrx = args.chrClmn
    pxpos = args.posClmn
    pxntA = args.refClmn
    pxntB = args.altrClmn
    
    
    inCDSfile = args.sequence
    inUCSCfile = args.gene_info
    fchrClmn = args.fchrClmn#chromosome column
    txStartClmn = args.txStartClmn#transcript start column
    txEndClmn = args.txEndClmn#transcript end column
    strandClmn = args.strandClmn#strand column
    geneNameClmn = args.geneNameClmn#gene name column
    cdsStartClmn = args.cdsStartClmn#coding sequence start column
    cdsEndClmn = args.cdsEndClmn#coding sequence end column
    startExsClmn = args.startExsClmn#exons start column
    endExsClmn = args.endExsClmn#exons end column
    
    inputCover = args.inputCover
    inputCover_type = args.inputCover_type
    gd_indivs_cover = args.gd_indivs_cover
    cvrgTreshold = args.cvrgTreshold
    pxchrxCov = args.chrClmnCvrg
    pxposCov = args.posClmnCvrg
    pxntACov = args.refClmnCvrg
    pxntBCov = args.altrClmnCvrg
    indvlsPrctTrshld = args.indvlsPrctTrshld


    #print inputCover, gd_indivs_cover, cvrgTreshold
    
    assert ((inputCover and gd_indivs_cover and cvrgTreshold>=0 and indvlsPrctTrshld>=0) or (inCDSfile and inUCSCfile))

    dPopsinSNPfPos = pos_dict(gd_indivs, inSNPf_type)
    #~ dPopsinSNPfPos.update({'ref':False})
    #~ 
    sPopsIntrst=set(dPopsinSNPfPos.keys())
    dChrdPosdPopsAlllsInit,seqref,chrx,startExs,endExs=rtrnFxdChrPos(inSNPf,dPopsinSNPfPos,pxchrx,pxpos,pxntA,pxntB)#~ print '1. Getting fixed alleles information...'
    #~ dENSEMBLTseqChrStEnEx,dChrdStrtEndENSEMBLT=rtrndENSEMBLTseq(inCDSfile,inUCSCfile)
    #~
    if  inputCover and gd_indivs_cover and cvrgTreshold>=0:
        dPopsinSNPfPos_cover=dict([(eachPop,False) for eachPop in dPopsinSNPfPos.keys()])
        dPopsinSNPfPos_cover.update(pos_dict(gd_indivs_cover, inputCover_type))
        dChrdPosdPopsAlllsInit,seqref,chrx,startExs,endExs=rtrnFxdChrPos(inputCover,dPopsinSNPfPos_cover,pxchrxCov,pxposCov,pxntACov,pxntBCov,dChrdPosdPopsAlllsInit,cvrgTreshold,indvlsPrctTrshld)
        rvrse=False
        dENSEMBLTseqChrStEnEx={'tmp':(seqref,chrx,startExs,endExs,rvrse)}
        dChrdStrtEndENSEMBLT={chrx:{(startExs[0],endExs[0]):'tmp'}}
    #~ 
    elif inCDSfile and inUCSCfile:
        dENSEMBLTseqChrStEnEx,dChrdStrtEndENSEMBLT=rtrndENSEMBLTseq(inCDSfile,inUCSCfile,fchrClmn,txStartClmn,txEndClmn,strandClmn,geneNameClmn,startExsClmn,endExsClmn,cdsStartClmn,cdsEndClmn)#~ print '2. Getting transcripts and exons information...'        
    #~ 
    dENSEMBLTChrPosdAlls,dChrPosdPopsAllls=rtrnFxdChrPosinCodReg(dChrdStrtEndENSEMBLT,dChrdPosdPopsAlllsInit)#~ print '3. Getting fixed alleles in exons...'
    #~ 
    dENSEMBLTPopsFasta,lchrStartexEndposAll=rtrnSeqVars(dENSEMBLTseqChrStEnEx,dENSEMBLTChrPosdAlls)#~ print '4. Getting fasta sequences of populations...'
    #~ 
    wrapSeqsFasta(dENSEMBLTPopsFasta,sPopsIntrst)
    #~ 

    ## get a list of output files
    files = glob.glob('*.phy')

    if len(files) == 0:
        with open(outfile, 'w') as ofh:
            print >> ofh, 'No output.'
    else:
        ## the first file becomes the output
        file = files.pop(0)
        shutil.move(file, outfile)

        ## rename the rest of the files
        for file in files:
            name = file[:-4]
            name = name.replace('_', '-')
            new_filename = 'primary_{0}_{1}_visible_txt_?'.format(outfile_id, name)
            os.rename(file, new_filename)

    return 0

if __name__ == '__main__':
    main()
