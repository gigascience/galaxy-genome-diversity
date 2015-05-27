#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       multiple_to_gd_genotype.py
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
import base64
import json
import os
import sys

def fold_line(line, maxlen=200, prefix="#"):
	"""
	format hader to a 200 char max.
	"""
	line_len = len(line)
	prefix_len = len(prefix)

	if line_len + prefix_len <= maxlen:
		return '%s%s' % (prefix, line)

	lines = []
	start_idx = 0
	offset = 0

	while start_idx < line_len - 1:
		last_idx = start_idx
		idx = start_idx
		start = start_idx

		while idx != -1 and idx < maxlen + prefix_len + offset - 1:
			last_idx = idx
			idx = line.find(',', start)
			start = idx+1

		if idx == -1:
			lines.append('%s%s' % (prefix, line[start_idx:]))
			break

		lines.append('%s%s' % (prefix, line[start_idx:last_idx+1]))
		start_idx = last_idx + 1
		offset = last_idx + 1

	return '\n'.join(lines)


def formthdr(lPops,dbkey,species):
	"""
	returns a formated metadata for a gd_genotype file from a paramSet
	dictionary and a list (lPops) of individuals
	"""
	clmnsVals=', '.join(['"%sG"'%(x+1) for x in range(len(lPops))])#"'1G', '2G', ..."
	indvdls='], ['.join(['"%s", %s'%(lPops[x],x+5) for x in range(len(lPops))])#['DU23M01 Duroc domestic breed Europe', 5], ['DU23M02 Duroc domestic breed Europe', 6], ...
	obj='{"rPos": 2, "column_names": ["chr", "pos", "A", "B", %s], "scaffold": 1, "pos": 2, "dbkey": "%s", "individuals": [[%s]], "ref": 1, "species": "%s"}'%(clmnsVals,dbkey,indvdls,species)
	json_value = json.loads(obj)
	hdr = fold_line(json.dumps(json_value, separators=(',',':'), sort_keys=True))
	#~ 
	return hdr


def formthdr_gdsnp(lPops,dbkey,species):
	"""
	returns a formated metadata for a gd_genotype file from a paramSet
	dictionary and a list (lPops) of individuals
	"""
	clmnsVals=', '.join(['"%sA","%sB","%sG","%sQ"'%((x+1),(x+1),(x+1),(x+1)) for x in range(len(lPops))])#"'1G', '2G', ..."
	indvdls='], ['.join(['"%s", %s'%(lPops[x],(((x+1)*4)+2)) for x in range(len(lPops))])#['DU23M01 Duroc domestic breed Europe', 5], ['DU23M02 Duroc domestic breed Europe', 9], ...
	obj='{"rPos": 2, "column_names": ["chr", "pos", "A", "B", "Q", %s], "scaffold": 1, "pos": 2, "dbkey": "%s", "individuals": [[%s]], "ref": 1, "species": "%s"}'%(clmnsVals,dbkey,indvdls,species)
	json_value = json.loads(obj)
	hdr = fold_line(json.dumps(json_value, separators=(',',':'), sort_keys=True))
	#~ 
	return hdr


def selAnc(SNPs):
	"""	
	returns the ancestral and derived snps, and an gd_genotype-encoded
	list of SNPs/nucleotides
	"""
	dIUPAC={'AC':'M','AG':'R','AT':'W','CG':'S','CT':'Y','GT':'K'}#'N':'N','A':'A','T':'T','C':'C','G':'G',
	dNtsCnts={}
	for eSNP in SNPs:
		if eSNP!='N':
			try:
				dNtsCnts[eSNP]+=1
			except:
				dNtsCnts[eSNP]=1
	#~ 
	nAlleles=len(dNtsCnts.keys())
	assert nAlleles<=3
	if nAlleles==3:
		nonCons=[x for x in dNtsCnts.keys() if x in set(['A','T','C','G'])]
		cons=[x for x in dNtsCnts.keys() if x not in nonCons]
		assert len(nonCons)==2 and len(cons)==1 and dIUPAC[''.join(sorted(nonCons))]==''.join(cons)
		ancd=nonCons[0]
		dervd=nonCons[1]
		if dNtsCnts[dervd]>dNtsCnts[ancd]:
			ancd,dervd=dervd,ancd
	elif nAlleles==2:
		cons=set(dNtsCnts.keys()).intersection(set(dIUPAC.values()))
		assert len(cons)<=1
		if len(cons)==0:
			ancd,dervd=dNtsCnts.keys()
			if dNtsCnts[dervd]>dNtsCnts[ancd]:
				ancd,dervd=dervd,ancd
		else:
			dervd=''.join(cons)
			ancd=''.join([x for x in dNtsCnts.keys() if x!=dervd])
	else:#<=1
		ancd=''.join(dNtsCnts.keys())
		dervd='N'
	#~ 
	lSNPsEncoded=[]
	for eSNP in SNPs:
		if eSNP=='N':
			lSNPsEncoded.append('-1')
		elif eSNP==ancd:
			lSNPsEncoded.append('2')
		elif eSNP==dervd:
			lSNPsEncoded.append('0')
		else:
			lSNPsEncoded.append('1')			
	#~ 
	return ancd,dervd,lSNPsEncoded



def from_csv(in_csv,outgd_gentp,dbkey,species):
	"""
	returns a gd_genotype file format from csv file (saved in excel).
	The input must consist of a set of rows with columns splited by a
	comma. The first row must contain the names of the individuals. For
	the other rows, the first of column must contain the chromosome and
	position of the snp. The other columns must contain any kind of
	fstat or genepop allele valid encoding, with at	most 2 alleles. Also,
	the user may input IUPAC valid nucleotide symbols. The program will
	assume that the mosts common is the ancestral.
	------------- The file starts here ---------------- 
	,Ind_1,Ind_2,Ind_3,Ind_4,Ind_5,Ind_6,Ind_7
	chr1 12334123,A,T,T,A,T,T,W
	chr2 1232654,C,G,G,C,N,S,N
	chr3    3356367,T,G,G,G,G,K,K
	chr4    95673,A,C,C,A,C,M,M
	chr5 45896,T,C,Y,Y,Y,C,T
	...
	or
	...
	chr6 2354,22,11,21,00,12,12,22
	------------- The file ends here -------------------
	"""
	infoLn=True
	slf=open(outgd_gentp,'w')
	for echl in open(in_csv,'r'):
		if echl.strip():
			if infoLn:
				lPops=[x for x in echl.splitlines()[0].split(',') if x.strip()]
				hdr=formthdr(lPops,dbkey,species)
				slf.write('%s\n'%hdr)
				infoLn=False
			else:
				lsplitd=echl.splitlines()[0].split(',')
				lchrpox,SNPs=lsplitd[0],[x for x in lsplitd[1:] if x.strip()]
				lchrpox='\t'.join(lchrpox.strip().split())
				if SNPs[0].isdigit():
					lSNPsEncoded=[]
					for snp in SNPs:
						cnt=0
						for ep in snp:
							ep=int(ep)
							assert ep<=2
							cnt+=ep
						cnt=4-cnt
						if cnt==4:
							frmtdSNP='-1'
						else:
							frmtdSNP=str(cnt)
						lSNPsEncoded.append(frmtdSNP)
					ancd,dervd='N','N'
				else:
					ancd,dervd,lSNPsEncoded=selAnc(SNPs)
				outfrmtdDat='%s\t%s\t%s\t%s'%(lchrpox,ancd,dervd,'\t'.join(lSNPsEncoded))
				slf.write('%s\n'%outfrmtdDat)
	#~ 
	slf.close()
	return 0
				

def from_fstat(in_fstat,outgd_gentp,dbkey,species):
	"""
	returns a gd_genotype file format from fstat file. Ignores pops
	structures and alleles other than the combinations of the alleles
	encoded by 0, 1, and 2 up to 9 digits. The first line contains 4 
	numbers separated by any number of spaces: number of samples, np,
	number of loci, nl, highest number used to label an	allele, nu
	(?=2), and 1 if the code for alleles is a one digit	number (1-2), a
	2 if code for alleles is a 2 digit number (01-02) or a 3 if code for
	alleles is a 3 digit number (001-002). Followed by nl lines: each
	containing the name of a locus, in the order they will appear in the
	rest of the file On line nl+2, a series of numbers as follow: 1 0102
	0101 0101 0201 0 0101 first number: identifies the sample to which
	the individual belongs second: the genotype of the individual at the
	first locus, coded with a 2 digits number for each allele third: the
	genotype at the second locus, until locus nl is entered. Missing
	genotypes are encoded with 0 (0001 or 0100 are not a valid format,
	so both alleles at a locus have to be known, otherwise, the genotype
	is considered as missing) no empty lines are needed between samples
	number of spaces between genotypes can be anything. Numbering of
	samples need not be sequential the number of samples np needs to be
	the same as the largest sample identifier. Samples need not to be
	ordered nu needs to be equal to the largest code given to an allele
	(even if there are less than nu alleles). Ancestral is taken as 01,
	derived 02. In all cases ancestral and derived SNPs are returned as
	N.
	------------- The file starts here ---------------- 
	7  5  2  1
	chr1 12334123
	chr2 1232654
	chr3    3356367
	chr4    95673
	chr5 45896
	   Ind_1   22 22 21 11 22
	   Ind_2   22 22 11 12 22
	   Ind_3   22 11 22 21 22
	   Ind_4   22 21 22 21 22
	   Ind_5   22 22 22 21 22
	   Ind_6   22 22 22 22 22
	   Ind_7   22 22 21 21 22
	------------- The file ends here -------------------
	"""
	dChrPos,lPops,lChrPos,addPop={},[],[],False
	clines=-1
	for echl in open(in_fstat,'r'):
		clines+=1
		if echl.strip():
			if clines==0:
				nSmpls,nLoci,nUsed,nDigs=[x for x in echl.splitlines()[0].split() if x.strip()]
				nLoci=int(nLoci)
			elif clines<=nLoci:
				addPop=True
				lchrpox='\t'.join(echl.strip().split())
				lChrPos.append(lchrpox)
			elif addPop:
				lsplitd=echl.splitlines()[0].split()
				pop,SNPs=lsplitd[0],[x for x in lsplitd[1:] if x.strip()]
				pop=pop.strip()
				for x in range(nLoci):
					snp=SNPs[x]
					cnt=0
					for ep in snp:
						ep=int(ep)
						assert ep<=2
						cnt+=ep
					cnt=4-cnt
					if cnt==4:
						frmtdSNP='-1'
					else:
						frmtdSNP=str(cnt)
					try:
						dChrPos[lChrPos[x]].append(frmtdSNP)
					except:
						dChrPos[lChrPos[x]]=[frmtdSNP]
				#~ 
				lPops.append(pop)
	#~ 
	hdr=formthdr(lPops,dbkey,species)
	outfrmtdDat=['%s\t%s\t%s\t%s'%(x,'N','N','\t'.join(dChrPos[x])) for x in lChrPos]
	#~ 
	slf=open(outgd_gentp,'w')
	slf.write('\n'.join([hdr,'\n'.join(outfrmtdDat)]))
	slf.close()
	return 0
	

def from_genepop(in_genepop,outgd_gentp,dbkey,species):
	"""
	returns a gd_genotype file format from genepop file . Ignores pops
	structures and alleles other than the combinations of the alleles
	encoded by 00, 01, and 02. The second line must contain the chromosome
	and position of the SNPs separated by an space or a tab. Each loci
	should be separated by a comma. Alternatively, they may be given one
	per line. They may be given one per line, or on the same line but
	separated by commas. The name of individuals are defined as
	everything on the left of a comma, and their genotypes following the
	same order of the SNPs in the second line. Ancestral is taken as 01,
	derived 02. In all cases ancestral and derived SNPs are returned as N
		------------- The file starts here ---------------- 
	Microsat on Chiracus radioactivus, a pest species 
		 chr1 23123, chr2 90394, chr3 90909, chr3 910909, chr4 10909
	POP 
	AA2, 0201 0111 0102 0000      0101 
	AA1, 0201 0201 0202 0000      0101 
	A10, 0201 0201 0101 0000      0101 
	A11, 0201 0202 0102 0000      0102 
	A12, 0202 0201 0101 0000      0101 
	A11, 0101 0101 0101 0000      0101 
	A12, 0202 0201 0201 0000      0101 
	A11, 0201 0201 0101 0000      0101 
	Pop
	AF1, 0000 0000 0000 0000      0101 
	AF2, 0201 0101 0102 0000      0101 
	AF3, 0202 0201 0202 0000      0101 
	AF4, 0201 0101 0000 0000      0101 
	AF5, 0201 0101 0202 0000      0101 
	AF6, 0101 0101 0102 0000      0101 
	AF7, 0201 0100 0000 0000      0101 
	AF8, 0101 0100 0000 0000      0201 
	AF9, 0201 0200 0000 0000      0101 
	AF10, 0101 0202 0202 0000      0101 
	pop 
	C211, 0101 0202 0202 0000      0101 
	C211, 0101 0101 0202 0000      0101 
	C21E, 0101 0102 0202 0000      0101 
	C21B, 0101 0101 0102 0000      0201 
	C21C, 0201 0101 0202 0000      0101 
	C21D, 0201 0101 0202 0000      0201 
	------------- The file ends here -------------------
	"""
	dChrPos,lPops,lChrPos,addPop,addDat={},[],[],False,True
	clines=-1
	for echl in open(in_genepop,'r'):
		clines+=1
		if echl.strip():
			if echl.strip() in set(['pop','POP','Pop']):
				addDat,addPop=False,True
			elif addDat and clines>0:
				lchrpox=['\t'.join(x.split()) for x in echl.split(',') if x.strip()]
				lChrPos.extend(lchrpox)
			elif addPop:
				pop,SNPs=echl.splitlines()[0].split(',')
				pop=pop.strip()
				SNPs=[x for x in SNPs.split() if x.strip()]
				for x in range(len(SNPs)):
					snp=SNPs[x]
					cnt=0
					for ep in snp:
						ep=int(ep)
						assert ep<=2
						cnt+=ep
					cnt=4-cnt
					if cnt==4:
						frmtdSNP='-1'
					else:
						frmtdSNP=str(cnt)
					try:
						dChrPos[lChrPos[x]].append(frmtdSNP)
					except:
						dChrPos[lChrPos[x]]=[frmtdSNP]
				#~ 
				lPops.append(pop)
	#~ 
	hdr=formthdr(lPops,dbkey,species)
	outfrmtdDat=['%s\t%s\t%s\t%s'%(x,'N','N','\t'.join(dChrPos[x])) for x in lChrPos]
	#~ 
	slf=open(outgd_gentp,'w')
	slf.write('\n'.join([hdr,'\n'.join(outfrmtdDat)]))
	slf.close()
	return 0
	

def from_vcf(inVCF,outgd_gentp,dbkey,species):
	"""
	returns a gd_genotype file format from vcf a file
	"""
	slf=open(outgd_gentp,'w')
	paramSet,addPop,adinfo=False,False,False
	lPops=[]
	for echl in open(inVCF,'r'):
		if echl.strip():
			if not paramSet:
				if echl.find('##')==0:
					pass
				elif echl.find('#')==0:
					paramSet={}
					all_params=[x for x in echl.split() if x.strip()]
					clmn=-1
					for eparam in all_params:
						clmn+=1
						if eparam=='#CHROM':
							paramSet['chr']=clmn
						elif eparam=='POS':
							paramSet['pos']=clmn
						elif eparam=='REF':
							paramSet['A']=clmn
						elif eparam=='ALT':
							paramSet['B']=clmn
						elif eparam=='QUAL':
							paramSet['qual']=clmn
						elif eparam=='FORMAT':
							paramSet['frmt']=clmn
							addPop=True
						elif addPop:
							lPops.append(eparam)
							paramSet[eparam]=clmn
					if paramSet:
						hdr=formthdr(lPops,dbkey,species)
						slf.write('%s\n'%hdr)
			else:
				all_vals=[x for x in echl.split() if x.strip()]
				frmt=[x for x in all_vals[paramSet['frmt']].split(':') if x.strip()]
				clmn=-1
				gtclmn,adclmn,qulclmn=False,False,False
				for p in frmt:
					clmn+=1
					if p=='GT':
						gtclmn=clmn
					elif p=='AD':
						adclmn=clmn
						adinfo=True
					elif p=='GQ':
						qulclmn=clmn
				#~
				if adinfo:
					outptInfo=[all_vals[paramSet['chr']],all_vals[paramSet['pos']],all_vals[paramSet['A']],all_vals[paramSet['B']],all_vals[paramSet['qual']]]
					for ePop in lPops:
						gntyp=all_vals[paramSet[ePop]].replace('|','/').split(':')[gtclmn].split('/')
						encdGntyp,adA,adB,qual='-1','0','0','-1'
						#~ 	
						if set(gntyp)!=set(['.']):
							encdGntyp=2-sum([int(x) for x in gntyp])
							if adclmn:
								try:
									adA,adB=all_vals[paramSet[ePop]].split(':')[adclmn].split(',')
								except:
									pass
							if qulclmn:
								try:
									qual=all_vals[paramSet[ePop]].split(':')[qulclmn]
								except:
									pass
						outptInfo.extend([adA,adB,str(encdGntyp),qual])
					slf.write('%s\n'%'\t'.join(outptInfo))
				#~ 
				else:
					outptInfo=[all_vals[paramSet['chr']],all_vals[paramSet['pos']],all_vals[paramSet['A']],all_vals[paramSet['B']]]
					for ePop in lPops:
						gntyp=all_vals[paramSet[ePop]].replace('|','/').split(':')[gtclmn].split('/')
						try:
							encdGntyp=2-sum([int(x) for x in gntyp])
						except:
							encdGntyp=-1
						outptInfo.append(str(encdGntyp))
					#~ 
					slf.write('%s\n'%'\t'.join(outptInfo))
	slf.close()
	#~ 
	if adinfo:
		hdr=formthdr_gdsnp(lPops,dbkey,species)
		slf=open('%stmp'%outgd_gentp,'w')
		slf.write('%s\n'%hdr)
		appnd=False
		for echl in open(outgd_gentp,'r'):
			if appnd:
				slf.write(echl)
			else:
				if echl[0]!='#':
					appnd=True
					slf.write(echl)
		slf.close()
		#~ 
		os.system('mv %stmp %s'%(outgd_gentp,outgd_gentp))
	#~ 
	return 0	
				

def main():
	#~ 
	parser = argparse.ArgumentParser(description='Returns the count of genes in KEGG categories and their statistical overrrepresentation, from a list of genes and an background file (i.e. plane text with ENSEMBLT and KEGG pathways).')
	parser.add_argument('--input',metavar='input TXT file',type=str,help='the input file with the table in VCF format.',required=True)
	parser.add_argument('--output',metavar='output TXT file',type=str,help='the output file with the table in gd_genotype format.',required=True)
	parser.add_argument('--dbkey',metavar='string',type=str,help='the input reference species dbkey (i.e. susScr3).',required=True)
	parser.add_argument('--species',metavar='string',type=str,help='the input reference species name (i.e. int).',required=True)
	parser.add_argument('--format',metavar='string',type=str,help='format of the input file (i.e. vcf, genepop, fstat, csv).',required=True)

	args = parser.parse_args()

	infile = args.input
	outgd_gentp = args.output
	dbkey = args.dbkey
	species = base64.b64decode(args.species)
	frmat = args.format

	#~ 
	if frmat=='vcf':
		from_vcf(infile,outgd_gentp,dbkey,species)
	elif frmat=='genepop':
		from_genepop(infile,outgd_gentp,dbkey,species)
	elif frmat=='fstat':
		from_fstat(infile,outgd_gentp,dbkey,species)
	elif frmat=='csv':
		from_csv(infile,outgd_gentp,dbkey,species)

	#~ 
	return 0

if __name__ == '__main__':
	main()

