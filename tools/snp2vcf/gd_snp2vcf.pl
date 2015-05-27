#!/usr/bin/perl -w
use strict;

#convert from gd_snp file to vcf file (with dbSNP fields)

#gd_snp table format:
#1. chr
#2. position (0 based)
#3. ref allele
#4. second allele
#5. overall quality
#foreach individual (6-9, 10-13, ...)
#a. count of allele in 3
#b. count of allele in 4
#c. genotype call (-1, or count of ref allele)
#d. quality of genotype call (quality of non-ref allele from masterVar)

if (!@ARGV) {
   print "usage: gd_snp2vcf.pl file.gd_snp[.gz|.bz2] -geno=8[,12:16,20...] -handle=HANDLE -batch=BATCHNAME -ref=REFERENCEID [-bioproj=XYZ -biosamp=ABC -population=POPID[,POPID2...] -chrCol=9 -posCol=9 ] > snpsForSubmission.vcf\n";
   exit;
}

my $in = shift @ARGV;
my $genoCols = '';
my $handle;
my $batch;
my $bioproj;
my $biosamp;
my $ref;
my $pop;
my $cr = 0; #allow to use alternate reference?
my $cp = 1;
my $meta;
my $offset = 0; #offset for genotype column, gd_snp vs gd_genotype indivs file
foreach (@ARGV) {
   if (/-geno=([0-9,]+)/) { $genoCols .= "$1:"; }
   elsif (/-geno=(.*)/) { $genoCols .= readGeno($1); }
   elsif (/-off=([0-9])/) { $offset = $1; }
   elsif (/-handle=(.*)/) { $handle = $1; }
   elsif (/-batch=(.*)/) { $batch = $1; }
   elsif (/-bioproj=(.*)/) { $bioproj = $1; }
   elsif (/-biosamp=(.*)/) { $biosamp = $1; }
   elsif (/-ref=(.*)/) { $ref = $1; } 
   elsif (/-population=(\S+)/) { $pop = $1; }
   elsif (/-chrCol=(\d+)/) { $cr = $1 - 1; }
   elsif (/-posCol=(\d+)/) { $cp = $1 - 1; }
   elsif (/-metaOut=(.*)/) { $meta = $1; }
}
if ($cr < 0 or $cp < 0) { die "ERROR the column numbers should be 1 based.\n"; }
 
#remove trailing delimiters
$genoCols =~ s/,:/:/g;
$genoCols =~ s/[,:]$//;

my @gnc = split(/,|:/, $genoCols);

if ($in =~ /.gz$/) { 
   open(FH, "zcat $in |") or die "Couldn't open $in, $!\n";
}elsif ($in =~ /.bz2$/) {
   open(FH, "bzcat $in |") or die "Couldn't open $in, $!\n";
}else {
   open(FH, $in) or die "Couldn't open $in, $!\n";
}
my @head = prepHeader();
if (@head) { 
   print join("\n", @head), "\n";
   #now column headers
   print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
   if (defined $pop) { 
      $pop =~ s/,$//;
      my $t = $pop;
      $t =~ s/,/\t/g;
      print "\tFORMAT\t$t";
   }
   print "\n";
}
while (<FH>) {
   chomp;
   if (/^#/) { next; }
   if (/^\s*$/) { next; } 
   my @f = split(/\t/);
   #vcf columns: chrom pos id ref alt qual filter info
   # info must have VRT=[0-9] 1==SNV 2=indel 6=NoVariation 8=MNV ...
   my $vrt = 1;
   if ($f[2] !~ /^[ACTG]$/ or $f[3] !~ /^[ACTG]$/) {
      die "Sorry this can only do SNV's at this time\n";
   }
   if (scalar @gnc == 1 && !defined $pop) { #single genotype column
      if (!defined $f[4] or $f[4] == -1) { $f[4] = '.'; }
      if ($f[$gnc[0]-1] == 2) { $vrt = 6; } #reference match
      if ($f[$gnc[0]-1] == -1) { next; } #no data, don't use
      print "$f[$cr]\t$f[$cp]\t$f[$cr];$f[$cp]\t$f[2]\t$f[3]\t$f[4]\t.\tVRT=$vrt\n";
      #TODO? put read counts in comment?
   }elsif ($pop) { #do as population
      my @cols;
      foreach my $gp (split(/:/,$genoCols)) { #foreach population
         my @g = split(/,/, $gp);
         my $totChrom = 2*(scalar @g);
         my $totRef = 0;
         foreach my $i (@g) { if (!defined $f[$i-1] or $f[$i-1] == -1) { $totChrom -= 2; next; } $totRef += $f[$i-1]; }
         if ($totChrom == $totRef) { $vrt = 6; }
         if ($totRef > $totChrom) { die "ERROR likely the wrong column was chosen for genotype\n"; }
         my $altCnt = $totChrom - $totRef;
         push(@cols, "$totChrom:$altCnt");
      }
      print "$f[$cr]\t$f[$cp]\t$f[$cr];$f[$cp]\t$f[2]\t$f[3]\t$f[4]\t.\tVRT=$vrt\tNA:AC\t", join("\t", @cols), "\n";
   }else { #leave allele counts off
      my $totChrom = 2*(scalar @gnc);
      my $totRef = 0;
      foreach my $i (@gnc) { if ($f[$i-1] == -1) { $totChrom -= 2; next; } $totRef += $f[$i-1]; }
      if ($totChrom == $totRef) { $vrt = 6; }
      print "$f[$cr]\t$f[$cp]\t$f[$cr];$f[$cp]\t$f[2]\t$f[3]\t$f[4]\t.\tVRT=$vrt\n";
   }
}
close FH or die "Couldn't close $in, $!\n";

if ($meta) {
   open(FH, ">", $meta) or die "Couldn't open $meta, $!\n";
   print FH "TYPE: CONT\n",
            "HANDLE: $handle\n",
            "NAME: \n",
            "FAX: \n",
            "TEL: \n",
            "EMAIL: \n",
            "LAB: \n",
            "INST: \n",
            "ADDR: \n",
            "||\n",
            "TYPE: METHOD\n",
            "HANDLE: $handle\n",
            "ID: \n",
            "METHOD_CLASS: Sequence\n",
            "TEMPLATE_TYPE: \n",
            "METHOD:\n",
            "||\n";
   if ($pop) {
      my @p = split(/,/, $pop);
      foreach my $t (@p) {
         print FH
            "TYPE: POPULATION\n",
            "HANDLE: $handle\n", 
            "ID: $t\n", 
            "POPULATION: \n", 
            "||\n";
      }
   }
   print FH "TYPE: SNPASSAY\n",
            "HANDLE: $handle\n",
            "BATCH: $batch\n",
            "MOLTYPE: \n",
            "METHOD: \n",
            "ORGANISM: \n",
            "||\n",
            "TYPE: SNPPOPUSE | SNPINDUSE\n",
            "HANDLE: $handle\n",
            "BATCH: \n",
            "METHOD: \n",
            "||\n";

   close FH or die "Couldn't close $meta, $!\n";
}

exit 0;

#parse old header and add or create new
sub prepHeader {
   my @h;
   $h[0] = '##fileformat=VCFv4.1';
   my ($day, $mo, $yr) = (localtime)[3,4,5];
   $mo++;
   $yr+=1900;
   $h[1] = '##fileDate=' . "$yr$mo$day";
   $h[2] = "##handle=$handle";
   $h[3] = "##batch=$batch";
   my $i = 4;
   if ($bioproj) { $h[$i] = "##bioproject_id=$bioproj";  $i++; }
   if ($biosamp) { $h[$i] = "##biosample_id=$biosamp"; $i++; }
   $h[$i] = "##reference=$ref";  ##reference=GCF_999999.99
   #$i++;
   #$h[$i] = '##INFO=<ID=LID, Number=1,Type=string, Description="Unique local variation ID or name for display. The LID provided here combined with the handle must be unique for a particular submitter.">'
   $i++;
   $h[$i] = '##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type,1 - SNV: single nucleotide variation,2 - DIV: deletion/insertion variation,3 - HETEROZYGOUS: variable, but undefined at nucleotide level,4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element,6 - NO VARIATON: sequence scanned for variation, but none observed,7 - MIXED: cluster contains submissions from 2 or more allelic classes (not used) ,8 - MNV: multiple nucleotide variation with all eles of common length greater than 1,9 - Exception">';
   #sometimes have allele freqs?
   if (defined $pop) {
      $i++;
      $h[$i] = "##FORMAT=<ID=NA,Number=1,Type=Integer,Description=\"Number of alleles for the population.\"";
      $i++;
      $h[$i] = '##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele count for each alternate allele.">';
      my @p = split(/,/, $pop);
      foreach my $t (@p) {
         $i++;
         $h[$i] = "##population_id=$t";
      }
   }
   #PMID?
##INFO=<ID=PMID,Number=.,Type=Integer,Description="PubMed ID linked to variation if available.">

   return @h;
}
####End

#read genotype columns from a file
sub readGeno { 
   my $list = shift @_;
   my @files = split(/,/, $list);
   my $cols='';
   foreach my $file (@files) {
      open(FH, $file) or die "Couldn't read $file, $!\n";
      while (<FH>) {
         chomp;
         my @f = split(/\s+/);
         if ($f[0] =~/\D/) { die "ERROR expect an integer for the column\n"; }
         $f[0] += $offset;
         $cols .= "$f[0],";
      }
      close FH;
      $cols .= ":";
   }
   $cols =~ s/,:$//;
   return $cols;
}
####End
