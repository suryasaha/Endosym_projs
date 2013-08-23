#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 9, 2012

use lib '/home/surya/bin/modules';
use SS;

use strict;
use warnings;
use Getopt::Long;


=head1 NAME

 illumgff2MFasta.pl - Write out fasta sequences for screened? regions listed in GFF

=head1 SYNOPSIS

  % illumgff2MFasta.pl --fna file1.fna --gff file2.gff --scr avoid.gff --minlen 500 --out out.fna
  
=head1 DESCRIPTION

 To get MFasta of regions of reference sequence illuminated by Bowtie2 mapping after screening?.
 
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --fna    <.fna>  Genome sequence in fasta format (required)
   --gff    <.gff>  Manually(using Artemis??) or automaticaly generaged GFF file (required)
   --scr    <.gff>  Screening GFF file of ribosomal genes, if available
   --minlen <300>   Minimum length of illuminated region, if screening
   --out    <.fna>  Name of output mfasta file

=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my($fna,$gff,$scr,$out,$mlen,@maskedseq,$seq,$i,$j,$rec,@temp,$ctr,$regs);

GetOptions (
	'fna=s' => \$fna,
	'gff=s' => \$gff,
	'scr:s' => \$scr,
	'minlen:i' => \$mlen,
	'out:s' => \$out) or (system('pod2text',$0), exit 1);

defined($fna) or (system('pod2text',$0), exit 1);
if (!(-e $fna)){print STDERR "$fna not found: $!\n"; exit 1;}
defined($gff) or (system('pod2text',$0), exit 1);
if (!(-e $gff)){print STDERR "$gff not found: $!\n"; exit 1;}
if ((defined $scr) && (!(-e $scr))){print STDERR "$scr not found: $!\n"; exit 1;}
if(defined $scr){$mlen ||= 300;}#non-ovlapping length > 300 bp (ML commn.)
$i=$gff; $i=~ s/.gff$//;
if(defined $scr){$out ||= "$i\.screened.mfasta";}
else{$out ||= "$i\.mfasta";}
unless(open(INGFF,"<$gff")){print "not able to open $gff\n\n";exit 1;}
if(defined $scr){ unless(open(INSCR,"<$scr")){print "not able to open $scr\n\n";exit 1;}}
unless(open(INFAS,"<$fna")){print "not able to open $fna\n\n";exit 1;}
unless(open(OUT,">$out")){print "not able to open $out\n\n";exit 1;}

# read in files
$seq='';
while($rec=<INFAS>){
	if($rec =~ />/){next;}
	$rec=~ s/\s*//g;#clean
	$seq=$seq.$rec;
}
close(INFAS);

for $i (1..length $seq){ $maskedseq[$i]=1;}#note position 0 is undefined
# get screening coords
if(defined $scr){
	while($rec=<INSCR>){
		if($rec =~ /#/){next;}
		@temp=split("\t",$rec);
		for $i($temp[3]..$temp[4]){ $maskedseq[$i]=0;}
	}
	close(INSCR);
}

#gff_seqname	artemis	exon	12136	12720	.	+	.	ID=CDS:12136..12720
#gff_seqname	artemis	exon	69917	70358	.	+	.	ID=CDS:69917..70358
#gff_seqname	artemis	exon	71700	72604	.	+	.	ID=CDS:71700..72604
# read in file
$regs=$ctr=0;
while($rec=<INGFF>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec);

	if(defined $scr){#check for overlap with ribosomal genes
		my ($st,$end);
		$j=0; chomp $temp[8];
		#find unmasked illuminated regions, print 1st regions if > mlen, go on till end
		for $i($temp[3]..$temp[4]){
			if(($maskedseq[$i]==0) && ($j==0)){#masked region
				$st=$end=0;
			}
			#elsif(($maskedseq[$i]==0) && ($j>0)){#masked region after unmasked
			elsif(($maskedseq[$i]==0) && (($j)>=$mlen)){#masked region after unmasked
				if($temp[8] ne ''){ print OUT "\>$temp[8]\|pos $st-$end \|length $j\n";}
				else{ print OUT "\>$ctr\|pos $st-$end \|length $j\n";}
				if(($temp[6] eq '+') || ($temp[6] eq '.')){
					print OUT substr($seq,$st-1,(($end-1)-($st-1))+1),"\n";#coordinate space to index space
				}
				elsif($temp[6] eq '-'){
					print OUT &SS::revcomp(substr($seq,$st-1,(($end-1)-($st-1))+1)),"\n";#coordinate space to index space
				}
				$st=$end=$j=0;#reset
				$regs++;
			}
			elsif($maskedseq[$i]==1) {#increment end of illum region
				if($j==0){ $st=$end=$i;}#for new sub-region after a masked section
				elsif($j>0){ 
					$end=$i;
					if(($end==$temp[4]) && (($j++)>$mlen)){
						if($temp[8] ne ''){ print OUT "\>$temp[8]\|pos $st-$end \|length $j\n";}
						else{ print OUT "\>$ctr\|pos $st-$end \|length $j\n";}
						if(($temp[6] eq '+') || ($temp[6] eq '.')){
							print OUT substr($seq,$st-1,(($end-1)-($st-1))+1),"\n";#coordinate space to index space
						}
						elsif($temp[6] eq '-'){
							print OUT &SS::revcomp(substr($seq,$st-1,(($end-1)-($st-1))+1)),"\n";#coordinate space to index space
						}
						$regs++;						
					}
				}
				$j++;
			}
		}
	}
	#else{#if no screening file
	elsif(($temp[4]-$temp[3]+1)>=$mlen){#if no screening file
		if($temp[8] ne ''){
			chomp $temp[8]; 
			print OUT "\>$temp[8]\|pos $temp[3]-$temp[4] \|length ",$temp[4]-$temp[3]+1,"\n";
		}
		else{ print OUT "\>$ctr:pos $temp[3]-$temp[4] \|length ",$temp[4]-$temp[3]+1,"\n";}

		if(($temp[6] eq '+') || ($temp[6] eq '.')){
			print OUT substr($seq,$temp[3]-1,(($temp[4]-1)-($temp[3]-1))+1),"\n";#coordinate space to index space
		}
		elsif($temp[6] eq '-'){
			print OUT &SS::revcomp(substr($seq,$temp[3]-1,(($temp[4]-1)-($temp[3]-1))+1)),"\n";#coordinate space to index space
		}
		$regs++;		
	}
	$ctr++;
}
close(INGFF); close(OUT);
print STDERR "$ctr GFF records processed\n";
print STDERR "$regs illuminated regions written to Fasta\n";

exit;