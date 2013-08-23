#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 9, 2012

use strict;
use warnings;
use Getopt::Long;
use Switch;

=head1 NAME

 covgff2illumgff.pl - Write out GFF file of illuminated regions

=head1 SYNOPSIS

  % covgff2illumgff.pl --len 300 --cov 1 --gap% 5.0 --in file2.gff --out out.illum.gff

=head1 DESCRIPTION

 To get GFF of regions of reference sequence illuminated by Bowtie2 mapping. Combining 
 regions that are separated by short distances.
 
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --len    <300>   Minimum length of illuminated region (required)
   --cov    <1>     Minimum coverage of illuminated region (required)
   --gap%   <5.0>   Maximum cumulative length of gaps in % of total length within a illuminated region (required)
   --in     <.gff>  genomeCoverageBed2GFF3.pl generated GFF file (required)
   --out    <.gff>  Name of output illuminated GFF file

=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my($mcov,$mgap,$gff,$out,$mlen,@maskedseq,$seq,$i,$j,$rec,@temp,$ctr,$st,$end,%ireg);

GetOptions (
	'len=i' => \$mlen,
	'cov=i' => \$mcov,
	'gap=f' => \$mgap,
	'in=s' => \$gff,
	'out:s' => \$out) or (system('pod2text',$0), exit 1);

defined($mlen) or (system('pod2text',$0), exit 1);
defined($mcov) or (system('pod2text',$0), exit 1);
defined($mgap) or (system('pod2text',$0), exit 1);
defined($gff) or (system('pod2text',$0), exit 1);
if (!(-e $gff)){print STDERR "$gff not found: $!\n"; exit 1;}
$i=$gff; $i=~ s/.gff$//;
$out ||= "$i\.illum\.gff";

unless(open(INGFF,"<$gff")){print "not able to open $gff\n\n";exit 1;}
#if(defined $scr){ unless(open(INSCR,"<$scr")){print "not able to open $scr\n\n";exit 1;}}
#unless(open(INFAS,"<$fna")){print "not able to open $fna\n\n";exit 1;}
unless(open(OUT,">$out")){print "not able to open $out\n\n";exit 1;}
print OUT "\#\#gff-version 3\n\#Min length: $mlen\n\#Min coverage: $mcov\n\#Max gap\%: $mgap\n";

#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	204	213	3	.	.	colour=9;
#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	214	214	5	.	.	colour=4;
#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	215	242	6	.	.	colour=4;
#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	243	257	9	.	.	colour=4;
#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	258	301	11	.	.	colour=16;
#gi|224466365|gb|CP000857.1|	BEDTools	misc_feature	302	306	10	.	.	colour=16;
# read in file
$ctr=$j=0;
$ireg{'st'}=$ireg{'end'}=$ireg{'1X'}=$ireg{'3X'}=0;
$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=$ireg{'gap'}=0;
$ireg{'segments'}=0;   
while($rec=<INGFF>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec);

	if(($ctr==0) && ($temp[5]>=$mcov)){#first record in file
		$ireg{'st'}=$temp[3];  $ireg{'end'}=$temp[4];
		switch($temp[5]){
			case {$_[0] < 3} {$ireg{'1X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 5} {$ireg{'3X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 10} {$ireg{'5X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 20} {$ireg{'10X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] >= 20} {$ireg{'20X'}=($temp[4]-$temp[3]+1)}
		}
		$ireg{'segments'}++;
	}
	elsif(($ireg{'st'}==0) && ($temp[5]>=$mcov)){#first record after a low coverage record where all vals reset to 0
		$ireg{'st'}=$temp[3];  $ireg{'end'}=$temp[4];
		switch($temp[5]){
			case {$_[0] < 3} {$ireg{'1X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 5} {$ireg{'3X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 10} {$ireg{'5X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] < 20} {$ireg{'10X'}=($temp[4]-$temp[3]+1)}
			case {$_[0] >= 20} {$ireg{'20X'}=($temp[4]-$temp[3]+1)}
		}
		$ireg{'segments'}++;
	}
	else{
		#(contiguous region or within gap tolerance) and above reqd coverage, just record and next
		#total length for gap tolerance calculated for all segments uptill this record, eligible following
		#records cannot be considered
		if(((($ireg{'end'}+1)==$temp[3]) || 
			(($temp[3]-$ireg{'end'}-1+$ireg{'gap'})<=(int(($temp[4]-$ireg{'st'}+1)*($mgap/100)))))
			&& ($temp[5] >= $mcov)){
#		if(((($ireg{'end'}+1)==$temp[3]) ||	(($temp[3]-$ireg{'end'}-1)<=5))	&& ($temp[5] >= $mcov)){#DEBUGGING	
			#increment gap, if any
			if(($ireg{'end'}+1)!=$temp[3]){ $ireg{'gap'}+=$temp[3]-$ireg{'end'}-1;}
			$ireg{'end'}=$temp[4];
			switch($temp[5]){
				case {$_[0] < 3} {$ireg{'1X'}+=($temp[4]-$temp[3]+1)}
				case {$_[0] < 5} {$ireg{'3X'}+=($temp[4]-$temp[3]+1)}
				case {$_[0] < 10} {$ireg{'5X'}+=($temp[4]-$temp[3]+1)}
				case {$_[0] < 20} {$ireg{'10X'}+=($temp[4]-$temp[3]+1)}
				case {$_[0] >= 20} {$ireg{'20X'}+=($temp[4]-$temp[3]+1)}
			}
			$ireg{'segments'}++;
		}
		#non-contigous region or disqualified contiguous region, write if length > mlen and reset
		elsif(($ireg{'end'}-$ireg{'st'}+1) >= $mlen){
			print OUT "$temp[0]\tBEDTools\tmisc_feature\t$ireg{'st'}\t$ireg{'end'}\t";
			#score (max 20xlen, min 1xlen, so scale from 1 to 20)
			$i=(($ireg{'1X'}*1) + ($ireg{'3X'}*3) + ($ireg{'5X'}*5) + ($ireg{'10X'}*10) + ($ireg{'20X'}*20));
			$i/=($ireg{'end'}-$ireg{'st'}+1);
#			#scale from 1 to 100, SOMETHING WRONG WITH MATH
#			#http://stackoverflow.com/questions/3460905/scale-a-list-of-numbers-to-be-between-1-0-and-1-0
#			$i = sprintf("%.2f",(((99*$i-(($ireg{'end'}-$ireg{'st'}+1)*1)) 
#				/ (($ireg{'end'}-$ireg{'st'}+1)*20) - (($ireg{'end'}-$ireg{'st'}+1)*1))) + 1);
			print OUT sprintf("%.2f",$i),"\t.\t.\tlength\=",($ireg{'end'}-$ireg{'st'}+1),"\; segments\=",$ireg{'segments'},
				"\; gaps\=",$ireg{'gap'},"\; 1\-2X\=",$ireg{'1X'},"\; 3\-4X\=",$ireg{'3X'},"\; 5\-9X\=",
				$ireg{'5X'},"\; 10\-19X\=",$ireg{'10X'},"\; 20\+X\=",$ireg{'20X'},"\;\n";
			$j++;
			#reset
			if ($temp[5]>=$mcov){
				$ireg{'st'}=$temp[3]; $ireg{'end'}=$temp[4];
				$ireg{'gap'}=0;	$ireg{'segments'}=1;
				switch($temp[5]){
					case {$_[0] < 3} {$ireg{'1X'}=($temp[4]-$temp[3]+1); $ireg{'3X'}=$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 5} {$ireg{'3X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 10} {$ireg{'5X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 20} {$ireg{'10X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'5X'}=$ireg{'20X'}=0;}
					case {$_[0] >= 20} {$ireg{'20X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'5X'}=$ireg{'10X'}=0;}
				}
			}
			else{#reinit to null for next eligible region
				$ireg{'st'}=$ireg{'end'}=$ireg{'1X'}=$ireg{'3X'}=0;
				$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=$ireg{'gap'}=0;
				$ireg{'segments'}=0; 
			}
		}
		#last region was dud but current record is eligible
		elsif(($temp[5]>=$mcov) && (($ireg{'end'}-$ireg{'st'}+1) < $mlen)){
				$ireg{'st'}=$temp[3]; $ireg{'end'}=$temp[4];
				$ireg{'gap'}=0;	$ireg{'segments'}=1;
				switch($temp[5]){
					case {$_[0] < 3} {$ireg{'1X'}=($temp[4]-$temp[3]+1); $ireg{'3X'}=$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 5} {$ireg{'3X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 10} {$ireg{'5X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'10X'}=$ireg{'20X'}=0;}
					case {$_[0] < 20} {$ireg{'10X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'5X'}=$ireg{'20X'}=0;}
					case {$_[0] >= 20} {$ireg{'20X'}=($temp[4]-$temp[3]+1); $ireg{'1X'}=$ireg{'3X'}=$ireg{'5X'}=$ireg{'10X'}=0;}
				}
		}
		else{#reinit to null for next eligible region
			$ireg{'st'}=$ireg{'end'}=$ireg{'1X'}=$ireg{'3X'}=0;
			$ireg{'5X'}=$ireg{'10X'}=$ireg{'20X'}=$ireg{'gap'}=0;
			$ireg{'segments'}=0;
		}
	}
	$ctr++;	
}
#write last illum region
if(($ireg{'end'}-$ireg{'st'}+1) >= $mlen){
	print OUT "$temp[0]\tBEDTools\tmisc_feature\t$ireg{'st'}\t$ireg{'end'}\t";
	#score (max 20xlen, min 1xlen, so scale from 1 to 20)
	$i=(($ireg{'1X'}*1) + ($ireg{'3X'}*3) + ($ireg{'5X'}*5) + ($ireg{'10X'}*10) + ($ireg{'20X'}*20));
	$i/=($ireg{'end'}-$ireg{'st'}+1);
#	#scale from 1 to 100, SOMETHING WRONG WITH MATH
#	#http://stackoverflow.com/questions/3460905/scale-a-list-of-numbers-to-be-between-1-0-and-1-0
#	$i = sprintf("%.2f",(((99*$i-(($ireg{'end'}-$ireg{'st'}+1)*1)) 
#		/ (($ireg{'end'}-$ireg{'st'}+1)*20) - (($ireg{'end'}-$ireg{'st'}+1)*1))) + 1);
	print OUT sprintf("%.2f",$i),"\t.\t.\tlength\=",($ireg{'end'}-$ireg{'st'}+1),"\; segments\=",$ireg{'segments'},
		"\; gaps\=",$ireg{'gap'},"\; 1\-2X\=",$ireg{'1X'},"\; 3\-4X\=",$ireg{'3X'},"\; 5\-9X\=",
		$ireg{'5X'},"\; 10\-19X\=",$ireg{'10X'},"\; 20\+X\=",$ireg{'20X'},"\;\n";
	$j++;
}

close(INGFF); close(OUT);
print STDERR "$ctr records processed\n";
if($j>0){print STDERR "$j region(s) written to GFF\n";}
else{print STDERR "No qualifying illuminated regions found.. empty GFF\n";}

exit;