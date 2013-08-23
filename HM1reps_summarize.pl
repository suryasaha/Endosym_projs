#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha  

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Proc::ProcessTable;

sub memory_usage {
	my $t = new Proc::ProcessTable;
	foreach my $got ( @{$t->table} ) {
		next if not $got->pid eq $$;
		return $got->size;
	}
}

my ($i,$j,$k,@rec,$row,%mhits,%mcov);

#Gene name, Maize line, ??, query_name, subj_name, query_length, subj_length, 
# 0              1       2       3        4               5			6
#hsp_rank, identity, positive, match_length, score, p_value, query_start, 
#    7        8       9             10     		11    12		13
#query_end, query_strand, subj_start, subj_end, subj_strand, query_seq, subj_seq, alignment 
# 14             15     	  16       17           18   	     19       20		21
#Microbial 

unless (@ARGV == 5){
	print "USAGE: $0 <Input .TAB> <Output xls> <Tag aln%> <Positive aln%> <Ovlap %>\n"; 
	print "Make sure no headers are present. Use 0.9 for 90%\n"; exit;
}
chomp $ARGV[4];
unless(open(REP,"<$ARGV[0]")){print "not able to open $ARGV[0]\n\n";exit;}

# parse thru file, add to hash->=array, 2 for W22, B73, Mo17 each in that order
while($row=<REP>){
	if($row =~ /#/){next;}
	@rec=split("\t",$row);
	for $i (0..20){$rec[$i]=~ s/\s*//g;}
	
	if(exists $mhits{$rec[3]}){
		#90%+ identity of match length and 90%+ of Maize tag length
		if(($rec[1] eq "W22") && ($rec[8]>=(0.9*$rec[10])) && ($rec[10]>=(0.9*$rec[6]))){
			#calculate overlap
			$i=0; for $j($rec[13]..$rec[14]){ if($mcov{$rec[3]}[0][$j]==1){$i++;}} 
			#incr ctrs and mark overlap
			$mhits{$rec[3]}[0]++;
			if($i>=($ARGV[4]*$rec[10])){ $mhits{$rec[3]}[1]++;}#ovlap if 50% or more
			for $i($rec[13]..$rec[14]){ $mcov{$rec[3]}[0][$i]=1;}
		}
		elsif(($rec[1] eq "B73") && ($rec[8]>=(0.9*$rec[10])) && ($rec[10]>=(0.9*$rec[6]))){
			#calculate overlap
			$i=0; for $j($rec[13]..$rec[14]){ if($mcov{$rec[3]}[1][$j]==1){$i++;}} 
			#incr ctrs and mark overlap
			$mhits{$rec[3]}[2]++;
			if($i>=($ARGV[4]*$rec[10])){ $mhits{$rec[3]}[3]++;}#ovlap if 50% or more
			for $i($rec[13]..$rec[14]){ $mcov{$rec[3]}[1][$i]=1;}
		}
		elsif(($rec[1] eq "Mo17") && ($rec[8]>=(0.9*$rec[10])) && ($rec[10]>=(0.9*$rec[6]))){
			#calculate overlap
			$i=0; for $j($rec[13]..$rec[14]){ if($mcov{$rec[3]}[2][$j]==1){$i++;}} 
			#incr ctrs and mark overlap
			$mhits{$rec[3]}[4]++;
			if($i>=($ARGV[4]*$rec[10])){ $mhits{$rec[3]}[5]++;}#ovlap if 50% or more
			for $i($rec[13]..$rec[14]){ $mcov{$rec[3]}[2][$i]=1;}
		}
	}
	else{
		#init all 3x2 ctrs, W22, B73, Mo17
		for $i (0..5){ $mhits{$rec[3]}[$i]=0;}
		#init all 3 arr in hash->2D arr, W22, B73, Mo17
		$mcov{$rec[3]}[0][0]=$rec[5];#length in 0 location for W22
		$mcov{$rec[3]}[1][0]=$rec[5];#B73
		$mcov{$rec[3]}[2][0]=$rec[5];#Mo17
		
		for $i (1..$rec[5]){ $mcov{$rec[3]}[0][$i]=0;}#NO coord to array space translation
		for $i (1..$rec[5]){ $mcov{$rec[3]}[1][$i]=0;}
		for $i (1..$rec[5]){ $mcov{$rec[3]}[2][$i]=0;}
		
		#store values if 90%+ identity of match length and 90%+ of Maize tag length
		if(($rec[1] eq "W22") && ($rec[8]>=($ARGV[3]*$rec[10])) && ($rec[10]>=($ARGV[2]*$rec[6]))){
			for $i ($rec[13]..$rec[14]){ $mcov{$rec[3]}[0][$i]=1;}#mark coverage
			$mhits{$rec[3]}[0]++;#incr hit ctr
		}
		elsif(($rec[1] eq "B73") && ($rec[8]>=($ARGV[3]*$rec[10])) && ($rec[10]>=($ARGV[2]*$rec[6]))){
			for $i ($rec[13]..$rec[14]){ $mcov{$rec[3]}[1][$i]=1;}#mark coverage
			$mhits{$rec[3]}[2]++;#incr hit ctr
		}
		elsif(($rec[1] eq "Mo17") && ($rec[8]>=($ARGV[3]*$rec[10])) && ($rec[10]>=($ARGV[2]*$rec[6]))){
			for $i ($rec[13]..$rec[14]){ $mcov{$rec[3]}[2][$i]=1;}#mark coverage
			$mhits{$rec[3]}[4]++;#incr hit ctr
		}
	}
}
close(REP);

#count %mcov from 1 to end, NOT 0 to end for coverage
unless(open(OUT,">$ARGV[1]")){print "not able to open $ARGV[1]\n\n";exit;}
print OUT "\t\t\tMicrobial genes to HAPMAP II Blast Summary\n\n";
print OUT "Cutoffs\nMaize tag aln len\t$ARGV[2]\nIdentity\t$ARGV[3]\nOverlap\t$ARGV[4]\n\n";
print OUT "\t\t\t\t\tW22\t\t\tB73\t\t\tMo17\nMicrobial query\tLength\tHits\tOvlap Hits\tCov(%)";
print OUT "\tHits\tOvlap Hits\tCov(%)\tHits\tOvlap Hits\tCov(%)\tHits\tOvlap Hits\tCov(%)\n";
@rec=(); @rec=keys (%mhits);
foreach $i(@rec){
	print OUT "$i\t$mcov{$i}[0][0]\t";
	#calculate and print totals
	$j=$mhits{$i}[0]+$mhits{$i}[2]+$mhits{$i}[4];
	$k=$mhits{$i}[1]+$mhits{$i}[3]+$mhits{$i}[5];
	print OUT "$j\t$k\t";
	$k=0; 
	for $j(1..$mcov{$i}[0][0]){ if(($mcov{$i}[0][$j]==1) 
		|| ($mcov{$i}[1][$j]==1) || ($mcov{$i}[2][$j]==1)){$k++;}}
	$j=sprintf("%.2f",($k/$mcov{$i}[0][0])*100); print OUT "$j\t";
	print OUT "$mhits{$i}[0]\t$mhits{$i}[1]\t";#W22
	$k=0; for $j(1..$mcov{$i}[0][0]){ if($mcov{$i}[0][$j]==1){$k++;}}
	$j=sprintf("%.2f",($k/$mcov{$i}[0][0])*100); print OUT "$j\t";
	print OUT "$mhits{$i}[2]\t$mhits{$i}[3]\t";#B73
	$k=0; for $j(1..$mcov{$i}[1][0]){ if($mcov{$i}[1][$j]==1){$k++;}}
	$j=sprintf("%.2f",($k/$mcov{$i}[1][0])*100); print OUT "$j\t";
	print OUT "$mhits{$i}[4]\t$mhits{$i}[5]\t";#Mo17
	$k=0; for $j(1..$mcov{$i}[2][0]){ if($mcov{$i}[2][$j]==1){$k++;}}
	$j=sprintf("%.2f",($k/$mcov{$i}[2][0])*100); print OUT "$j\n";
}
close(OUT);

my($user_t,$system_t,$cuser_t,$csystem_t);	($user_t,$system_t,$cuser_t,$csystem_t) = times;
print STDERR "System time for process: $system_t\n"; print STDERR "User time for process: $user_t\n";
print STDERR "Process id=",$$,"\n"; print "Memory used(MB)=", memory_usage()/1024/1024, "\n";

exit;