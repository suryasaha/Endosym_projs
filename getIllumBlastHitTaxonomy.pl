#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Apr 17, 2012

use strict;
use warnings;
use Getopt::Long;
eval {
	require Bio::SearchIO;
};
use Bio::SearchIO; 
use Switch;
use lib '/home/surya/bin/modules';
use NCBI_Taxonomy;
use NCBI;

=head1 NAME

 getIllumBlastHitTaxonomy.pl - Write out report of Taxonomy of Blastn hits

=head1 SYNOPSIS

  % getIllumBlastHitTaxonomy.pl --minid 0.9 --in blastn.out --out blastn.taxonomy.txt

=head1 DESCRIPTION

 Produce a matrix of taxonomy hierarchy with counts for each remote blast hit for a illum 
 region. Need hits to have Accessions or GI numbers in 2nd position before space in name in 
 report. A very specific purpose to this tool. May not apply for general taxonomy analysis.
 
=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. 
 Some options are mandatory (see below).

   --minid  <0.9>   Minimum percent identity shared between query and hit for illuminated 
                    region without gaps. A float value (0.0-1.0, required)
   --in     <.out>  NCBI blastn generated report (required)
   --out    <.txt>  Name of output report

=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my($mid,$rep,$out,$in,$i,$j,$k,@temp,$result,$hit,$hsp);

GetOptions (
	'minid=f' => \$mid,
	'in=s' => \$rep,
	'out:s' => \$out) or (system('pod2text',$0), exit 1);
defined($mid) or (system('pod2text',$0), exit 1);
if (($mid<0.0)||($mid>1.0)){print STDERR "Minimum percent identity shared between query and hit is invalid. Need to be between 0.0 to 1.0.\n";}
defined($rep) or (system('pod2text',$0), exit 1);
if (!(-e $rep)){print STDERR "$rep not found: $!\n"; exit 1;}
$out ||= "$rep\.taxonomy\.txt";
unless(open(OUT,">$out")){print "not able to open $out\n\n";exit 1;}

$in = new Bio::SearchIO(-format => 'blast', -file   => $rep);

while($result = $in->next_result){
	## $result is a Bio::Search::Result::ResultI compliant object
	if($result->no_hits_found()){next;}
	
	#get blast report
	if($result->num_hits>0){
		my($qhspident,$gi,$tid,@tax,%phy,%cls,%ord,%fam,%gen,%spc,$hcnt,$qhcnt);
		$hcnt=$qhcnt=0;
		while($hit = $result->next_hit){
			#get total length of all hsps
			$qhspident=0;
			while($hsp = $hit->next_hsp()){
				$qhspident+=$hsp->num_identical;#illum region
			}
			$i=$qhspident/$result->query_length();
			if($i >= $mid){#if frac ident qualifies
				#get taxonomy info and print it
				@temp=split(/\|/,$hit->name);
				$gi=&NCBI::Accession2GI($temp[1]); chomp $gi;
				#@temp=&NCBI_Taxonomy::GI2TaxonomyLevels($gi);
				#BUG What if no taxonomy record exists for GI??
				$tid=&NCBI_Taxonomy::GI2TaxID($gi); chomp $tid;
				@tax=&NCBI_Taxonomy::TaxID2TaxonomyLevels($tid);
				for $j (@tax){
					switch($j){
						case {$j->[1] eq 'phylum'} {
							if (exists $phy{$j->[0]}) {$phy{$j->[0]}++}
							else{ $phy{$j->[0]}=1;}
						}
						case {$j->[1] eq 'class'} {
							if (exists $cls{$j->[0]}) {$cls{$j->[0]}++}
							else{ $cls{$j->[0]}=1;}
						}
						case {$j->[1] eq 'order'} {
							if (exists $ord{$j->[0]}) {$ord{$j->[0]}++}
							else{ $ord{$j->[0]}=1;}
						}
						case {$j->[1] eq 'family'} {
							if (exists $fam{$j->[0]}) {$fam{$j->[0]}++}
							else{ $fam{$j->[0]}=1;}
						}
						case {$j->[1] eq 'genus'} {
							if (exists $gen{$j->[0]}) {$gen{$j->[0]}++}
							else{ $gen{$j->[0]}=1;}
						}
						case {$j->[1] eq 'species'} {
							if (exists $spc{$j->[0]}) {$spc{$j->[0]}++}
							else{ $spc{$j->[0]}=1;}
						}
					}
				}
				$qhcnt++;
			}
			$hcnt++;
		}
		print OUT "Query\t",$result->query_name,"\nLength\t",$result->query_length,"\nHits\t",$hcnt,"\nQual Hits\t",$qhcnt,"\n";
		print OUT "Phylum"; while (($j,$k) = each %phy){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "Class"; while (($j,$k) = each %cls){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "Order"; while (($j,$k) = each %ord){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "Family"; while (($j,$k) = each %fam){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "Genus"; while (($j,$k) = each %gen){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "Species"; while (($j,$k) = each %spc){print OUT "\t$j\t$k";} print OUT "\n";
		print OUT "\n\n";
	}
}
close (OUT);
exit;
