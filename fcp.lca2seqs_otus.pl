#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Dec 13, 2012

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use lib '/home/surya/bin/modules';
use SS;
eval {
	require Bio::DB::Fasta;
	require Switch;
};
use Switch;

=head1 NAME

fcp.lca2seqs_otus.pl - Convert FCP LCA output to seqs_otu format readable by QIIME make_otus.py script 

=head1 SYNOPSIS

  % orthomcl.getUniqProts.Pv.pl --lca lca.out --seq seq.fas --format fasta --tlevel genus 
  
=head1 DESCRIPTION

This script reads in output of FCP LCA 3 and outputs otu member lists in following format 
0 	seq1 	seq5 	 
1 	seq2 	  	 
2 	seq3 	  	 
3 	seq4 	seq6 	seq7
wt fasta files for all OTU members and representative members and also taxonomy file

=head1 Version

1. 

=head1 TODO

Convert all FCP output files based on cmd line param switch

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --lca    <.out>     Output of FCP LCA algorithm (required)
   --seq    <.fna>     Sequence file in Fasta, Genbank etc format (required)
   --format <x>        Format of sequence file. Should be readable by Bio::SeqIO. Default is fasta
   --tlevel <>         Taxonomy level to use for OTU selection (suggested: genus)
                       domain/phylum/class/order/family/genus/species/strain (required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

sub write_otu_rep_seq{
	# Params: otu number, seqname, sequence
	# 1, MICHAELJACKSON:1:1:14:85#0/1, TGTGTGGCCATTGGTCATCGACCAGAGGCTCATACAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG)
	# no return
	unless(open(OTUF,">>otu_rep_set.fas")){print "not able to open otu_rep_ses.fas\n\n";exit 1;}
	print OTUF '>',$_[0],' ',$_[1],"\n",$_[2],"\n"; close (OTUF);		
}

sub write_otu_seq{
	# Params: taxonomy level, otu name, open/append, seqname, sequence
	# class, 1, open/append, MICHAELJACKSON:1:1:14:85#0/1, TGTGTGGCCATTGGTCATCGACCAGAGGCTCATACAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG)
	# no return
	if($_[2] eq 'open'){
		unless(open(OTUF,">${_[0]}_otu_fasta/${_[1]}.fas")){print "not able to open otu_fasta/${_[1]}.fas\n\n";exit 1;}
		print OTUF '>',$_[3],"\n",$_[4],"\n"; close (OTUF);
	}
	elsif($_[2] eq 'append'){
		unless(open(OTUF,">>${_[0]}_otu_fasta/${_[1]}.fas")){print "not able to open otu_fasta/${_[1]}.fas\n\n";exit 1;}
		print OTUF '>',$_[3],"\n",$_[4],"\n"; close (OTUF);		
	}
	return;
}

my ($rep,$seq,$format,$tlevel);
GetOptions (
	'lca=s' => \$rep,
	'seq=s' => \$seq,
	'format:s' => \$format,
	'tlevel=s' => \$tlevel) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($rep) or (system('pod2text',$0), exit 1);
if (!(-e $rep)){print STDERR "$rep not found: $!\n"; exit 1;}
defined($seq) or (system('pod2text',$0), exit 1);
if (!(-e $seq)){print STDERR "$seq not found: $!\n"; exit 1;}
$format ||= 'fasta';
if(($tlevel ne 'domain') && ($tlevel ne 'phylum') && ($tlevel ne 'class') && ($tlevel ne 'order') && ($tlevel ne 'family')
&& ($tlevel ne 'genus') && ($tlevel ne 'species') && ($tlevel ne 'strain')){
	print STDERR "tlevel needs to be domain/phylum/class/order/family/genus/species/strain\n"; exit 1;
}
#$oname ||= 'Pvirid';

my ($i,$j,$k,@temp,@temp1,$rec,$ctr);

# whole section can be pushed into a function later to handle all FCP output reports 
# only input reqd is report and fasta name
unless(open(INREP,"${rep}")){print "not able to open ${rep}\n\n";exit 1;}
$rec=<INREP>;#header row

# LCA report
#PRESLEY_0001:1:73:867:1796#0/1	Bacteria;Cyanobacteria;Cyanobacteria (class);Chroococcales;Chroococcales (family);Synechococcus;Synec        hococcus sp. CC9605;Synechococcus sp. CC9605;
#PRESLEY_0001:1:73:603:642#0/2	unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified        ;
#MENDEL_0001:1:18:37:2024#0/1	Archaea;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;
#PRESLEY_0001:1:82:1458:1912#0/2unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassifi        ed;

# Taxonomy file
#334 PC.607_1128	Root;k__Bacteria;p__Tenericutes;c__Erysipelotrichi;o__Erysipelotrichales;f__Erysipelotrichaceae	1.000
#235 PC.634_190	Root;k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__	1.000
#225 PC.481_373	Root;k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__	1.000

my (%otuNames,%otuMembers,%otuTaxonomy,$seqName,$db,$nofMembers,$nofSeqs);
# data description
# otuNames: {name}=otu number
# otuMembers: {otu number}= seq1	seq2	seq3	....
# otuTaxonomy: {otu number}= domain;phylum;class;order;family;genus;species;strain 

#$seqFileObj = Bio::SeqIO->new('-file' => "<$seq",'-format' => $format );
# making DB of the seq file, forcing reindex
$db = Bio::DB::Fasta->new($seq, '-reindex'=>1);
if(!($db)){ die "These was a problem creating database from $seq\n";}

&SS::mk_dir("${tlevel}_otu_fasta");
$ctr=$nofMembers=$nofSeqs=0;#otu name ctr and others
while($rec=<INREP>){
	@temp=(); @temp=split("\t",$rec);
	$seqName=$temp[0];#name
	@temp1=split(';',$temp[1]);
	$nofSeqs++;
	# chk qualifying read, add to otu hashes, write fasta
	# skip for unclassified domain, since remaining levels will not be classified anyways
	if(($temp1[0] ne 'unclassified')){
		switch($tlevel){
			case {$_[0] eq 'domain'} {
				if (exists $otuNames{$temp1[0]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[0]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[0]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[0]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[0]}=$ctr;
					$otuMembers{$otuNames{$temp1[0]}}=$seqName; $otuTaxonomy{$otuNames{$temp1[0]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));#first seq is rep
					&write_otu_seq($tlevel, $otuNames{$temp1[0]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
			case {$_[0] eq 'phylum'}{
				$temp1[1]=~ s/ (phylum)$//;
				if (exists $otuNames{$temp1[1]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[1]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[1]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[1]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[1]}=$ctr;
					$otuMembers{$otuNames{$temp1[1]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[1]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[1]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}		
			case {$_[0] eq 'class'} { 
				$temp1[2]=~ s/ (class)$//;
				if (exists $otuNames{$temp1[2]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[2]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[2]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[2]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[2]}=$ctr;
					$otuMembers{$otuNames{$temp1[2]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[2]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[2]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
			case {$_[0] eq 'order'} { 
				$temp1[2]=~ s/ (order)$//;
				if (exists $otuNames{$temp1[3]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[3]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[3]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[3]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[3]}=$ctr;$otuMembers{$otuNames{$temp1[3]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[3]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[3]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}		
			case {$_[0] eq 'family'} { 
				$temp1[2]=~ s/ (family)$//;
				if (exists $otuNames{$temp1[4]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[4]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[4]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[4]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[4]}=$ctr;
					$otuMembers{$otuNames{$temp1[4]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[4]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[4]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
			case {$_[0] eq 'genus'} { 
				$temp1[2]=~ s/ (genus)$//;
				if (exists $otuNames{$temp1[5]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[5]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[5]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[5]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[5]}=$ctr;
					$otuMembers{$otuNames{$temp1[5]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[5]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[5]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
			case {$_[0] eq 'species'} { 
				$temp1[2]=~ s/ (species)$//;
				if (exists $otuNames{$temp1[6]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[6]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[6]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[6]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[6]}=$ctr;
					$otuMembers{$otuNames{$temp1[6]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[6]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[6]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
			case {$_[0] eq 'strain'} { 
				$temp1[2]=~ s/ (strain)$//;
				if (exists $otuNames{$temp1[7]}){#add seqname
					$i=$otuMembers{$otuNames{$temp1[7]}}; $j=$i."\t".$seqName; $otuMembers{$otuNames{$temp1[7]}}=$j;
					#append fas 
					&write_otu_seq($tlevel, $otuNames{$temp1[7]}, 'append', $seqName, $db->seq($seqName)); $nofMembers++;
				}
				else{#init hash items, add taxonomy
					$otuNames{$temp1[7]}=$ctr;
					$otuMembers{$otuNames{$temp1[7]}}=$seqName; 
					$temp[1]=~ s/ \(\w+\)//g;# clean (XXX)
					$otuTaxonomy{$otuNames{$temp1[7]}}=$temp[1];
					#create fas and rep fas
					&write_otu_rep_seq($ctr, $seqName, $db->seq($seqName));
					&write_otu_seq($tlevel, $otuNames{$temp1[7]}, 'open', $seqName, $db->seq($seqName)); $nofMembers++;
				}
			}
		}
		$ctr++;
	}
}
close(INREP);
unlink "${seq}.index";

# write seqs_otus file
@temp= sort{ $a <=> $b} (keys %otuMembers);#0,1,2...
unless(open(OTUS,">${rep}.seqs_otus.txt")){print "not able to open ${rep}.seqs_otus.txt\n\n";exit 1;}
foreach $i (@temp){
	print OTUS "$i\t${otuMembers{$i}}\n"
}
close(OTUS);

# write taxonomy file
@temp= sort{ $a <=> $b} (keys %otuTaxonomy);#0,1,2...
unless(open(OTUT,">${rep}.seqs_rep_tax_assignments.txt")){print "not able to open ${rep}.seqs_rep_tax_assignments.txt\n\n";exit 1;}
foreach $i (@temp){
	print OTUT $i,' ';
	@temp1=split("\t",$otuMembers{$i});
	print OTUT $temp1[0],"\tRoot\;";#first seq is rep
	
	# taxonomy of rep seq is presumed to be taxonomy of OTU
	# can be wrong in cases where taxonomy at levels<tlevel are different for other members 
	@temp1=split(';',$otuTaxonomy{$i});
	foreach $j (0..7){
		switch($j){
			case {$_[0]==0} {print OTUT 'k__',$temp1[$j],';';}
			case {$_[0]==1} {print OTUT 'p__',$temp1[$j],';';}
			case {$_[0]==2} {print OTUT 'c__',$temp1[$j],';';}
			case {$_[0]==3} {print OTUT 'o__',$temp1[$j],';';}
			case {$_[0]==4} {print OTUT 'f__',$temp1[$j],';';}
			case {$_[0]==5} {print OTUT 'g__',$temp1[$j],';';}
			case {$_[0]==6} {print OTUT 's__',$temp1[$j],';';}
			case {$_[0]==7} {print OTUT 's__',$temp1[$j],';';}
		}
	}
	print OTUT "\t1.000\n";
}
close(OTUT);

print STDERR "Taxonomy level: ",$tlevel,"\n";
print STDERR "Number of OTUs: ",scalar(keys %otuNames),"\n";
print STDERR "Number of Seqs in OTUs: ",$nofMembers,"\n";
print STDERR "Avg number of seqs per OTU: ",&SS::round2($nofMembers/(scalar(keys %otuNames))),"\n";
print STDERR "Number of unclassified seqs: ",$nofSeqs-$nofMembers,"\n\n";

&SS::mem_used();
&SS::runtime();
exit;