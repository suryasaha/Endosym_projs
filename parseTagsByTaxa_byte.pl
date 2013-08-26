#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 24, 2012
# Tassel pipeline file

use strict;
use warnings;
use Getopt::Long;

#unless(open(IN,"<:encoding(UTF-8)","merged_chr1-4.tbt.byte.20lines")){print "not able to open \n\n";exit 1;}
#unless(open(IN,"<","merged_chr1-4.tbt.byte.20lines")){print "not able to open \n\n";exit 1;}
unless(open(IN,"<:encoding(UTF-8)","merged_chr1-4.tbt.byte.20lines")){print "not able to open \n\n";exit 1;}

#All values not UNICODE-8!!
while(<IN>){
	print $_,"\n";
}
close(IN);

exit;
