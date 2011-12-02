#!/usr/bin/perl -w
#sqd2sto.pl

use strict;
use Class::Struct;

use vars qw ($opt_v);  # required if strict used
use Getopt::Std;
getopts ('v');

my $easel = "easel/miniapps";

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  sqd2sto2.pl [options] <sqdfile> <stofile> \n\n";
        print "options:\n";
        print "-v    :  be verbose\n";
        exit;
}
my $sqdfile  = shift;
my $stofile  = shift;

if (! -f $sqdfile)  { die "$sqdfile doesn't exist"; }

my $verbose = 0;
if ($opt_v) { $verbose = 1; }

sqd2sto($sqdfile, $stofile);


if ($verbose) { system "more $stofile\n"; }



sub sqd2sto {
    my ($file, $outfile) = @_;
    
    my $nsq = 0;
    my $name;
    my $src;
    my $includess = 0;
    my $sq;
    my $ss;
    my $nline;


    open (FILE, "$file") || die; 
    while (<FILE>) {
	if    (/^NAM\s+(\S+)$/) { $name = $1; $nsq ++; $sq = ""; $ss = ""; $nline = 0; }
	elsif (/^SRC\s+(.+)$/)  { $src .= $1; }
	elsif (/^SEQ\s+\+SS$/)  { $includess = 1; }
	elsif (/^\+\+$/)        { 

	    print_sto($outfile, $name, $src, $sq, $ss);
	    $includess = 0; 
	    $sq = ""; 
	    $ss = ""; 
	    $nline = 0; 
	}
	elsif (/^(\S+)\s*$/)       { 
	    $nline ++;
	    
	    if ($includess == 0) { $sq .= $1; } 
	    else {
		if ($nline%2 == 0) { $ss .= $1; } 
		else               { $sq .= $1; } 
	    }
	}
	
    }
    close (FILE);

    print "NSQ = $nsq\n";
    
}
	   
sub print_sto {
    my ($stofile, $name, $src, $sq, $ss) = @_;

    open (STO, ">$stofile") || die; 

    printf STO "# STOCHOLM 1.0\n";
    printf STO "$name      \t$sq\n";
    printf STO "\#=GR $name\t$ss\n";
    printf STO "//\n";

    close (STO);
}
