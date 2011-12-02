#!/usr/bin/perl -w
#sqd2sto.pl

use strict;
use Class::Struct;

use vars qw ($opt_f $opt_v);  # required if strict used
use Getopt::Std;
getopts ('fv');

my $easel = "easel/miniapps";
my $esl_reformat = "$easel/esl-reformat ";
my $esl_wussify = "$esl_reformat --wussify ";

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  sqd2sto2.pl [options] <sqdfile> <stofile> \n\n";
        print "options:\n";
        print "-f    :  produce also a fasta file\n";
        print "-v    :  be verbose\n";
        exit;
}
my $sqdfile  = shift;
my $stofile  = shift;

if (! -f $sqdfile)  { die "$sqdfile doesn't exist"; }
if (  -f $stofile)  { system("rm $stofile\n"); }

my $verbose = 0;
if ($opt_v) { $verbose = 1; }

my $dofasta = 0;
if ($opt_f) { $dofasta = 1; }

sqd2sto($sqdfile, $stofile, $dofasta);


if ($verbose) { system "more $stofile\n"; }



sub sqd2sto {
    my ($file, $outfile, $dofasta) = @_;
    
    my $nsq = 0;
    my $name;
    my $src;
    my $includess = 0;
    my $badss     = 0;
    my $sq;
    my $ss;
    my $nline;
    my $line;
    

    open (FILE, "$file") || die; 
    while (<FILE>) {
	if (/^NAM\s+(.+)$/) { 
	    $name = $1;
	    $name =~ s/ /_/g;
 
	    $nsq ++; 
	    $sq = ""; 
	    $ss = ""; 
	    $nline = 0;  
	    if ($verbose) { print "NAM:$name\n"; }
	}
	elsif (/^SRC\s+(.+)$/)     { $src  = $1; }
	elsif (/^DES\s+(.+)$/)     { $name .= "_$1"; }
	elsif (/^SEQ\s+\+SS\s*$/)  { $includess = 1; }
	elsif (/^\+\+$/)        { 
	    
	    print_sto($outfile, $name, $src, $sq, $ss, $dofasta);
	    $includess = 0; 
	    $badss = 0;
	    $sq = ""; 
	    $ss = ""; 
	    $nline = 0; 
	}
	elsif (/^$/)       {
	    if ($includess == 1) { $badss = 1; if ($nline==1) {my $x = 0; while ($x++ < length($line)) { $ss .= "."; } if ($verbose) { print "2SS:$ss\n";   }  } }
	}
	else {
	    if (/^\s*(\D+\s*.*)$/)       { $line = $1; }
	    if (/^\s*\d+\s+(\D+\s*.*)$/) { $line = $1; }
	    
	    $line =~ s/ //g;
	    $line =~ s/\t//g;
	    $line =~ s/\n//g;
	    
	    $nline ++;
	    
	    if    ($includess == 0) { $sq .= $line; if ($verbose) { print "1SQ:$sq\n"; } } 
	    elsif ($badss     == 1) { 
		$sq .= $line;                                           if ($verbose) { print "2SQ:$sq\n"; }
		my $x = 0; while ($x++ < length($line)) { $ss .= "."; } if ($verbose) { print "2SS:$ss\n"; } 
	    }
	    elsif ($badss     == 0) {
		if ($nline%2 == 0) { $ss .= $line; if ($verbose) { print "3SS:$ss|line$nline\n"; } }
		else               { $sq .= $line; if ($verbose) { print "3SQ:$sq|line$nline\n"; } } 
	    }
	}
    }
    close (FILE);
    
    if ($verbose) { print "NSQ = $nsq\n"; }
    
}
	   
sub print_sto {
    my ($stofile, $name, $src, $sq, $ss, $dofasta) = @_;
 
    if ($verbose) { print "SQ:$sq|\n"; print "SS:$ss|\n"; }

    if (length($ss) != length($sq)) { printf "bad lengths %d versus %d\n", length($ss), length($sq); die; }
    if ($sq =~ /\./ || $sq =~ /</ || $sq =~ />/) { printf "bad seq\n", $sq; die; }

    open (STO, ">>$stofile") || die; 

    printf STO "# STOCKHOLM 1.0\n\n";
    printf STO "$name          $sq\n";
    printf STO "\#=GR $name SS $ss\n";
    printf STO "//\n";

    close (STO);

    # use wuss notation
    my $stof = "stofile.wss";
    system("$esl_wussify stockholm $stofile > $stof\n");
    system("mv $stof $stofile\n");

    if ($dofasta) {
	my $fastafile = "$stofile";
	if ($fastafile =~ /^(\S+).sto$/) { $fastafile = "$1.fa"; }
	if (  -f $fastafile)  { system("rm $fastafile\n"); }
	system("$esl_reformat fasta $stofile > $fastafile\n");

    }
 }
