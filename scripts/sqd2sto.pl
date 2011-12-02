#!/usr/bin/perl -w
#sqd2sto.pl

use strict;
use Class::Struct;

use vars qw ($opt_v);  # required if strict used
use Getopt::Std;
getopts ('v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  sqd2sto2.pl [options] <sqdfile>\n\n";
        print "options:\n";
        print "-v    :  be verbose\n";
        exit;
}
my $sqdfile  = shift;
if (! -f $sqdfile)  { die "$sqdfile doesn't exist"; }

