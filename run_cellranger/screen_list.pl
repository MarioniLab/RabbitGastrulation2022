#!/usr/bin/perl -w

use warnings;

$n_args = @ARGV;
if (($n_args != 2) && ($n_args != 3)) {
    print "Filter a fasta file for/against a list of entry names.\n";
    print "Usage: ./screen_list.pl <list> <fasta file> <<keep?>>\n";
    print "If a third argument is given the list entries will be kept,\n";
    print "otherwise they will be rejected.\n";
    exit;
}

my $invert_sense = 0;
if ($n_args == 3) {$invert_sense = 1;}

#
# 7/08/2005
# Modified the variable i to remove trailing white spaces
# -Jasmyn
#
open(F,$ARGV[0]) || die "Couldn't open file $ARGV[0]\n";
my %bad_ones = (); 
while (my $i = <F>) {
    chomp $i;
#    print "[$i]\n";    #For testing purposes.
    $i =~ s/\s*$//g;
    $bad_ones{$i} = 1;
}
close F;

#
# 6/29/2004
# Modified to allow for files that are larger than 2 GB.
# -Harris
#
#open(F,$ARGV[1]) || die "Couldn't open file $ARGV[1]\n";
open(F,"cat $ARGV[1] |") || die "Couldn't open file $ARGV[1]\n";
my $good_one = 1;
#my $line_counter=0;
#my $read_counter=0;
while (my $i = <F>) {
#    $line_counter++;
    if ($i =~ />/) {
#	$read_counter++;
	my ($id) = $i =~ />(\S+)/;
	if (exists($bad_ones{$id})) {$good_one = 0;}
	else {$good_one = 1;}
    }

    if (($invert_sense == 0) &&  ($good_one == 1)) {print $i;}
    elsif (($invert_sense == 1) &&  ($good_one == 0)) {print $i;}
}
close F;

#print "Line Count:\t$line_counter\n";
#print "Read Count:\t$read_counter\n";
