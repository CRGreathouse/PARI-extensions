#!/usr/bin/perl -w

use strict;
use warnings;
my $defined = 0;
my $output = 0;
my $cur;
my $offset;

while (<>) {	# For each line of the b-file
	$_ =~ s/\s*#.*//;	# Remove comments with optional whitespace.
	if (m/^(-?[0-9]+) (-?[0-9]+)$/) {
		if ($defined == 0) {
			$defined = 1;
			$cur = $1;
			$offset = $1;
			if ($output) {
				print "v=[$2";
			}
		} else {
			$cur++;
			if ($cur != $1) {
				print "Expected line #$cur; got line #$1\n";
				exit;
			}
			if ($output) {
				print ",$2";
			}
		}
	} else {
		if (m/^(-?[0-9]+) -?[0-9]*\r$/) {
			print "File contains Windows CR/LF pairs.\n";
			exit;
		}
		if (m/./) {
			print "Malformed line: $_\n";
			exit;
		}
	}
}

if ($output) {
	print "];\n\n";
}

my $lines = $cur - $offset + 1;
print "Found $lines valid lines, from $offset to $cur.\n";
