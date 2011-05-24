#!/usr/bin/perl -w

$defined = 0;
$output = 0;

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

$lines = $cur - $offset + 1;
print "Found $lines valid lines, from $offset to $cur.\n";
