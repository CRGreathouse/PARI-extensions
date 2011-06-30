#!/usr/bin/perl

#$/ = '';  		# read paragraph rather than line
undef $/;		# read file rather than line

print "psp=[";
while (<>) {
	$a = $_;
	$a =~ s/\n\d+ /,/g;
	$a =~ s/^1 //;
	$a =~ s/\n+$//;
	print $a;
}
print "];\n";

