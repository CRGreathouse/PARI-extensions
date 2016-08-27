#!/usr/bin/perl -T
use strict;
use warnings;

use Capture::Tiny 'capture';
use Array::Utils qw(:all);
use Regexp::Assemble;	# only needed if writing regex

# Set these to the locations of these files on your system.
my $src = '../src';
my $whatnow = "$src/gp/whatnow.h";
my $paridesc = "$src/desc/pari.desc";
my $gp2c = '/usr/local/bin/gp2c';


# Get raw output from gp2c.
sub gp2c {
	$ENV{PATH} = "/bin:/usr/bin";	# Otherwise taint mode won't let system run
	my ($stdout, $stderr, $exit) = capture {
		system( $gp2c, ('-l') );
	};
	die if $exit != 0 || $stderr ne '';
	return split("\n", $stdout);
}


# Get a list of functions from gp2c.
sub getFunctions {
	my @func = ();
	my @operator = ();	# currently ignored
	for my $func (gp2c()) {
		if ($func =~ /^([^ ]+) ([^ ]+) (.*)$/i) {
			my $GP = $1;
			my $PARI = $2;	# currently ignored
			my $prototype = $3;	# currently ignored
			if ($GP =~ /^[a-z][a-z0-9_]*$/i) {
				push @func, $GP;
			} else {	
				push @operator, $GP;
			}
		} elsif ($func =~ /^[\s\r\n]*$/) {
		} else {
			print "Unrecognized gp2c line: '$func'\n";
		}
	}
	return sort(@func);
}


sub getRemoved {
	my @removed = ();
	open my $in, '<', $whatnow or die "Could not find src/gp/whatnow.h";
	while(my $line = <$in>) {
		if ($line =~ /^\{"[^" ]+", *_SAME\}/) {
		} elsif ($line =~ /^\{"([a-zA-Z][a-zA-Z0-9_]*)", *_REMOV\}/) {
			push @removed, $1;	# Explicitly marked as removed
		} elsif ($line =~ /\{"([a-zA-Z][a-zA-Z0-9_]*)","([^" ]+)","[^"]*","[^"]*"\}/) {
			my $oldname = $1;
			my $newname = $2;
			if ($oldname eq $newname) {
				# No need to do anything -- only the arguments changed
			} else {
				push @removed, $oldname;
			}
		}
	}
	close $in;
	return sort(@removed);
}


sub getObsolete {
	my @obsolete = ();
	my $name = '';
	open my $in, '<', $paridesc or die "Could not find src/desc/pari.desc";
	while(my $line = <$in>) {
		if ($line =~ /^Function: ([a-zA-Z][a-zA-Z0-9_]*)$/) {	# Function
			$name = $1;
		} elsif ($line =~ /^Function: ([^a-zA-Z]+)$/) {	# Operator
			$name = $1;
		} elsif ($line =~ /^Function: (_[a-zA-Z0-9_]+)$/) {	# gp2c internal function
			$name = $1;
		} elsif ($line =~ /^Function: (_\.[a-zA-Z][a-zA-Z0-9_]*)$/) {	# member function
			$name = $1;
		} elsif ($line =~ /^Function:/) {	# weird, skipped
			# Currently this includes only 'O(_^_)' and '_^s'
			$name = '';
		} elsif ($line =~ /^Obsolete: (.+)$/) {
			my $obsoleteDate = $1;	# currently ignored
			push @obsolete, $name if $name ne '';
		}
	}
	return sort(@obsolete);
}


my @func = getFunctions();
my @removed = getRemoved();
my @obsolete = getObsolete();

# If it's in both, one function was removed and another with the same name added.
@removed = array_minus(@removed, @func);

my @bad = sort(unique(@removed, @obsolete));
my @good = sort(array_minus(@func, @obsolete));

print 'Good functions: ' . join(', ', @good) . "\n";
print 'Bad functions: ' . join(', ', @bad) . "\n";
print 'In total there are ' . scalar(@good) . ' good and ' . scalar(@bad) . " bad functions.\n";

#my $ra = Regexp::Assemble->new;
#$ra->add(@func);
#my $r = $ra->re;
#$r =~ s/\(\?:/\(/g;
#print "$r\n";

