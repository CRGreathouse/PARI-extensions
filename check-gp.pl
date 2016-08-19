#!/usr/bin/perl -T
use strict;
use warnings;

use enum qw(Normal BlockComment LineComment String);
use File::Temp 'tempfile';
use Capture::Tiny 'capture';

# Checks PARI programs for balanced prarentheses and other issues.
# Usage: ./pari.pl | tee fixme.html

#my @opname = qw(_||_ _&&_ _===_ _==_ _!=_ _>=_ _>_ _<=_ _<_ _-_ _+_ _<<_ _>>_ _%_ _\\/_ _\\_ _/_ _*_ _^_ __ _-- _++ _-=_ _+=_ _<<=_ _>>=_ _%=_ _\\/=_ _\\=_ _/=_ _*=_ +_ -_ !_ _! _' _~ [_.._] [_|_<-_,_] [_|_<-_,_;_] % %# #_ _[_,_]);
my $date = qr/(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) (0[1-9]|[12][0-9]|3[01]) 20[01][0-9]/;
my $name = qr/(?:[A-Za-z.' -]+|_[A-Za-z.' -]+_)/;
my $part = qr/[a-z0-9]+(?:-[a-z0-9]+)/;
my $email = qr/\([a-z0-9_.+-]\(AT\)$part(?:\.$part)+\)/i;

sub codeErr {
	my $code = shift;
	my ($fh, $filename) = tempfile();
	print $fh "$code\n";
	close $fh;
	$ENV{PATH} = "/bin:/usr/bin";	# Otherwise taint mode won't let system run
	my ($stdout, $stderr, $exit) = capture {
		system( '/usr/local/bin/gp2c', ('-d', '-G', $filename) );
	};
	unlink $filename;
	return undef if $exit == 0;
	return fmtErr($stderr, $filename);
}

sub fmtErr {
	my $err = shift;
	my $file = shift;
	chomp $err;
	$err =~ s/Error:\Q$file:/Line /g;
	$err =~ s/\s*\Q$file:/Line /g;
	$err =~ s/\n([0-9]+ e|E)rrors? found: aborting\.\.\.//g;
	$err = html($err);
	$err =~ s/\n/<br>/g;
	return $err;
}

sub skipCode {
	my $code = shift;
#return 1 unless $code =~ /Michel Marcus/;	# Testing
	return 1 if $code eq '';
	return 1 if $code =~ /^See [Ll]inks?\.?$/;
	return 1 if $code =~ /^See .+ links?\.?$/;
	return 1 if $code =~ /^See [Ll]inks? section\.?$/;
	return 1 if $code =~ /^See A[0-9]{6}\.?$/;
#	return 1 if $code =~ /\[From $name, $date\]$/;
#	return 1 if $code =~ /-\s*$name, $date$/;
#	return 1 if $code =~ /-\s*$name $email, $date$/;
	return 1 if $code =~ /Michael Somos/;
	return 0;
}

# Returns (bad, total, skipped)
sub checkCode {
	my $code = shift;
	return (0, 0, 1) if skipCode($code);
	my $Anum = shift;
	my $errParse = parseCode($code);
	my $errgp2c = codeErr($code);

	return (0, 1, 0) unless defined $errParse or defined $errgp2c;

	$code = formatCode($code);
	my $err;
	$errParse = '<span title="error from Perl parser" class="perl">' . html($errParse) . '</span>' if defined $errParse;
	$errgp2c = '<span title="error from gp2c" class="gp2c">' . $errgp2c . '</span>' if defined $errgp2c;
	if (defined($errParse) and defined($errgp2c)) {
		$err = $errParse . "</span><br>\n" . $errgp2c;
	} elsif (defined $errParse) {
		$err = $errParse;
	} else {
		$err = $errgp2c;
	}
	print "<tr><td><a href='https://oeis.org/edit?seq=$Anum'>$Anum</a></td>";
	print "<td style='width: 20%'>$err</td>";
	print "<td><code class='language-pari-gp'>$code</code></td></tr>\n";
	return (1, 1, 0);
}

# Parse the given code, removing comments and replacing strings with "".
# If errors are found, return an error string. Otherwise, pass cleaned
# code to checkClean where additional tests are performed.
sub parseCode {
	my $code = shift;
	my @stack = ();
	my $status = Normal;
	my $level = 0;	# Brace level, cannot be > 1
	my $last = '';	# previous character;
	my $clean = '';
	for my $c (split //, $code) {
		if ($status == Normal) {
			if ($c eq '"') {
				$status = String;
			} elsif ($c eq '(') {
				push @stack, ')';
			} elsif ($c eq '[') {
				push @stack, ']';
			} elsif ($c eq '{') {
				return "Nested braces not allowed" if $level > 0;
				$level = 1;
			} elsif ($c =~ /\)|\]/) {
				if (scalar(@stack) == 0) {
					my $parens = $clean;
					$parens =~ tr/()[]//cd;
					return "Stack underflow, unexpected $c, had $parens$c...";
				}
				my $d = pop @stack;
				return "Mismatched delimiters, expected $d, found $c" unless $d eq $c;
				$level = 0 if $c eq '}';
			} elsif ($c =~ /\}/) {
				return "Closing brace without matching open brace" if $level == 0;
				$level = 0;
			} elsif ($c eq '*' && $last eq '/') {
				$status = BlockComment;
			} elsif ($c eq '\\' && $last eq '\\') {
				$status = LineComment;
			} elsif ($c eq "\n") {
				# If $level > 0, you are inside braces and newlines are permitted.
				if (scalar(@stack) > 0 && $level == 0) {
					my $parens = $clean;
					$parens =~ tr/()[]//cd;
					return "Newline with unmatched delimiters, had $parens...";
				}
			}

			# Check new status, in case it changed above.
			if ($status == BlockComment || $status == LineComment) {
				chop $clean;	# Drop comments
			} elsif ($status == String) {
				$clean .= '""';	# Hack: replace all strings with empty strings
			} else {
				$clean .= $c;
			}
		} elsif ($status == BlockComment) {
			if ($c eq '/' && $last eq '*') {
				$status = Normal;
			}
		} elsif ($status == LineComment) {
			if ($c eq "\n") {
				$status = Normal;
			}
		} elsif ($status == String) {
			if ($last eq '\\') {
			} elsif ($c eq '"') {
				$status = Normal;
			} elsif ($c eq "\n") {	# Are there any other characters not allowed in strings?
				return "Run-away string";
			}
		}

		$last = $c;
	}
	if (scalar(@stack) > 0) {
		my $parens = $clean;
		$parens =~ tr/()[]{}//cd;
		return "Missing closing $stack[0], had $parens" if scalar(@stack) == 1;
		return "Unclosed delimiters, had $parens, expected " . join('', reverse @stack);
	}
	return 'Block comment still open' if $status == BlockComment;
	return 'Unterminated string' if $status == String;
	return checkClean($clean) if $status == Normal || $status == LineComment;
	return 'Bug in pari.pl, unknown status';
}


sub formatCode {
	my $code = shift;
	my @stack = ();
	my $status = Normal;
	my $last = '';
	my $formatted = '';
	for my $c (split //, $code) {
		if ($status == Normal) {
			if ($c eq '"') {
				$status = String;
			} elsif ($c eq '(') {
				push @stack, ')';
			} elsif ($c eq '[') {
				push @stack, ']';
			} elsif ($c eq '{') {
				push @stack, '}';
			} elsif ($c =~ /[\)\]\}]/) {
				# Ignore underflow and mismatched delimiters.
				pop @stack;
			} elsif ($c eq '*' && $last eq '/') {
				$status = BlockComment;
			} elsif ($c eq '\\' && $last eq '\\') {
				$status = LineComment;
			}
			
			# Check new status, in case it changed above.
			if ($status == BlockComment) {
				chop $formatted;
				$formatted .= "<span class='comment'>/*";
			} elsif ($status == LineComment) {
				chop $formatted;
				$formatted .= "<span class='comment'>\\\\";
			} elsif ($status == String) {
				$formatted .= '<span class="string">"';
			} else {
				$formatted .= html($c);
			}
		} elsif ($status == BlockComment) {
			if ($c eq '/' && $last eq '*') {
				$status = Normal;
				$formatted .= "/</span>";
			} else {
				$formatted .= html($c);
			}
		} elsif ($status == LineComment) {
			if ($c eq "\n") {
				$status = Normal;
				$formatted .= "</span><br>\n";
			} else {
				$formatted .= html($c);
			}
		} elsif ($status == String) {
			if ($last eq '\\') {
				chop $formatted;	 # FIXME: This doesn't work properly if the last thing in $formatted is a <span>...</span>.
				$formatted .= "<span class='escape'>\\" . html($c) . "</span>";
			} elsif ($c eq '"') {
				$status = Normal;
				$formatted .= '"</span>';
			} elsif ($c eq "\n") {	# Are there any other characters not allowed in strings?
				# Run-away string, terminate
				$status = Normal;
				$formatted .= "</span>";
			} else {
				$formatted .= html($c);
			}
		}
		
		$last = $c;
	}

	# Handle unclosed delimiters
	#$formatted .= '</span>' x (scalar(@stack));
	$formatted .= '</span>' if $status != Normal;

	# Flag deprecated/removed code
	# This is a hack, it dosn't use proper parsing, but it's good enough here.
	$formatted =~ s/((^|[^|\s])\s*)\|(\s*($|[^|\s]))/$1<span class="error">|<\/span>$3/g;
	$formatted =~ s/((?<!&amp;)\s*)&amp;(\s*(?!&amp;))/$1<span class="error">&amp;<\/span>$2/g;

	return $formatted;
}


# Given delimiter-balanced code without comments or nonempty strings, check for errors.
# If any are found, return an error string. Otherwise, return undef.
# At the moment the only checks are for & and |.
sub checkClean {
return undef;
	my $code = shift;
	
	my $dep = '';
	$dep .= '|' if $code =~ /(^|[^|\s])\s*\|\s*($|[^|\s])/;
	$dep .= '&' if $code =~ /(^|[^&\s,])\s*&\s*($|[^&\s])/;	# Note: a single & is a pointer (address-of) operator. Right now these can appear only in 
	
	if ($dep ne '') {
		# http://pari.math.u-bordeaux.fr/archives/pari-announce-11/msg00001.html
		# http://pari.math.u-bordeaux.fr/cgi-bin/gitweb.cgi?p=pari.git;a=commitdiff;h=c1f32470f80cd544e67c7603344e3b40920f686c
		return "Use of operator |, deprecated since 2011 and removed in 2012; use || instead" if $dep eq '|';
		return "Use of operator &, deprecated since 2011; use && instead" if $dep eq '&';
		return "Use of operators | and &, deprecated since 2011; use || and && instead";
	}
	
	return undef;
}

sub html {
	my $txt = shift;
	$txt =~ s/&/&amp;/g;
	$txt =~ s/</&lt;/g;
	$txt =~ s/>/&gt;/g;
	$txt =~ s/\n/<br>\n/g;
	return $txt;
}

sub commify {
	my $text = reverse shift;
	$text =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}

print "<!doctype html>
<html><head>
<meta charset='utf-8'>
<title>GP bugs in the OEIS</title>
<style>
	table {border-collapse:collapse}
	table td {border: 1px dotted #444}
	code {color: #111}
	.comment {color: #090}
	.string {color: #006}
	.escape {color: #622}
	.error {color: #f00}
	.perl {color: #003}
	.gp2c {color: #610}
</style>
</head><body>
<h1>GP bugs in the OEIS</h1>
<table>
<tr><th>Sequence</th><th style='width: 20%'>Error code</th><th>Code</th></tr>\n";

open FILE, '<', 'cat25';
my $bad = 0;	# Instances of bad code discovered.
my $totalPARI = 0;
my $skipped = 0;
my $seq = '';
my $curCode = '';
my $inPARI = 0;
while (<FILE>) {
	if(/^%([A-Za-z]) (A......) (.*)/) {
		my ($type, $Anum, $txt) = ($1, $2, $3);
		if ($type eq 'o') {
			if ($txt =~ /^\(([^()]+)\)\s*(.*)$/) {
				checkCode($curCode, $seq);
				$curCode = '';
				my $lang = $1;
				my $code = $2;
				if ($lang eq 'PARI') {
					$curCode = $code;
					$seq = $Anum;
					$inPARI = 1;
				} else {
					$inPARI = 0;
				}
			} elsif ($inPARI) {
				$curCode .= "\n$txt";
			}
		} elsif ($inPARI) {
			# ($bad, $totalPARI, $skipped) += checkCode($curCode, $seq);
			my ($a, $b, $c)= checkCode($curCode, $seq);
			$bad += $a;
			$totalPARI += $b;
			$skipped += $c;

			$curCode = '';
			$inPARI = 0;
		}
	}
}
my $percent = sprintf("%.1f%%", 100.0 * $bad / $totalPARI);
$bad = commify($bad);
$totalPARI = commify($totalPARI);
$skipped = commify($skipped);
print "</table>
<p>Found $bad PARI programs with apparent errors out of $totalPARI ($percent). $skipped were skipped.</p>
</body></html>\n";

