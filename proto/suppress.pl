#!/usr/bin/perl -w

use Getopt::Long;
use File::Basename;

my $prog = basename($0);
my $tmpFile = "/tmp/$prog.$$";
my $out = "$prog.out";
my $debug = 0;

my $usage = "
  Program to add a SuppressWarnings header to .java proto output files

Usage:

  $prog input [...]

Options:

  input     Java proto file(s)
  --help    Print this help and exit
  --unchecked
  --unused
  --deprecation

";

my $help;
GetOptions(
	"help" => \$help,
    "unchecked" => \$unchecked,
    "unused" => \$unused,
    "deprecation" => \$deprecation,
);

die $usage if $help;

my @warning;
push @warning, "unchecked" if ($unchecked);
push @warning, "unused" if ($unused);
push @warning, "deprecation" if ($deprecation);

exit(0) unless @warning;

$warning = '@SuppressWarnings({';
for ($i=0; $i<=$#warning; $i++) {
    $warning .= ', ' if $i;
    $warning .= '"' . $warning[$i] . '"';
}
$warning .= "})\n";

@ARGV or die $usage;

my @files;
for (@ARGV) {
	if (m/\*/) {
		push @files, glob "$_";
	} else {
		push @files, $_;
	}
}

for $input (@files)
{
	open (IN, $input) or die "Failed to open '$input': $!\n";
	my @file;
	$looking = 1;
	while (<IN>)
	{
        if ($looking && m/ class /)
        {
            $looking = 0;
            push @file, $warning;
        }
        push @file, $_;
	}
	close IN;
	
	next if $looking;

	open (OUT, ">$input") or die "Failed to open '>$input': $!\n";
	print OUT @file;
	close OUT;
}

