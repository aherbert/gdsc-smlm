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

  -help     Print this help and exit

";

my $help;
GetOptions(
	"help" => \$help,
);

die $usage if $help;

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
            push @file, 
                sprintf('@SuppressWarnings({"unchecked", "unused"})%s', "\n");
        }
        push @file, $_;
	}
	close IN;
	
	next if $looking;

	open (OUT, ">$input") or die "Failed to open '>$input': $!\n";
	print OUT @file;
	close OUT;
}

