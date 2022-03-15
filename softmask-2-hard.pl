#!/usr/bin/perl 
# Usage: perl softmask-2-hard.pl SoftmaskedSequences.fasta HardmaskedSequences.fasta
### ~~~ be wise: this obviously only works one way, so don't throw away SoftmaskedSequences.fasta. ~~~
# pseudo-code: 
# Read file in line by line 
# if line starts with > skip it 
# otherwise replace lowercase acgt with uppercase NNNN
# write line to outfile 
#
my $infile = $ARGV[0];
my $outfile = $ARGV[1]; 

open (FASTA, $infile) or 
    die "I cannot open $infile \n"; 

open (OUT, ">>$outfile") or
    die "I cannot write to $outfile \n";

my @LINES = <FASTA>; 

foreach my $line (@LINES) {
    chomp $line;
    if($line =~ /^>/) { pass; }
    else { 
        $line =~ s/[agct]/N/g;
    }
print OUT $line . "\n";
}

close FASTA;
close OUT;


