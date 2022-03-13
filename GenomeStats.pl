#!/opt/perl/bin/perl

die print "usage: <GenomeStats.pl> <infile>" unless scalar(@ARGV) == 1;
$infile = $ARGV[0];

open (FILE, $infile); 


#@genome = <FILE>;

$contigs=0; # total number of contigs
$countNucleotides = 0; #total number of all nucleotides including N ( this equals the length of the genome)
$countA = 0; # Total number of As
$countC = 0; # etc.
$countG = 0;
$countT = 0;
$countN = 0;

%contigLens = ( # Going to build a hash that holds the name of the contigs and their length
#    "contigName" => "contigLength"
);


while ($line=<FILE>) { #read each element (line) in the file one at a time
    if ($line =~ /^>/) {    #then I have found a header!
    #print "$contigLength\n";
    $contigLength = 0; # Set the length of contig to 0, because we JUST started a new contig
    $contigs++; # Add one to the count of contigs
    #What else do I want to do with the header? -- Grab the name of the contig!
    ($name) = $line =~ m/>(.*)\n/; # Put the name of the contig into the variable name
    #print "Contig: $name \t";
    } else { # Not a header -- must be DNA sequence
    #Then I should figure out the length, and count nucleotides!
    chomp $line; 
#    print $line;
    $contigLength = $contigLength + length($line); # As long as we haven't found a new contig yet, keep adding the length of this line to this contig's length
    $contigLens{$name} = $contigLength; # In the hash, set the value of this contig's name to the current length
    $A = $line =~ tr/A/A/; # Count As
    $countA += $A;
    $C = $line =~ tr/C/C/; #etc.
    $countC += $C;
    $G = $line =~ tr/G/G/;
    $countG += $G;
    $T = $line =~ tr/T/T/;
    $countT += $T;
    $N = $line =~ tr/N/N/;
    $countN += $N;
    $countNucleotides = $countNucleotides + length($line);
    }

    
}
#print "$contigLength\n"; # Have to deal with the very last contig, which will not be followed by a contig
$contigLens{$name} = $contigLength; # Set the very last contig name = ultimate contiglength in the hash

print "Number of contigs:\t$contigs\n";
print "Total Nucleotides:\t$countNucleotides\n";
print "Total A:\t$countA\n";
print "Total C:\t$countC\n";
print "Total G:\t$countG\n";
print "Total T:\t$countT\n";
print "Total N:\t$countN\n";

@lengths_array = values %contigLens; # put all of the different lengths of contigs into an array
@sorted_lengths = reverse sort {$a<=>$b} @lengths_array; # So that it is smallest to largest not largest to smallest
print "Min Contig:\t $sorted_lengths[-1]\n";
print "Max Contig:\t $sorted_lengths[0]\n";

$N25 = N_statistic (25, $countNucleotides, \@sorted_lengths); 
$N50 = N_statistic (50, $countNucleotides, \@sorted_lengths); 
$N75 = N_statistic (75, $countNucleotides, \@sorted_lengths); 
$average = $countNucleotides / $contigs;
$median = $sorted_lengths[$contigs/2];
print "Average length:\t$average\n";
print "Median length:\t$median\n";
print "N25 = $N25\n";
print "N50 = $N50\n";
print "N75 = $N75\n";

@lengths_array = ();
foreach $line(@lengths) {
    ($num) = $line =~ m/\t(\d+)\n/;
    push (@lengths_array, $num);

}



############## BEHOLD THE SUBROUTINES #####################

sub N_statistic {
    my ($cutoff, $total, $contig_length_ref) = @_;
    my $i = 0;
    my $sum = 0;

    while ($sum < $total * $cutoff/100) {
        $sum += $$contig_length_ref[$i];
        $i++;
    }
    $N_cutoff = $$contig_length_ref[$i-1];
}


############ Trashcan ###########


# to get N50 - First you have to know how long the genome is -- which is $countNucleotides!
# Then you need to know what 50% of that is -- which should be half of $countNucleotides
# 
# $halfLen = $countNucleotides*.5; # value for the N50
# #print "N50 target is $halfLen\n";
# $quarterLen = $countNucleotides*.25; # value for the N25
# #print "N25 target is $quarterLen\n";
# $threequarterLen = $countNucleotides *.75; # value for the N75
# #print "N75 target is $threequarterLen\n";

# Move through the sorted_length array - adding the contig length to a total until we get to something that is >= $halfLen
# Then we report the length of the contig that we just added as the N50 contig length
# $contigSum = 0;
# 
# foreach $value(@sorted_lengths) { # Go through the sorted values one by one
#     $contigSum = $contigSum + $value; # Add value to the sum of values "so far"
#     if ($contigSum > $threequarterLen) { # 
#         if ($N75) {next;} else {$N75 = $value;}
#     } elsif ($contigSum > $halfLen) {
#         if ($N50) {next;} else {$N50 = $value;}
#         }
#       elsif ($contigSum > $quarterLen) {
#         if ($N25) {next;} else {$N25 = $value;}
#         }
#     else
#     {next;}
#     }


