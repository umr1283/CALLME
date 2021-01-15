#!/usr/bin/perl -

use strict;

use File::Spec;
use Getopt::Long;
use IO::File;
use Pod::Usage;


# initialize user-definable parameters to default values
my $output_dir='';
my $output_format='PennCNV';
my $all_samples='';
#my $suffix='';
# set user defined parameters

GetOptions('output_directory=s' => \$output_dir,
	   'output_format=s' => \$output_format,   
           'samples_list=s' => \$all_samples
            )or pod2usage(1);

$output_dir=File::Spec->canonpath($output_dir);

unless((-d $output_dir)) {
    print STDERR "Error: Cannot access directory for storing output: " . 
	"$output_dir\n";
    exit(0);
}


# Generate normalized data in the specified output format
print "\nGenerate normalized data in the specified output format\n";
my $sample_file=File::Spec->canonpath($all_samples);
open(SAMPLES,"$sample_file") || 
    die "Error: Cannot open sample list file: $sample_file\n";
my $sample_header=<SAMPLES>;
chomp($sample_header);
my $ok_sample_header="Assay\tFilename\tIGV_index";
if($sample_header ne $ok_sample_header) {
    print STDERR "Error: Unexpected file format in sample list file $sample_file: " .
	"The header is not $ok_sample_header\n";
    exit(0);
}
my @samples;
my $nsamp=0;
while(my $line=<SAMPLES>) {
    chomp($line);
    my @items=split(/\t/,$line);
    push(@samples,$items[0]);
    $nsamp++;
}
close(SAMPLES);

# A single file per sample
# In addition a file samples_names.txt is generated which is used by BAFsegmentation

if($output_format eq 'PennCNV' || $output_format eq 'QuantiSNP' || 
   $output_format eq 'BAFsegmentation') {
    my $progress = 1;
    foreach my $sample (@samples) {
        print "  Processing sample ".$sample." (".$progress."/".$nsamp.")\n"; 
   	my $filename=File::Spec->canonpath("$output_dir/${sample}_tQN.txt");
	my $fh=IO::File->new("$filename");
	if(!defined($fh)) {
	    print "WARNING: Cannot open sample file: $filename\n";
	}
	else{
	    $progress++;
	    my $outfile=File::Spec->canonpath("$output_dir/${sample}_tQN_${output_format}.txt");
	    open(OUTPUT,">$outfile");
	    # PennCNV
	    my $headerline="Name\tChr\tPosition\t${sample}.B Allele Frequency\t" .
		"${sample}.Log R Ratio";
	    print OUTPUT "$headerline\n"; 
	    my $header=<$fh>;
	    chomp($header);
	    my $ok_header="Name\tChr\tPosition\ttQN B Allele Frequency\ttQN Log R Ratio\ttQN X\ttQN Y";
	    if($header ne $ok_header) {
		print STDERR "Error: Incorrect column headers in: $filename\n" .
		    "expected: $ok_header\n";	    
		exit(0);
	    }
	    while(my $line=<$fh>) {
		chomp($line);
		my @items=split("\t",$line,-1);
		if($items[3] ne "NA" && $items[4] ne "NA") {
		    # Penn CNV and BAFsegmentation format
		    my $outline="$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]";
		    if($output_format eq  'QuantiSNP') {
			$outline="$items[0]\t$items[1]\t$items[2]\t$items[4]\t$items[3]";
		    }
		    print OUTPUT "$outline\n";
		}
	    }
	    close(OUTPUT);
	    close($fh);
	}
    }
}

print "Analysis finished\n";



# End of main program

# Documentation

=pod

=head1 NAME

tQN_normalize_samples.pl - normalize samples using tQN

=head1 SYNOPSIS

tQN_normalize_samples.pl --beadchip=<name> [--input_directory=<name>] [--output_directory=<name>] [--output_format=<name>] 

=head1 OPTIONS

=over 8

=item B<--beadchip=<name>>

Specify which BeadChip the whole genome genotyping data was generated
on. Supported BeadChips include "humanhap550" and "humancnv370-duo". For a
complete list of currently supported BeadChips check for which
BeadChips there are cluster files in the "lib" directory.

=item B<--input_directory=<name>>

Specify a directory with files for individual samples with BeadStudio
data. The directory should contain a file "sample_names.txt" with a
list of samples to analyse. Files in this directory can be generated
using "split_beadstudio_samples.pl". Default is "extracted".

=item B<--output_directory=<name>>

Specify a directory to store the generated files in. Default is "normalized".

=item B<--output_format=<name>>

Specify the format of the generated normalized data. Alternatives are "BeadStudio", "PennCNV", "QuantiSNP", and "BAFsegmentation". Default is "BeadStudio".

=back

=head1 AUTHORS 

Markus Ringner

Please report bugs to markus.ringner@med.lu.se

=cut

# End of documentation
