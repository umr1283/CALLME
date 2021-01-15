#!/usr/bin/perl 

# Copyright (C) 2008 Markus Ringner
# 
# This file is part of tQN,
# http://baseplugins.thep.lu.se/wiki/se.lu.onk.IlluminaSNPNormalization
# 
# tQN is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
# 
# tQN is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Class Discoverer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307 USA

use strict;
#use warnings;

use File::Spec;
use Getopt::Long;
use IO::File;
use Pod::Usage;

######################################################

# NOTE!!! Change here depending on how you run R on your system

###
# Mac OS X and Linux
my $tQNdir = '/ep10/disks/SANA8/boris/home/bin/tQN-1.1.2';
#my $R_command="R --vanilla --no-save --slave < ".$tQNdir."/tQN.R";
###
# Windows
# Note that we are using Rscript, which is a part of the R distribution.
#my $R_windows=File::Spec->canonpath('C:/"Program Files"/R/R-2.7.0/bin/Rscript');
#my $R_command="$R_windows --vanilla tQN.R";

#######################################################

# initialize user-definable parameters to default values
my $beadchip='';
my $input_dir='extracted';
my $output_dir='normalized';
my $output_format='PennCNV';
my $samples_lists=$input_dir.'/list.txt';
my $all_samples='extracted/sample_names.txt';
my $logdir ='';
my $script_generate="/ep10/disks/SANA8/boris/home/scripts/CALLME/tQN_generate_normalized_data.pl";
#my $suffix='';
# set user defined parameters

GetOptions('beadchip=s' => \$beadchip,
	   'input_directory=s' => \$input_dir,
	   'output_directory=s' => \$output_dir,
	   'output_format=s' => \$output_format,   
           'samples_lists=s' => \$samples_lists,
           'all_samples=s' => \$all_samples,
	   'tQN_dir=s' => \$tQNdir
            )or pod2usage(1);

$input_dir=File::Spec->canonpath($input_dir);
mkdir $output_dir;
$logdir = $output_dir."/logs";
mkdir $logdir;
$output_dir=File::Spec->canonpath($output_dir);

unless((-d $input_dir)) {
    print STDERR "Error: Cannot access directory for reading input: " . 
	"$input_dir\n";
    exit(0);
}

unless((-d $output_dir)) {
    print STDERR "Error: Cannot access directory for storing output: " . 
	"$output_dir\n";
    exit(0);
}

# There should be a cluster file for the beadchip
my $clusterfile=File::Spec->canonpath($tQNdir."/".${beadchip}."_tQN_clusters.txt");
if($beadchip ne '') {
    unless((-r $clusterfile)) {
	print STDERR "Error: There is no clusterfile for the beadchip: $beadchip\n";
	exit(0);
    }
}
else {
    print STDERR "Error: No BeadChip type is specified\n";
    pod2usage(1);
}

unless($output_format eq 'PennCNV' || 
       $output_format eq 'QuantiSNP' || 
       $output_format eq 'BeadStudio' || 
       $output_format eq 'BAFsegmentation') {
    print STDERR "Error: Unsupported output format: $output_format\n";
    exit(0);
}

print "\tLauching R  commands ... \n\n";

my @samples_list_file = ();
my $proc_index = 1;
open(FHI,$samples_lists) || die("cannot open $samples_lists"); 
while(my $sampl_file=<FHI>){
    $sampl_file =~ s/\012//;
    $sampl_file =~ s/\015//;
    push @samples_list_file, $sampl_file;
    my $R_command=" nohup R --vanilla --no-save --slave analysis=tQN_norm bead.platform=".$beadchip."_tQN_clusters.txt input.path=".$input_dir."/ output.path=".$output_dir."/ samples.file=".$sampl_file." cluster.path=".$tQNdir."/ < ".$tQNdir."/tQN.R >& ".$logdir."/normalization_".$proc_index.".out &";
    # Run R script to generate normalized data
    system("$R_command");
    $proc_index++;
}
print "\tWaiting for R processes to finish ... \n\n";
my $nprocs = 1;
while($nprocs > 0){
system("ps ux > tmpfile");
open(FHps,"tmpfile") || die("cannot open file tmpfile");
   my $lin = <FHps>;
   $nprocs = 0;
   while($lin=<FHps>){
     chomp($lin); 
     my @elts2 = split(/\s+/,$lin);
     if(scalar(@elts2)>13){
       #my $analysis = '';
       #$analysis=$elts2[15];
       if($elts2[14] eq 'analysis=tQN_norm'){
            $nprocs ++; 
       }
     }
   }
   close(FHps);
}
system("rm -f tmpfile");

# Generate normalized data in the specified output format
print "\tGenerating normalized data in the specified output format...\n";
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
my %IGV_index;
my $toto=1;
while(my $line=<SAMPLES>) {
    chomp($line);
    my @items=split(/\t/,$line);
    if(-r $output_dir."/".$items[0]."_tQN.txt"){
	push(@samples,$items[0]);
	$IGV_index{$items[0]}=$items[2];
    }
    else{
	print "WARNING : quantile normalization failed for sample ".$items[0].", because more than 20% of its data is missing\n";
    }
}
close(SAMPLES);

# A single file per sample
# In addition a file samples_names.txt is generated which is used by BAFsegmentation


my $samples_norm_file=File::Spec->canonpath("$output_dir/sample_names.txt");
open(SAMPLES_NORM,">$samples_norm_file");
print SAMPLES_NORM "$ok_sample_header\n";
foreach my $sample (@samples) {
 print SAMPLES_NORM "${sample}\t${sample}_tQN_${output_format}.txt\t" .
	"$IGV_index{$sample}\n";
}
close(SAMPLES_NORM);

my $toto=1;
foreach my $samp_list (@samples_list_file){
    system("nohup perl ".$script_generate." --output_directory=".$output_dir." --output_format=".$output_format." --samples_list=".$samp_list." >& ".$output_dir."/logs/generate".$toto.".out &");
    $toto++;
}

print "\tWaiting for all processes to finish ...\n";
waitFor($script_generate);

sub waitFor(){

    my $keyword = $_[0];
    my $nproc = 1;
    while($nproc > 0){
	system("ps ux > tmpfile");
	open(FHps,"tmpfile") || die("cannot open file tmpfile");
	my $lin = <FHps>;
	$nproc = 0;
	while($lin=<FHps>){
	    chomp($lin);
	    if($lin =~ m/$keyword/g){
	       $nproc ++; 
	    }
	}
	close(FHps);
	unlink "tmpfile";
    }
}









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
