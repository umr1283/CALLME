#!/usr/bin/perl 


use strict;
use File::Spec;
use Getopt::Long;
use IO::File;
use Pod::Usage;

my $input_dir='extracted';
my $geno_dir='genotypes';
my $output_file='clusters';
my $samples_lists='';
# set user defined parameters
my $libdir='';
my $logdir='logs';
GetOptions('input_directory=s' => \$input_dir,
	   'geno_directory=s' => \$geno_dir,
	   'output_file=s' => \$output_file,
           'samples_lists=s' => \$samples_lists,
	   'lib_dir=s' => \$libdir,
           'log_dir=s' => \$logdir
	   )or pod2usage(1);


mkdir $logdir;
$input_dir=File::Spec->canonpath($input_dir);

unless((-d $input_dir)) {
    print STDERR "Error: Cannot access directory for reading input: " . 
	"$input_dir\n";
    exit(0);
}


print "Launching R  commands ... \n\n";
my $proc_index = 1;
open(FHI,$samples_lists) || die("cannot open $samples_lists"); 
while(my $sampl_file=<FHI>){
    $sampl_file =~ s/\012//;
    $sampl_file =~ s/\015//;     
    my $R_command="nohup R --vanilla --no-save --slave analysis=cluster_centers input.path=".$input_dir."/ geno.path=".$geno_dir."/ out.file=".$output_file.$proc_index." samples.file=".$sampl_file." < ".$libdir."/cluster_centers.R >& ".$logdir."/cluster_centers_".$proc_index.".out &";
    # Generate parameter file for the R-script performing the normalization
    #open(FILE,">tQN_parameters.txt"); 
    #print FILE "bead.platform<-\"${beadchip}_tQN_clusters.txt\"\n";
    #print FILE "input.path<-\"$input_dir/\"\n"; 
    #print FILE "output.path<-\"$output_dir/\"\n";
    #print FILE "samples.file<-\"$samples_file\"\n";
    #close(FILE);
    # Run R script to generate normalized data
    system("$R_command");
    $proc_index++;
}

