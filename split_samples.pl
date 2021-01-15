#!/usr/bin/perl 

use strict;
use File::Spec;
use Getopt::Long;
use IO::File;
use Pod::Usage;

# initialize user-definable parameters to default values

my $dir='';
my $beadstudio_file='';
my $rsOnly="true";
my $ncols=3;
my $samples_to_process='all';
my $xcol = 0;
my $ycol = 1;
my $genocol=2;

# set user defined parameters
GetOptions('data_file=s' => \$beadstudio_file,
	   'output_directory=s' => \$dir,
	   'rs_only=s' => \$rsOnly,
           'ncols=s' => \$ncols,
           'samples=s' => \$samples_to_process,
           'xcol=s' => \$xcol,
           'ycol=s' => \$ycol,
           'genocol=s' => \$genocol
           ) or pod2usage(1);


my $extract_dir = $dir."/normalization/";
my $extract_geno_dir =  $dir."/genotypes/";
mkdir $extract_dir;
mkdir $extract_geno_dir;


$beadstudio_file=File::Spec->canonpath($beadstudio_file);
$extract_dir=File::Spec->canonpath($extract_dir);

if($beadstudio_file ne '') {
    unless((-r $beadstudio_file)) {
	print STDERR "Error: Cannot read data file: $beadstudio_file\n";
	exit(0);
    }
}
else {
    print STDERR "Error: No file with data specified!\n";
    pod2usage(1);
}

unless((-d $dir)) {
    print STDERR "Error: Cannot access directory for storing output: " . 
	"$dir\n";
    exit(0);
}

# Main program - Split multiple sample file into one file per sample

# Read header

open(BEADSTUDIO,"$beadstudio_file");
my $beadstudio_error="Error: BeadStudio data in $beadstudio_file " . 
    "is not in the correct format";
my $header=<BEADSTUDIO>;
#chomp($header);
my @header=split(/\t/,$header,-1);
if($header[0] ne 'Name' || $header[1] ne 'Chr' || $header[2] ne 'Position') {
    print STDERR "$beadstudio_error\n";
    print STDERR "The first three column headers should be: " .
	"Name, Chr, Position\n";
    exit(0);
}
#close(BEADSTUDIO);

my @samples;
my $beadstudio_format='normalization';
my $ndup=0;
for(my $i=3;$i<scalar(@header);$i+=$ncols) {    
    if($header[($i+$xcol)]!~/\.X/ && $header[($i+$xcol)]!~/\.B Allele Freq/) {
	print STDERR "$beadstudio_error\n";
	print STDERR "Column " . ($i+$xcol+1) . " header, $header[($i+$xcol)]," .
	    " should be <sample>.X or <sample>.B Allele Freq\n";
	exit(0);
    } 
    if($header[$i]=~/\.B Allele Freq/) {
	$beadstudio_format='segmentation';
    }
    my $sample_name=$header[($i+$xcol)];
    $sample_name=~s/\.X//;
    $sample_name=~s/\.B Allele Freq//;
    my $found = 0;
    foreach my $s (@samples){
	if($sample_name eq $s){
	    $found=1;
	    $ndup++;
        }
    }
    if($found>0){
	print "\t\tWARNING : sample name ".$sample_name." is duplicated\n";
        $sample_name = $sample_name."duplicated_".$ndup;
    }
    push(@samples,$sample_name);
}


my $name="sample_names.txt";
my $filename=File::Spec->canonpath("$extract_dir/$name");
my $fh1=IO::File->new(">$filename");
if(!defined($fh1)) {
    print STDERR "Error: Cannot open file $filename for writing\n";
    exit(0);
}

print $fh1 "Assay\tFilename\tIGV_index\n";

my $igv_index=1;
foreach my $sample (@samples) {
	print $fh1 "$sample" .
	    "\t${sample}_extracted.txt" .
	    "\t$igv_index" .
	    "\n";
	$igv_index++;
}
close($fh1);


my $start=3;
my $stop=scalar(@header);
my $sampleStart = 0;
my $sampleStop = scalar(@samples)-1;
if($samples_to_process ne 'all'){
    my @sp = split(/-/,$samples_to_process);
    $start=$sp[0]*$ncols+3;
    $stop=$sp[1]*$ncols+3;
    $sampleStart=$sp[0];
    $sampleStop=$sp[1];
    
}


print "\t\tProcessing ".($sampleStop-$sampleStart)."/".scalar(@samples)." samples for sample ".$sampleStart." to ".$sampleStop."\n";


my $decimal_error=0;
my $first = 0;  
print "\t\tInitializing outputfiles \n";
    
my @file_handlesXY;
my @file_handlesGeno;
my $i=0;
foreach my $sample (@samples) {
    if($i >= $sampleStart & $i < $sampleStop){
	my $filename=File::Spec->canonpath("$extract_dir/${sample}_extracted.txt");
	my $fh=IO::File->new(">$filename");
	if(!defined($fh)) {
	    print $sample."\n";
	    print STDERR "Error: Cannot open file $filename for writing\n";
	    exit(0);
	}
	my $out_header="Name\tChr\tPosition\tX\tY";
	print $fh "$out_header\n";
	push(@file_handlesXY,$fh);
	$filename = File::Spec->canonpath("$extract_geno_dir/${sample}_extracted.txt");
        $fh=IO::File->new(">$filename");
	if(!defined($fh)) {
	   print $sample."\n";
	   print STDERR "Error: Cannot open file $filename for writing\n";
	   exit(0);
	 }
	 $out_header="Name\tChr\tPosition\tGType";
	 print $fh "$out_header\n";
	 push(@file_handlesGeno,$fh);
     }
     $i++;
}
    
print "\t\tProcessing data \n";

my $k=0;
my $step = 1;
my $nok=0;
while(my $lin = <BEADSTUDIO>) {
     my $snp_ok = 'nok';
     my $snp;     
     my $i = 1;
     if($rsOnly eq "true"){
	 if(substr($lin,0,2) eq "rs"){
	     $snp_ok='ok';  
	 }
     }
     else{
	 $snp_ok='ok';
     }
     if($snp_ok eq 'ok'){ 
	 $nok++;
	 chomp($lin);
	 my @line=split(/\t/,$lin);
	 my $snp = $line[0];
         if(scalar(@header) != scalar(@line)) {
	       print STDERR "$beadstudio_error\n";
	       print STDERR "Not the same number of columns in all rows\n";
	       exit(0);
	 }
         my $fh_index = 0;
	 for(my $i=$start;$i<$stop;$i+=$ncols) { 
	     my $fh=$file_handlesXY[$fh_index];
     	     my $data_line="$line[($i+$xcol)]\t$line[($i+$ycol)]";
	     if($data_line =~ /,/) {
		$decimal_error=1;
		$data_line =~ s/,/./g;
	      }
	      print $fh "$line[0]\t$line[1]\t$line[2]\t$data_line\n";
	      $fh=$file_handlesGeno[$fh_index]; 
	      my $geno=$line[($i+$genocol)];
	      print $fh "$line[0]\t$line[1]\t$line[2]\t$geno\n";
	      $fh_index++;     	   	   		                
	 }
      }
      $k++;
      if($k == $step*10000){
	    print "\t\t\t".$k." lines processed \n";
	    $step++;
      }
}
close(BEADSTUDIO);
print "\t\t=>".$nok." SNPs retained\n";
foreach my $fh (@file_handlesXY) {
    close($fh);
}
foreach my $fh (@file_handlesGeno) {
    close($fh);
}
  
if($decimal_error) {
    print "Warning: commas (,) were found in you BeadStudio data columns. " . 
	"These have been interpreted as decimal points and " . 
	"have been replaced with points (.) in the generated files. " . 
	"BAFsegmentation and tQN requires decimal points to be points (.) and not commas (,).\n"
}

# Documentation

=pod

=head1 NAME

split_samples.pl - split multiple samples in a data file into separate files for use with BAFsegmentation or tQN

=head1 SYNOPSIS

split_samples.pl --data_file=<name> [--output_directory=<name>] [--prefix=<name>]

=head1 OPTIONS

=over 8

=item B<--data_file=<name>>

Specify a file with data from BeadStudio.

=item B<--output_directory=<name>>

Specify a directory to store the generated files in. Default is "extracted".

=item B<--prefix=<name>>

Specify a prefix for the generated file with sample names: <prefix>_sample_names.txt. Default is to use no prefix.


=back

=head1 AUTHORS 

Markus Ringner

Please report bugs to markus.ringner@med.lu.se

=cut

# End of documentation
