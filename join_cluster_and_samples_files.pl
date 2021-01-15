#!/usr/bin/perl 

use strict;
use Pod::Usage;
use Getopt::Long;

# initialize user-definable parameters to default values

my $file_list_clusters='';
my $file_list_samples='';
my $out_file_clusters = '';
my $out_file_samples = '';


# set user defined parameters
GetOptions('file_list_clusters=s' => \$file_list_clusters,
	   'file_list_samples=s' => \$file_list_samples,
	   'out_file_clusters=s' => \$out_file_clusters,
	   'out_file_samples=s' => \$out_file_samples
	   ) or pod2usage(1);

#print "######### MERGING CLUSTER FILES ##########\n\n";
my @files;
open(FILES,$file_list_clusters) || die("cannot open file ".$file_list_clusters);
while(my $line=<FILES>){
    $line =~ s/\n//g;
    push @files, $line;
}
close(FILES);
my $k=0;
my $header;
my %AAT=();
my %ABT=();
my %BBT=();
my %AAR=();
my %ABR=();
my %BBR=();
my %AAn=();
my %ABn=();
my %BBn=();
foreach my $f(@files){
    open(FHI,$f) || die("cannot open file ".$f);
    $header = <FHI>;
    while(my $line=<FHI>){
	chomp($line);
	my @elts = split(/\t/,$line);
	if($k==0){
	    $AAT{$elts[0]} = 0;
	    $ABT{$elts[0]} = 0;
	    $BBT{$elts[0]} = 0;
	    $AAR{$elts[0]} = 0;
	    $ABR{$elts[0]} = 0;
	    $BBR{$elts[0]} = 0;
	    $AAn{$elts[0]} = 0;
	    $ABn{$elts[0]} = 0;
	    $BBn{$elts[0]} = 0;
	}
	if($elts[1] ne 'NA'){
	    $AAT{$elts[0]} = $AAT{$elts[0]} + $elts[1];
	    $AAn{$elts[0]} = $AAn{$elts[0]} + 1;
	}
	if($elts[2] ne 'NA'){
	    $ABT{$elts[0]} = $ABT{$elts[0]} + $elts[2];
	    $ABn{$elts[0]} = $ABn{$elts[0]} + 1;
	}
	if($elts[3] ne 'NA'){
	    $BBT{$elts[0]} = $BBT{$elts[0]} + $elts[3];
	    $BBn{$elts[0]} = $BBn{$elts[0]} + 1;
	}
	if($elts[4] ne 'NA'){
	    $AAR{$elts[0]} = $AAR{$elts[0]} + $elts[4];
	}
	if($elts[5] ne 'NA'){
	    $ABR{$elts[0]} = $ABR{$elts[0]} + $elts[5];
	}
	if($elts[6] ne 'NA'){	
	    $BBR{$elts[0]} = $BBR{$elts[0]} + $elts[6];
	}
    }
    $k++;
    #print $k."/".scalar(@files)." processed\n";
    close(FHI);
}
#print "Writing results\n";

open(CLUST,">".$out_file_clusters) || die("cannot open file ".$out_file_clusters);
print CLUST $header;
foreach my $snp (keys(%AAT)){
       print CLUST $snp."\t";
       if($AAn{$snp} <= 0){
	    print CLUST "NA\t";
        }
        else{
	    print CLUST $AAT{$snp}/$AAn{$snp}."\t"; 
	}
        if($ABn{$snp} <= 0){
	    print CLUST "NA\t";
        }
        else{
	    print CLUST $ABT{$snp}/$ABn{$snp}."\t"; 
	}
        if($BBn{$snp} <= 0){
	    print CLUST "NA\t";
        }
        else{
	    print CLUST $BBT{$snp}/$BBn{$snp}."\t"; 
	}    
	if($AAn{$snp} <= 0){
	    print CLUST "NA\t";
        }
        else{
	    print CLUST $AAR{$snp}/$AAn{$snp}."\t"; 
	}
        if($ABn{$snp} <= 0){
	    print CLUST "NA\t";
        }
        else{
	    print CLUST $ABR{$snp}/$ABn{$snp}."\t"; 
	}
        if($BBn{$snp} <= 0){
	    print CLUST "NA\n";
        }
        else{
	    print CLUST $BBR{$snp}/$BBn{$snp}."\n"; 
	}
}
close(FHO);

if($file_list_samples ne ''){
    print "######### MERGING SAMPLE FILES ##########\n\n";
    @files=();
    open(FILES,$file_list_samples) || die("cannot open file ".$file_list_samples);
    while(my $line=<FILES>){
	$line =~ s/\n//g;
	push @files, $line;
    }
    close(FILES);
    my $index = 1;
    open(FHI,$files[0])|| die("cannot open file ".$files[0]);
    $header = <FHI>;
    close(FHI);
    open(SAMPLES,">".$out_file_samples) || die("cannot open file ".$out_file_samples);
    print SAMPLES $header;
    foreach my $f(@files){
	open(FHI,$f) || die("cannot open file ".$f);
	my $line = <FHI>;
	while($line=<FHI>){
	    my @elts = split(/\t/,$line);
	    print SAMPLES $elts[0]."\t".$elts[1]."\t".$index."\n";
	}
	$index++;
	close(FHI);
    }
    close(SAMPLES);
}
