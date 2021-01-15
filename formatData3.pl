#!/usr/local/bin/perl


main();

sub main {

$norm_dir = $ARGV[0];
$genotype_index = $ARGV[1];
$name_index = $ARGV[2];
$genoFilesList = $ARGV[3];
#$snpListFile = $ARGV[3]; 
$outdir = $ARGV[4];

mkdir $outdir;

@datafiles = ();
$nfiles = 0;
open(FHI,$genoFilesList) || die ("cannot open file ".$genoFilesList);
while($line=<FHI>){
  $line =~ s/\012//;
  $line =~ s/\015//;
  push @datafiles, $line;  
  $nfiles ++;
}
close(FHI);
print $nfiles." data files to format\n\n";

$k=0;
$step=10;
@indivs = ();
for $fileIn (@datafiles){
  #print "  - Formatting ".$fileIn." : ".$k."/".$nfiles." \n";
  %geno_data = ();
  @elts = split(/\//,$fileIn);
  $iid = $elts[scalar(@elts)-1];
  $iid =~ s/\.txt//g;
  $iid =~ s/_extracted//g;
  push @indivs, $iid;  
  mkdir $outdir."/".$iid;
  mkdir $outdir."/".$iid."/rawData";
  open(FHI,$fileIn)|| die ("cannot open file ".$fileIn);
  $line = <FHI>;
  while($line=<FHI>){
    $line =~ s/\012//;
    $line =~ s/\015//;
    @elts = split(/\t+/,$line);
    $geno_data{$elts[$name_index]}=$elts[$genotype_index];   
  }
  close(FHI);
  open(FHO,"> ".$outdir."/".$iid."/rawData/".$iid.".txt")|| die ("cannot open file ".$outdir."/".$iid."/rawData/".$iid.".txt");;  
  open (FHI, $norm_dir."/".$iid."_tQN_PennCNV.txt.adjusted")|| die ("cannot open file ". $norm_dir."/".$iid."_tQN_PennCNV.txt.adjusted");
  $line = <FHI>;
  print FHO "Name\tChr\tPosition\tLog R Ratio\tB Allele Freq\tGType\n";
  while($line = <FHI>){
    $line =~ s/\012//;
    $line =~ s/\015//;
    @elts = split(/\t+/,$line);
    if($geno_data{$elts[0]} ne ""){
      print FHO $elts[0]."\t".$elts[1]."\t".$elts[2]."\t".$elts[4]."\t".$elts[3]."\t".$geno_data{$elts[0]}."\n";
    }
  }
  close(FHO);
  close(FHI);
  $k++;
  $pct = $k/scalar(@indivs)*100;
  if($pct >= $step){
      print "\t".$step." % done\n";
      $step=$step+10;
  } 
}
close(FHO1);
print "Analysis finished\n";
}
