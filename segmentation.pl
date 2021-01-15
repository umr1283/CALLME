#!/usr/local/bin/perl

#use LWP::UserAgent;
#use DBI;
use Getopt::Long;
use Pod::Usage;


$absoluteDirPath = '';
$fileIndivs = ''; # Input file
$maxProc = 50;
$maxTime = 20;
$Tcutoff = 9;
$minSegLength = 75;
$aAlpha = 0.8;
$sdBAFthreshold = 1000;
$sdLRRthreshold = 1000;
$libDir="";

GetOptions('absolute_dir_path=s' => \$absoluteDirPath,
	   'indivs_file=s' => \$fileIndivs,
	   'max_procs=s' => \$maxProc,
	   'T=s' => \$Tcutoff,   
           'min_seg_length=s' => \$minSegLength,
           'aAlpha=s' => \$aAlpha,
           'sd_LRR_thres=s' => \$sdLRRthreshold,
           'sd_het_BAF_thres=s' => \$sdBAFthreshold,
	   'lib_dir=s' => \$libDir
            )or pod2usage(1);

$scriptSeg = $libDir."/launchSegmentation.R";
$scriptExtract = $libDir."/extractResults.R";

mkdir $absoluteDirPath."/logs";
$nprocs = 0;
@indivs = ();
$k=0;
$step=10;
open (FHI, $fileIndivs)|| die ("cannot open file ".$fileIndivs);
while($line = <FHI>){
  $line =~ s/\012//;
  $line =~ s/\015//;
  push @indivs, $line;
}


for $ind (@indivs){
   $cmd = "nohup R --slave dir=".$absoluteDirPath."/".$ind." analysis=MADseg Tcutoff=".$Tcutoff." minSegLength=".$minSegLength." aAlpha=".$aAlpha." < ".$scriptSeg." >& ".$absoluteDirPath."/logs/segmentation_".$ind.".out &";
  if($nprocs<$maxProc){
    system($cmd);
    $nprocs++;
    $pct = $k/scalar(@indivs)*100;
    if($pct >= $step){
	print "\t".$step." % done\n";
	$step=$step+10;
    }	
    $k++;
  }
  else{ 
    while($nprocs>=$maxProc){
      system("ps ux > tmpfile");
      open(FHps,"tmpfile") || die("cannot open file tmpfile");
      $lin = <FHps>;
      $nprocs = 0;
      while($lin=<FHps>){
        $lin =~ s/\012//;
        $lin =~ s/\015//;
        if($lin =~ m/analysis\=MADseg/gi){
            $nprocs ++; 
        }
      }
      close(FHps);
    }
    system($cmd);
    $pct = $k/scalar(@indivs)*100;
    if($pct >= $step){
	print "\t".$step." % done\n";
	$step = $step+10;
    }	
    $k++;
  }
}
close(FHI);

print "\n\tWaiting for all processes to finish ... \n\n";
$nprocs = 1;
while($nprocs > 0){
system("ps ux > tmpfile");
open(FHps,"tmpfile") || die("cannot open file tmpfile");
   $lin = <FHps>;
   $nprocs = 0;
   while($lin=<FHps>){
     $lin =~ s/\012//;
     $lin =~ s/\015//;
     if($lin =~ m/analysis\=MADseg/gi){
          $nprocs ++; 
     }
   }
   close(FHps);
}

system("rm -f tmpfile");


print "\n\tExtracting results ... \n"; 
$cmd = "R --slave dir=".$absoluteDirPath." indivsFile=".$fileIndivs." T=".$Tcutoff." minSegLength=".$minSegLength." aAlpha=".$aAlpha." sdBAFthreshold=".$sdBAFthreshold." sdLRRthreshold=".$sdLRRthreshold."< ".$scriptExtract." > ".$absoluteDirPath."/logs/segmentation_extraction.out ";
system($cmd);

