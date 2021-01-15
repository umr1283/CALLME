#!/usr/local/bin/perl

#use LWP::UserAgent;
#use DBI;
use Getopt::Long;
use Pod::Usage;

$absoluteDirPath = '';
$fileIndivs = ''; # Input file
$maxProc = 50;
$libDir = '';

GetOptions('absolute_dir_path=s' => \$absoluteDirPath,
	   'indivs_file=s' => \$fileIndivs,
	   'max_procs=s' => \$maxProc,
	   'lib_dir=s' => \$libDir
            )or pod2usage(1);

mkdir $absoluteDirPath."/logs";
$scriptCheck = $libDir."/checkData.R";
$scriptExtract = $libDir."/mergeDataChecks.R";

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
  $cmd = "nohup R --slave dir=".$absoluteDirPath."/".$ind." analysis=checkData < ".$scriptCheck." >& ".$absoluteDirPath."/logs/check_".$ind.".out &";
  if($nprocs<$maxProc){
    system($cmd);
    $nprocs++;
    $pct = $k/scalar(@indivs)*100;
    if($pct >= $step){
	print "\t".$step." % done\n";
	$step=$step+10;
    }   
    $k++;
    #print "Number of processes running : ".$nprocs."\n";    
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
    	      if($lin =~ m/analysis\=checkData/gi){
		  $nprocs ++;   
	      }
	  }
	  close(FHps);
	  #print "Number of processes running : ".$nprocs."\n";     
      }
      system($cmd);
      $pct = $k/scalar(@indivs)*100;
      if($pct >= $step){
	  print "\t".$step." % done\n";
	  $step=$step+10;
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
     if($lin =~ m/analysis\=checkData/gi){
          $nprocs ++; 
     }
   }
   close(FHps);
}

system("rm -f tmpfile");


print "\n\tExtracting results ... \n"; 
$cmd = "R --slave dir=".$absoluteDirPath." indivsFile=".$fileIndivs." < ".$scriptExtract." > ".$absoluteDirPath."/logs/check_extraction.out ";
system($cmd);


