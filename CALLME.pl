#!/usr/bin/perl 

use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use POSIX;

# initialize user-definable parameters to default values

$dataFile='';
$nMarkers=0;
$rsOnly='false';
$cluster_centers='true';
$chip='';
$nprocs=1;
$T=2;
$aAlpha=0.2;
$lenProbe=75;
$libDir='.';
$sdHetBAF = 0.05;
$sdLRR = 0.33;
$remove='true';
# set user defined parameters
GetOptions('data_file=s' => \$dataFile,
	   'n_markers=s' => \$nMarkers,
	   'rs_only=s' => \$rsOnly,
	   'cluster_centers=s' => \$cluster_centers,
	   'chip=s' => \$chip,
	   'n_procs=s' => \$nprocs,
	   'T=s' => \$T,
	   'aAlpha=s' => \$aAlpha,
	   'min_seg_length=s' => \$lenProbe,
	   'sd_het_BAF_thres=s' => \$sdHetBAF,
	   'sd_LRR_thres=s' => \$sdLRR,
	   'remove_temp_files=s' => \$remove
           ) or pod2usage(1);

if (! -r $dataFile){
    die("cannot read ".$dataFile." (--data_file)");
}
if($nMarkers==0){
    die("you did not specify the maximum number of markers (--n_markers)");
}
if($chip eq '' & $cluster_centers ne "true"){
    die("you did not specify the chip (--chip)");
}


$dir=abs_path($dataFile);
@sp=split(/\//,$dir);
$filename=$sp[scalar(@sp)-1];
$dir =~ s/\/$filename//;
$libDir=abs_path($libDir);



print "############### ANALYSIS BEGAN ".localtime()." ##############\n\n";


system("top -b -n1  | head > ".$dir."/topfile");

open(FHI,$dir."/topfile") || die("cannot open file ".$dir."/topfile");
$line = <FHI>;
$line = <FHI>;
$line = <FHI>;
$line = <FHI>;
@elts = split(/, /,$line);
$mem_available = $elts[2];
$mem_available =~ s/k free//g;
$mem_available = POSIX::ceil($mem_available/1000);
close(FHI);
unlink $dir."/topfile";

system("head -n1 ".$dataFile." > ".$dir."/header");
open(FHI,$dir."/header") || die("cannot open file ".$dir."/header");
$line = <FHI>;
close(FHI);
unlink $dir."/header";


@elts = split(/\t/,$line);
$nInd = (scalar(@elts) - 3)/3;
#$nInd=11900;
$mem_required = 42950.1764636666 + 0.0625582763528102*$nInd*$nMarkers;
$mem_required = POSIX::ceil($mem_required*(1.1)/1000);


$option = 1;
print "\n\n******* ".$nInd." samples to process *******\n\n"; 
print "\t".$mem_available." MB of RAM available\n";
print "\tApproximately ".$mem_required." MB of RAM needed for analysis\n\n";

if($mem_required > $mem_available){
    print "\tWARNING : RAM required is higher than 80% of RAM available\n";
    print "\t\tWhat do you want to do? \n";
    print "\t\t\t - Still try to stock all data in the RAM (stock)\n";
    print "\t\t\t - Use a slower but less RAM-consuming program (less)\n";
    print "\t\t\t - Quit (quit)\n";
    $answer=<>;
    if($answer =~ m/quit/gi){
	die();
    }
    if($answer =~ m/less/gi){
	$option=2;
    }	
}


print "******* Extracting data *******\n\n";


#print $libDir."/split_samples.exe --data ".$dataFile." --nInd ".$nInd." --nSNPmax ".$nMarkers." --rsOnly ".$rsOnly." --nBlock ".$nprocs." --outDir ".$dir."\n";
if($option==1){
    $myCommandLY = $libDir."/split_samples.exe --data ".$dataFile." --nInd ".$nInd." --nSNPmax ".$nMarkers." --rsOnly ".$rsOnly." --nBlock ".$nprocs." --outDir ".$dir;
    #print "--> $myCommandLY";
    $exit_status=system($myCommandLY);
    if($exit_status != 0){
    	die("An error occured during data extraction");
    }
}
else{
    $finished=0;
    $start=0;
    $stop=250;
    if($stop>=$nInd){
	$stop=$nInd;
	$finished++;
    }
    while($stop<=$nInd & $finished<2){
	 print"perl ".$libDir."/split_samples.pl --data_file=".$dataFile." --output_directory=".$dir." --rs_only=".$rsOnly." --samples=".$start."-".$stop."\n";
	 system("perl ".$libDir."/split_samples.pl --data_file=".$dataFile." --output_directory=".$dir." --rs_only=".$rsOnly." --samples=".$start."-".$stop);
	 if($exit_status != 0){
	     die("An error occured during data extraction");
	 }
	 $start=$stop;
	 $stop=$start+250;
	 if($stop>=$nInd){
	     $stop=$nInd;
	     $finished++
	 }	 
    }
}

#question();
$nfilesPerList = POSIX::ceil($nInd/$nprocs);
system("R --vanilla --no-save --slave dir=".$dir."/normalization nFilesPerList=".$nfilesPerList." < ".$libDir."/split_samples.R > tmp");
unlink "tmp";

if($cluster_centers eq "true"){
    print "\n\n******* Computing cluster centers *******\n\n";
    mkdir $dir."/cluster_centers";
    system("perl ".$libDir."/cluster_centers.pl --input_directory=".$dir."/normalization --geno_directory=".$dir."/genotypes --output_file=".$dir."/cluster_centers/clusters --samples_lists=".$dir."/normalization/lists.txt --lib_dir=".$libDir." --log_dir=".$dir."/cluster_centers/logs");
    waitFor('analysis=cluster_centers');
    checkLogs($dir."/cluster_centers/logs","cluster");
    open(FHO,"> ".$dir."/cluster_centers/list.txt") || die("cannot open file ".$dir."/cluster_centers/list.txt for writing");
    for($i=1;$i<=$nprocs;$i++){
	print FHO $dir."/cluster_centers/clusters".$i."\n";
    }
    close(FHO);
    if(! -d $libDir."/lib_tQN"){
	mkdir  $libDir."/lib_tQN";
    }
    system("perl ".$libDir."/join_cluster_and_samples_files.pl --file_list_clusters=".$dir."/cluster_centers/list.txt --out_file_clusters=".$libDir."/lib_tQN/my_chip_tQN_clusters.txt");
    print "\n => custom cluster centers are available here :".$libDir."/lib_tQN/my_chip_tQN_clusters.txt\n";
    $chip="my_chip";
    if($remove =~ m/true/gi){
	system("rm -rf ".$dir."/cluster_centers");
    }
}

print "\n\n******* Performing quantile normalization *******\n\n";

system("perl ".$libDir."/tQN_normalize_samples.pl --beadchip=".$chip." --input_dir=".$dir."/normalization --output_dir=".$dir."/normalized --samples_lists=".$dir."/normalization/lists.txt --all_samples=".$dir."/normalization/sample_names.txt --tQN_dir=".$libDir."/lib_tQN");
checkLogs($dir."/normalized/logs","normalization");
checkLogs($dir."/normalized/logs","generate");

#question();

if($remove =~ m/true/gi){
    system("rm -rf ".$dir."/normalization");
}

print "\n\n******* Performing GC wave correction *******\n\n";

system("R --vanilla --no-save --slave dir=".$dir."/normalized  nFilesPerList=".$nfilesPerList." libDir=".$libDir." < ".$libDir."/launch_genomic_wave.R > tmp");
waitFor("genomic_wave.pl");


checkLogs($dir."/normalized/GC_corrected/logs","correction");
unlink "tmp";

#question();

print "\n\n******* Formatting data for segmentation *******\n\n";

system("R --vanilla --no-save -slave dir1=".$dir."/normalized/GC_corrected dir2=".$dir."/genotypes nFilesPerList=".$nfilesPerList." outdir=".$dir."/segmentation < ".$libDir."/launchFormatting.R > tmp");
waitFor("formatData3.pl");
checkLogs($dir."/segmentation/logs","formatData");
unlink "tmp";
unlink $dir."/segmentation/indivs.txt";
if($remove =~ m/true/gi){
    system("rm -rf ".$dir."/normalized");
    system("rm -rf ".$dir."/genotypes");
}
opendir(DIR, $dir."/segmentation") || die("Cannot open directory ".$dir."/segmentation"); 
@indivs = readdir(DIR);
close(DIR);

open(FHO,"> ".$dir."/segmentation/indivs.txt") || die("cannot open file ".$dir."/segmentation/indivs.txt");
for $ind (@indivs){
    if($ind ne "logs" & $ind ne "." & $ind ne ".." & $ind !~ m/res/g & $ind !~ m/data/g){
	print FHO $ind."\n";
    }
}
close(FHO);

if($remove =~ m/true/gi){
    system("rm -rf ".$dir."/normalized");
    system("rm -rf ".$dir."/genotypes");    
}

print "\n\n******* Calculating quality metrics on data  *******\n\n";

system("perl ".$libDir."/checkData.pl --absolute_dir_path=".$dir."/segmentation/ --indivs_file=".$dir."/segmentation/indivs.txt --max_procs=".$nprocs." --lib_dir=".$libDir);
checkLogs($dir."/segmentation/logs","check");


print "\n\n******* Performing segmentation  *******\n\n";

system("perl ".$libDir."/segmentation.pl --absolute_dir_path=".$dir."/segmentation/ --indivs_file=".$dir."/segmentation/indivs.txt --max_procs=".$nprocs." --T=".$T." --min_seg_length=".$lenProbe." --aAlpha=".$aAlpha." --sd_LRR_thres=".$sdLRR." --sd_het_BAF_thres=".$sdHetBAF." --lib_dir=".$libDir);

checkLogs($dir."/segmentation/logs","segmentation");
checkLogs($dir."/segmentation/logs","extraction");

system("cp -r ".$dir."/segmentation/results_T_".$T."_aAlpha_".$aAlpha."_minSegLength_".$lenProbe." ".$dir) ;
system("cp ".$dir."/segmentation/data_check_all.txt ".$dir) ;
system("cp ".$dir."/segmentation/indivs.txt ".$dir) ;

print "\n\n############### ANALYSIS FINISHED ".localtime()." ##############\n\n";


sub waitFor(){
    $keyword = $_[0];
    $nproc = 1;
    while($nproc > 0){
	system("ps ux > tmpfile");
	open(FHps,"tmpfile") || die("cannot open file tmpfile");
	$lin = <FHps>;
	$nproc = 0;
	while($lin=<FHps>){
	    chomp($lin);
	    if($lin =~ m/$keyword/g){
	       $nproc ++; 
	    }
	}
	#print $nprocs."\n";
	close(FHps);
	unlink "tmpfile";
    }
}
sub checkLogs(){
    $logdir = @_[0];
    $keyword = @_[1];
    opendir(DIR, $logdir) || die("Cannot open directory ".$logdir); 
    @logs = readdir(DIR);
    close(DIR);
    for $log (@logs){
	if($log ne '.' & $log ne '..' & $log =~ m/$keyword/g){
	    #print $log."\n";
	    open(LOG,$logdir."/".$log) || die("Cannot open file ".$log);
	    $OK=0;	    
	    $line = <LOG>;
	    $content = $line;
	    while($OK eq 0 & $line ne ""){
		#print "ok : ".$OK."\n";
		if($line =~ m/Analysis finished/g){
		    $OK = 1;
		}
		$line=<LOG>;
                $content=$content.$line;
	    }
	    if($OK eq 0){
		print "NOTICE : An error occured  :\n";
		open $content."\n\n";
	    }
	    close(LOG);
	}
    }
}

sub question(){
    print "\tContinue (y or n) ?\n";
    $ans = <>;
    chomp($ans);
    if($ans ne "y"){
	die();
    }
}
