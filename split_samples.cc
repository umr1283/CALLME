#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
using namespace std;
string toString(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

typedef struct{
  int  *Inds;
  int  nInds;
  int   part;
  int   nSNP;
  int   nInd;
  float  *X;
  float  *Y;
  string *GEN;
  string *indivNames;
  string *SNP;
  string *CHR;
  string *POS;
  char   *tmpDir;
} parm;


string* getMemCPU(){
  string *res = new string[2];
  int i;
  string psfile;
  char *PSfile;  
  string line;  
  char current[20]; 
  i=system("ps ux | grep nSNPmax > psfiletemp");
  psfile="psfiletemp";
  PSfile = new char[psfile.size()+1];
  PSfile[psfile.size()]=0;
  memcpy(PSfile,psfile.c_str(),psfile.size());
  ifstream yx (PSfile,ios::in);
  if(!yx){
    cerr<<"\tCannot open the file!"<<psfile<<endl;
    exit(1);    
  }else{
    //while(!yx.eof()){
      getline(yx,line);
      stringstream ss(line);
      ss>>current;
      //cout<<current<<"\t";
      ss>>current;
      //cout<<current<<"\t";
      ss>>current;
      //cout<<current<<"\t";
      ss>>current;
      //cout<<current<<"\t";
      res[0]=current;
      ss>>current;
      //cout<<current<<endl;
      res[1]=current;
      //}
  }
  i=system("rm -f psfiletemp");
  return res;
}
void writeFiles(int *Inds, int nInds, int part,int nSNP,int nInd,
		float  *X, float  *Y, string *GEN, string *indivNames,
		string *SNP, string *CHR, string *POS, char *tmpDir){
  int i, indiv;
  int iSNP;
  string FileGeno, FileXY;
  char *filegeno, *filexy;
  ofstream outFile;
  string base = tmpDir;
  //cout<<"\tEntering writeFiles"<<endl;
  for(i=0;i<nInds;i++){
    indiv = Inds[i]-1;
    // Write Genotypes
    FileGeno = base + "/genotypes/"+indivNames[indiv] + "_extracted.txt";
    filegeno = new char[FileGeno.size()+1];
    filegeno[FileGeno.size()]=0;
    memcpy(filegeno,FileGeno.c_str(),FileGeno.size());
    outFile.open (filegeno); 
    outFile<<"Name\tChr\tPosition\tGType"<<endl;
    for(iSNP=0;iSNP<nSNP;iSNP++){
      outFile<<SNP[iSNP]<<"\t"<<CHR[iSNP]<<"\t"<<POS[iSNP]<<"\t"<<GEN[indiv+nInd*iSNP]<<endl;
    }
    outFile.close();
    // Write XY
    FileXY = base + "/normalization/" + indivNames[indiv] + "_extracted.txt";
    filexy = new char[FileXY.size()+1];
    filexy[FileXY.size()]=0;
    memcpy(filexy,FileXY.c_str(),FileXY.size());
    outFile.open (filexy); 
    outFile<<"Name\tChr\tPosition\tX\tY"<<endl;
    for(iSNP=0;iSNP<nSNP;iSNP++){
      outFile<<SNP[iSNP]<<"\t"<<CHR[iSNP]<<"\t"<<POS[iSNP]<<"\t"<<X[indiv+nInd*iSNP]<<"\t"<<Y[indiv+nInd*iSNP]<<endl;
    }
    outFile.close();    
  }

}
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}
int makeDir(string dir){
    char *dirchar = new char[dir.size()+1];
    dirchar[dir.size()]=0;
    memcpy(dirchar,dir.c_str(),dir.size());
    int res=0;
    if(access(dir.c_str(),0)!=0){
      res=mkdir(dirchar,0777);
    }
    return res;
}

void * worker(void *arg){
  parm           *p = (parm *) arg;
  writeFiles(p->Inds,p->nInds, p->part, p->nSNP,p->nInd,p->X, p->Y, p->GEN,p->indivNames,p->SNP,p->CHR,p->POS,p->tmpDir);  
}


int main(int argc, char *argv[]){
  cout<<"\tAnalysis started "<<currentDateTime()<<endl;
  string datafile   = "none";
  string directory  = ".";
  char *tmpDir;
  int nInd          =      0;
  int nSNPmax       =      0;
  int nBlock        =      8; 
  string rsOnly = "false";
  // Read arguments
    for (int i = 1; i <= argc; ) { 
    if (i + 1 == argc | argc <= 2) {	
      exit(1);
    }      
    else if (i + 1 < argc) {
      string sw = argv[i];
      if (sw == "--data"){
	datafile = argv[i + 1];
	i = i + 2;
      }            
      else if (sw == "--nInd") {
	nInd = atoi(argv[i + 1]);
	i = i + 2;
      }
      else if (sw == "--nSNPmax") {
	nSNPmax = atoi(argv[i + 1]);
	i = i + 2;
      }
      else if (sw == "--rsOnly") {
	rsOnly = argv[i + 1];
	i = i + 2;
      }
      else if (sw == "--nBlock") {
	nBlock = atoi(argv[i + 1]);
	i = i + 2;
      }
      else if (sw == "--outDir") {
	directory = argv[i + 1];
	tmpDir = new char[directory.size()+1];
	tmpDir[directory.size()]=0;
	memcpy(tmpDir,directory.c_str(),directory.size());
	i = i + 2;
      }
      else {
 	++i;
      }
    }
    else break;
  }
  if(datafile=="none"){
    cerr<<"\tAn input file must be specified\n"<<endl;
    exit(1);
  }
  if(nSNPmax==0){
    cerr<<"\tA maximum number of SNP must be specified\n"<<endl;
    exit(1);
  }
  if(nInd==0){
    cerr<<"\tAn sample size must be specified\n"<<endl;
    exit(1);
  }
  string base = tmpDir; 
  int n = makeDir(base);
  if(n!=0){
     cerr<<"\tCould not create "<<tmpDir<<"\n"<<endl;
     exit(1);
   }
  n=makeDir(base+"/genotypes");
  if(n!=0){
     cerr<<"\tCould not create "<<base<<"/genotypes \n"<<endl;
     exit(1);
  }  
  n=makeDir(base+"/normalization");
  if(n!=0){
    cerr<<"\tCould not create "<<base<<"/normalization \n"<<endl;
    exit(1);
  }

  int nThread;
  int part,ipart;
  int cond         = nInd%(nBlock);
  if(cond>0){
    nThread = nBlock + 1;
  }else{
    nThread = nBlock;
  }
  int *BlockLength = new int[nThread];
  int sizeBase = (nInd-cond)/nBlock;  
  int sizeComp = cond;
  if(cond>0){
    for(part = 0;part<nBlock;part++){
      BlockLength[part] = sizeBase;
    }
    BlockLength[nThread-1] = sizeComp;
  }else{
    for(part = 0;part<nThread;part++){
      BlockLength[part] = sizeBase;
    }
  }
  
  int **IndsBlock   = new int*[nThread];
  for(part = 0;part<nBlock;part++){
    IndsBlock[part] = new int[BlockLength[part]];
    for(ipart=0;ipart<BlockLength[part];ipart++){
      IndsBlock[part][ipart] = 1+ipart + part * BlockLength[part];
    }
  }
  if(cond>0){
    part = nThread - 1;
    IndsBlock[part] = new int[BlockLength[part]];
    for(ipart=0;ipart<BlockLength[part];ipart++){      
      IndsBlock[part][ipart] = 1 + ipart + IndsBlock[part-1][sizeBase-1];
    }
  }
  
  string *mem;  
  //mem=getMemCPU();
  //cout<<"\t\tMemory used before initialization:"<<mem[1]<<" KB ("<<mem[0]<<" %)"<<endl;
  // Build internal data
  string *SNP = new string[nSNPmax];
  string *POS = new string[nSNPmax];
  string *CHR = new string[nSNPmax];

  // Matrices
  string *GEN = new string[nInd * nSNPmax];
  float  *X   = new float [nInd * nSNPmax];
  float  *Y   = new float [nInd * nSNPmax];
  
  string tit;
  //mem=getMemCPU();
  //cout<<"\t\tMemory used after initialization:"<<mem[1]<<" KB ("<<mem[0]<<" %)"<<endl;
 

  char *filename = new char[datafile.size()+1];
  filename[datafile.size()]=0;
  memcpy(filename,datafile.c_str(),datafile.size());
  ifstream yx (filename,ios::in);
  
  string line;
  string title;
  string name;
  string *indivNames = new string[nInd];
  char fl[20];
  int count = 0;
  int step=1;
  int nSNP  = 0;
  string isrs;
  int i;
  int j;
  int alreadyFound;
  int nDuplicates=1;
  int index = 1;
  string fileSampleNames;
  char *filesamples;
  ofstream outFile;
  cout<<"\t\tReading data ..."<<endl;
  if(!yx){
    cerr<<"\tCannot open the file!"<<datafile<<endl;
    exit(1);    
  }else{
    if(!yx.eof()){
      getline(yx,title);
      //cout<<title;
      stringstream is(title);
      is>>name;
      is>>name;
      is>>name;
      fileSampleNames = base + "/normalization/sample_names.txt";
      filesamples = new char[fileSampleNames.size()+1];
      filesamples[fileSampleNames.size()]=0;
      memcpy(filesamples,fileSampleNames.c_str(),fileSampleNames.size());
      outFile.open (filesamples); 
      outFile<<"Assay\tFilename\tIGV_index"<<endl;
      for(i=0;i<nInd*3;i+=3){
	is>>name;
	name.replace(name.find("."),name.length(),"");
	alreadyFound=0;
	for(j=0;j<nInd;j++){
	  if(indivNames[j] == name){
	    alreadyFound=1;
	  }
	}
	if(alreadyFound==0){
	  indivNames[i/3] = name;
	}
	else{
	  indivNames[i/3] = name+"_duplicate"+toString(nDuplicates);
	  nDuplicates++;
	}
	outFile<<indivNames[i/3]<<"\t"<<base+"/normalization/"+indivNames[i/3]+"_extracted.txt\t"<<index<<endl;
	is>>name;
	is>>name;
	index++;
      }
      outFile.close();
      cout<<"\t\tNOTICE : "<<nDuplicates-1<<" duplicated sample names found"<<endl;
    }
    while(!yx.eof()){
      getline(yx,line);
      isrs = line.substr (0,2);
      count++;
      //cout<<"\tReading line "<<count<<endl;
      if(isrs=="rs" | rsOnly=="false" | rsOnly=="False" | rsOnly=="FALSE" ){
	stringstream is(line);	
	nSNP++;
	if(nSNP >= nSNPmax){
	   cerr<<"\tYou specified a number of SNPs that is inferior to the number of SNPs found in the file\n"<<endl;
	   exit(1);
	}
	is>>SNP[nSNP];
	is>>CHR[nSNP];
	is>>POS[nSNP];
	for(i=0;i<nInd;i++){
	  is>>fl; X[i+nInd*nSNP] = atof(fl);
	  is>>fl; Y[i+nInd*nSNP] = atof(fl);
	  is>>GEN[i+nInd*nSNP] ;
	}
	//if(count==step*30000){
	// step++;
	// mem=getMemCPU();
	// cout<<"\t\tMemory used after "<<count<<" line : "<<mem[1]<<" KB ("<<mem[0]<<" %)"<<endl;
	//}
      }
    }
  }
  // From command line...
  cout<<"\t\tDatafile: "<<datafile<<endl;
  cout<<"\t\tSample size: "<<nInd<<endl;
  cout<<"\t\tNumber of SNPs: "<<nSNP<<endl;
  cout<<"\t\ttmpDir = "<<tmpDir<<endl;
  //mem=getMemCPU();
  //cout<<"\t\tMemory used after reading data: "<<mem[1]<<" KB ("<<mem[0]<<" %)"<<endl;
  
  
  //   // Write files
  pthread_t      *threads;
  parm           *arg;
  threads = (pthread_t *) malloc(nThread * sizeof(pthread_t));
  arg=(parm *)malloc(sizeof(parm)*nThread);
  // Spawn thread 
  for (part = 0; part < nThread; part++){
    arg[part].Inds       = IndsBlock[part];
    arg[part].nInds      = BlockLength[part];
    arg[part].part       = part;
    arg[part].nSNP       = nSNP;
    arg[part].nInd       = nInd;
    arg[part].X          = X;
    arg[part].Y          = Y;
    arg[part].GEN        = GEN;
    arg[part].indivNames = indivNames;
    arg[part].SNP        = SNP;
    arg[part].CHR        = CHR;
    arg[part].POS        = POS;
    arg[part].tmpDir     = tmpDir;
    pthread_create(&threads[part], NULL, worker, (void *)(arg+part));
    //  cout<<"\tOK : "<<tmpDir<<endl;
  }
  for (part = 0; part < nThread; part++){
    pthread_join(threads[part], NULL);    
  }
  
  free(arg);

  // Free memory
  if(SNP){
    delete [] SNP;
  }
  if(POS){
    delete [] POS;
  }
  if(CHR){
    delete [] CHR;
  }
  if(GEN){
    delete [] GEN;
  }
  if(X){
    delete [] X;
  }
  if(Y){
    delete [] Y;
  }
  if(IndsBlock){
    delete [] IndsBlock;
  }
  cout<<"\tAnalysis finished "<<currentDateTime()<<endl;

};

