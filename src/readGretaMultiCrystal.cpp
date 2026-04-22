#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <iomanip>
#include "Riostream.h"
#include <vector>
#include <signal.h>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <arpa/inet.h>
#include <map>
#include <ncurses.h>

#include <unistd.h>
#define ISATTY isatty
#define FILENO fileno

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "Globals.h"
#include "SortingStructures.h"


//#include "GRETA.h"

using namespace std;

#define CFD_INT_LEN 4
#define CFD_DELAY 4
#define CFD_FRACTION 4
#define AVG_TR_LENGTH 108
#define AVG_TR_STRIDE 110
#define TR_SCALE 10000

#define EB_DIFF_TIME 800

#define NUM_CHAN 40

#define DEBUG 0

/*****************************************************/
void PrintHelpInformation();

void GetData(FILE* inf, Int_t lastLength);
void SkipData(FILE *inf, UShort_t junk[]);
void LookForGoodData(FILE *inf, UShort_t junk[], Int_t lastEvt);
/*****************************************************/


/* Utility functions for interrupting the sort cleanly, and for terminal progress bar. */

Int_t gotSignal;
void breakhandler(int dummy) {
  cout << "Got break signal.  Aborting cleanly..." << endl;
  gotSignal = 1;
}

void progressB(int pct) {
  if (ISATTY(FILENO(stdout))) {
    string bar;
    struct winsize uk;
    if (ioctl(0, TIOCGWINSZ, &uk) != 0) {
      exit(1);
    }
    int wdt = uk.ws_col - 20;
    if (wdt < 5) { wdt = 5; }
    for (int i=0; i<wdt; i++) {
      if (i<(pct*wdt/100)) {
	bar.replace(i, 1, "=");
      } else if (i==(pct*wdt/100)) {
	bar.replace(i, 1, ">");
      } else {
	bar.replace(i, 1, " ");
      }
    }
    cout << "\r";
    cout << "[" << bar << "] "; 
    cout .width(3);
    cout << pct << "% complete" << flush;
  }
}

/* Function for 64-bit numbers endian-ness stuff */
uint64_t ntoh64(uint64_t input) {
  uint64_t rval;
  uint8_t *data = (uint8_t *)&rval;
  data[0] = input >> 56;
  data[1] = input >> 48;
  data[2] = input >> 40;
  data[3] = input >> 32;
  data[4] = input >> 24;
  data[5] = input >> 16;
  data[6] = input >> 8;
  data[7] = input >> 0;
  // printf("%x\n", rval);
  return rval;
}

uint64_t hton64(uint64_t input) {
  return (ntoh64(input));
}

int main (int argc, char *argv[]) {
  
  /* Some CTRL-C interrupt handling stuff... */
  gotSignal = 0;
  signal(SIGINT, breakhandler);

  controlVariables *ctrl = new controlVariables();
  ctrl->Initialize();
  
  if (argc < 3) {
    PrintHelpInformation(); exit(1);
  }
  if (ctrl->InterpretCommandLine(argc, argv) != 1) {
    exit(-1);
  }

  gret = new GRETA();
  gret->Initialize(ctrl);

  gret->rot.ReadMatrix("crmat.dat");
  
  //  gretaWaveformMsg wform;

  float ccEnergy[4];
  float segEnergy[36];
  float bl[40];
  vector<float> pz;
  float sumA = 0.; float sumB = 0.;
  int16_t trReOrder[NUM_CHAN][512]; 

  for (Int_t i=0; i<121; i++) {
    for (Int_t j=0; j<40; j++) {
      gret->gain[i][j] = 0.08; // Nominal for segments etc.
      gret->offset[i][j] = 0.0;
      if (j==37) {
	gret->gain[i][j] = 0.032;
      } else if (j==38) {
	gret->gain[i][j] = 0.016;
      } else if (j==39) {
	gret->gain[i][j] = 0.24;
      }
    }
  }

  if (ctrl->specifyCalibration) {
    gret->readGRETACalibration(ctrl->calibrationFile);
  }

  /* I think I don't need this anymore */  
  int mult = 0;
  float cfdTime = -1.;
  int data4net[36] = {0};
  float averageTrace[40][4096];
  int averageTraceINT[40][4096];
  float traceGain[40] = {0.0};
  
  /* Initialize ROOT output */
  TFile *fOut = new TFile(ctrl->rootFile.Data(), "RECREATE");
  TTree *gdata = new TTree("gdata", "gdata");
  if (!ctrl->calibrationRun) {
    gdata->Branch("g3", "g3OUT", &(gret->g3Out));
    gdata->Branch("g2", "g2OUT", &(gret->g2Out));
  } else {
    for (Int_t i=0; i<121; i++) {
      for (Int_t j=0; j<4; j++) {
	gret->ccSpec[i][j] = new TH1F(Form("xtal%dccSpec%d", i, j), Form("xtal%dccSpec%d", i, j), 10000, 0, 3000/gret->gain[i][j+36]);
      }
      for (Int_t j=0; j<36; j++) {
	gret->segSpec[i][j] = new TH1F(Form("xtal%dsegSpec%d", i, j), Form("xtal%dsegSpec%d", i, j), 10000, 0, 3000/gret->gain[i][j]);
      }
    }
  }
    
  unsigned char buf[65536];
  UShort_t junk[8192];

  /* Check file length for progress */
  struct stat fileStatus;
  Int_t multParts = 0;
  if (stat(ctrl->inputFile.Data(), &fileStatus) != 0) {
    cout << "Seems this is a file that is written in pieces.  Scanning all pieces... \n";
    multParts = 1;
  } 

  Int_t nFiles = 0;
  vector<Int_t> runNumList;
  vector<Int_t> runNumListEndNum;
  
  if (multParts && ctrl->holeNum >= 1) { // One Quad Only
    nFiles = 4;
    runNumList.push_back(ctrl->holeNum*4);
    runNumList.push_back(ctrl->holeNum*4 - 1);
    runNumList.push_back(ctrl->holeNum*4 - 2);
    runNumList.push_back(ctrl->holeNum*4 - 3);
  }
    
  int64_t bytesInFile = 0;

  if (multParts == 0) {
    bytesInFile = (int64_t) fileStatus.st_size;
    cout << "Input data file size is " << (Float_t)bytesInFile/1024./1024./1024. << " GB\n\n";
  } else if (multParts == 1) {
    if (nFiles > 0) {
      for (Int_t n=0; n<nFiles; n++) {
	for (Int_t endNum = 1; endNum<1000; endNum++) {
	  TString fName = ctrl->inputFile;
	  fName += Form("%d_%d", runNumList[n], endNum);	
	  if (stat(fName.Data(), &fileStatus) == 0) {
	    runNumListEndNum.push_back(endNum);
	    bytesInFile += (int64_t)fileStatus.st_size;
	    endNum = 1000;
	    printf("%s\n", fName.Data());
	  }
	}
      }
    } else {
      for (Int_t n=1; n<121; n++) {
	for (Int_t endNum = 1; endNum<1000; endNum++) {
	  TString fName = ctrl->inputFile;
	  fName += Form("%d_%d", n, endNum);	
	  if (stat(fName.Data(), &fileStatus) == 0) {
	    runNumList.push_back(n);
	    runNumListEndNum.push_back(endNum);
	    bytesInFile += (int64_t)fileStatus.st_size;
	    nFiles++;
	    endNum = 1000;
	    printf("%s\n", fName.Data());
	  }
	}
      }
    }
    cout << "Total input files (" << nFiles << " of them) size is " << (Float_t)bytesInFile/1024./1024./1024. << "GB\n\n";
  }
  
  int64_t bytesRead = 0;

  FILE *inf;

  Int_t nRuns = 1;
  if (multParts==1) { nRuns = nFiles; }
  
  long long int firstTS = 0;
  
  for (Int_t mm = 0; mm<nRuns; mm++) {
    if (nRuns == 1) {
      inf = fopen(ctrl->inputFile.Data(), "r");
      printf("\nOpened file - %s\n", ctrl->inputFile.Data());
    } else {
      TString fName = ctrl->inputFile;
      fName += Form("%d_%d", runNumList[mm], runNumListEndNum[mm]);
      inf = fopen(fName.Data(), "r");
      printf("\nOpened file - %s\n", fName.Data());
    }
        
    Int_t siz = 0;
    siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
    bytesRead += sizeof(struct routingHdr);
    
    long long int lastTS = 0;  long long int currTS = 0;
    long long int deltaEvent = 0;

    Int_t lastEvtLength = 0;
    
    Int_t TSreports = 0;
    
    while (siz && !gotSignal) {
      rHeader.seqnum = ntohs(rHeader.seqnum);
      rHeader.timestamp = (int64_t)ntoh64((uint64_t)rHeader.timestamp);
      if (firstTS == 0) { firstTS = rHeader.timestamp; }
      rHeader.length = ntohs(rHeader.length);
      rHeader.checksum = ntohs(rHeader.checksum);
      
      if (DEBUG) {
	printf("Routing Header: \n");
	printf("    Version:  0x%x\n", rHeader.version);
	printf("    Flags:    %i\n", rHeader.flags);
	printf("    Type:     %i\n", rHeader.type);
	printf("    SubType:  %i\n", rHeader.subtype);
	printf("    Length:   %i\n", rHeader.length);
	printf("    SeqNum:   %i\n", rHeader.seqnum);
	printf("    TS:       %lld\n", rHeader.timestamp);
	printf("    Checksum: %lli\n", rHeader.checksum);
	cin.get();
      } 
      
      if (rHeader.timestamp < lastTS && TSreports<10) {
	printf("TS out of order: last TS %lld, current %lld\n", lastTS, rHeader.timestamp);
	TSreports++;
      }
      lastTS = rHeader.timestamp;
      
      if ((rHeader.timestamp != 0) && (currTS == 0) ) {
	currTS = rHeader.timestamp;
      }
      
      deltaEvent = (Float_t)(rHeader.timestamp - currTS);
      if (DEBUG) {
	printf(" DeltaT from last TS: %lld\n", deltaEvent);
      }
      
      if (abs(deltaEvent) < EB_DIFF_TIME) {
	
	GetData(inf, lastEvtLength);
	//printf("Returned from GetData\n");
	bytesRead += rHeader.length;
	
      } else {
	
	if (!ctrl->calibrationRun) { gdata->Fill(); }
	gret->Reset();
	currTS = rHeader.timestamp;
	deltaEvent = (Float_t)(rHeader.timestamp - currTS);
	
	GetData(inf, lastEvtLength);
	bytesRead += rHeader.length;
      }

      lastEvtLength = rHeader.length;
      
      if (bytesRead%2000 == 0) { progressB(100*bytesRead/bytesInFile); }    
      siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
      bytesRead += sizeof(struct routingHdr);
      
    }
  } /* Loop over run segments */
  
  if (!ctrl->calibrationRun) { gdata->Write(); }
  if (ctrl->calibrationRun) {
    for (Int_t i=1; i<121; i++) {
      for (Int_t j=0; j<4; j++) {
      	gret->ccSpec[i][j]->Write();
      }
      for (Int_t j=0; j<36; j++) {
	gret->segSpec[i][j]->Write();
      }
    }
  }
  fOut->Write();
  fOut->Close();

  //std::cout << "\n --> nG2 " << gret->ng2 << std::endl;
  //long long int delta = (rHeader.timestamp-firstTS);
  //Float_t lengthInS = (Float_t)delta / 100000000.;
  //std::cout << "\n  Duration = " << Float_t(delta)/100000000. << " s" << std::endl;
  //std::cout << "\n Rate = " << (Float_t)gret->ng2/lengthInS << endl; 

  printf("\nDone and done.\n\n");
  return 0;
}

void GetData(FILE* inf, Int_t lastLength) {
  UShort_t junk[8192];

  // printf("In GetData - %d\n", rHeader.type);

  switch (rHeader.type) {
  case 3:
    {
      gret->getMode3(inf, rHeader.length, rHeader.subtype, rHeader.type);
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhTS = rHeader.timestamp;
    }
    break;
  case 4:
    {
      gret->getMode3(inf, rHeader.length, rHeader.subtype, rHeader.type);
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhTS = rHeader.timestamp;
    }
    break;
  case 2:
    {
      gret->getMode2(inf, rHeader.length, rHeader.subtype);
      gret->ng2++;
      gret->g2Out.xtals[gret->g2Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g2Out.xtals[gret->g2Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g2Out.xtals[gret->g2Out.crystalMult()-1].rhTS = rHeader.timestamp;
    }
    break;
  default:
    {
      cout << "Routing Header type not recognized: " << rHeader.type << endl;
      // SkipData(inf, junk);
      LookForGoodData(inf, junk,lastLength);
    }
  }

}

void SkipData(FILE *inf, UShort_t junk[]) {
  Int_t siz = fread(&junk, 1, ntohs(rHeader.length), inf);
  if (siz != ntohs(rHeader.length)) {
    cout << "SkipData(): Failed.\n";
    cout << endl;
  }
}

void LookForGoodData(FILE *inf, UShort_t junk[], Int_t lastEvt) {
  printf("Last Event Length - %d\n", lastEvt);
  
  fseek(inf, -1*(sizeof(struct routingHdr)), SEEK_CUR);
  Int_t success = 0;
  
  for (Int_t i=0; i<lastEvt*4; i++) {
    if (!success) {
      Int_t siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
      rHeader.seqnum = ntohs(rHeader.seqnum);
      rHeader.timestamp = (int64_t)ntoh64((uint64_t)rHeader.timestamp);
      rHeader.length = ntohs(rHeader.length);
      rHeader.checksum = ntohs(rHeader.checksum);
      //printf("%d---", i);
      if (DEBUG) {
	printf("Routing Header in Look4GoodData: \n");
	printf("    Version:  0x%x\n", rHeader.version);
	printf("    Flags:    %i\n", rHeader.flags);
	printf("    Type:     %i\n", rHeader.type);
	printf("    SubType:  %i\n", rHeader.subtype);
	printf("    Length:   %i\n", rHeader.length);
	printf("    SeqNum:   %i\n", rHeader.seqnum);
	printf("    TS:       %lld\n", rHeader.timestamp);
	printf("    Checksum: %lli\n", rHeader.checksum);
      }
      if (rHeader.version == 2 && rHeader.flags == 0 && rHeader.type == 4) {
	printf("Recovered after bad header - %d\n", i);
	GetData(inf, lastEvt);
	success = 1;
      } else {
	fseek(inf, -(sizeof(struct routingHdr)-1), SEEK_CUR);
      }
    }
  }
  if (!success) {  printf("Failed to find the start of another good event... just giving up now.\n"); }
  // cin.get();
}

void PrintHelpInformation() {
  printf("\n");
  printf("Usage: readGreta <Usage Flags> -f <InputFileWithPath> -rootFile <ROOTOutputName>\n");
  printf("     Valid usage flags: -readCal <calFileName>\n");
  printf("                        -hole <X=1 to 120>  (for analyzing the files from a single Quad only if data is taken crystal-wise)");
  printf("\n");
}

