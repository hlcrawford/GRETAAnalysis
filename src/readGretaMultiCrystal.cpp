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
#include "TVector3.h"
#include "TClass.h"
#include "TSystem.h"

#include "Globals.h"
#include "SortingStructures.h"

//#include "GRETA.h"

using namespace std;

/** Write one coincidence to gdata: full pre-Doppler g2Out + coincTS (tracking is pass 2: trackGreta). */
static void flushCoincidenceToGdata(GRETA *gret, g2OUT *&g2OutBranch, Long64_t &coincTS, TTree *gdata)
{
  if (!g2OutBranch) g2OutBranch = new g2OUT();
  *g2OutBranch = gret->g2Out;
  coincTS = 0;
  if (g2OutBranch->crystalMult() > 0)
    coincTS = g2OutBranch->xtals[0].timestamp;
  gdata->Fill();
}

#define CFD_INT_LEN 4
#define CFD_DELAY 4
#define CFD_FRACTION 4
#define AVG_TR_LENGTH 108
#define AVG_TR_STRIDE 110
#define TR_SCALE 10000

#define EB_DIFF_TIME 210

#define NUM_CHAN 40

#define DEBUG 0

/*****************************************************/
void PrintHelpInformation();

void GetData(FILE* inf, Int_t lastLength, Int_t lastSeqNum);
void SkipData(FILE *inf, UShort_t junk[]);
void LookForGoodData(FILE *inf, UShort_t junk[], Int_t lastEvt, Int_t lastSeqNum);
/*****************************************************/


/* Utility functions for interrupting the sort cleanly, and for terminal progress bar. */

Int_t gotSignal;
void breakhandler(int dummy) {
  cout << "Got break signal.  Aborting cleanly..." << endl;
  gotSignal = 1;
}

void progressB(int pct) {
  string bar;
  struct winsize uk;
  if (ISATTY(FILENO(stdout))) {
    if (ioctl(0, TIOCGWINSZ, &uk) != 0) {
      return;
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

/** Defined in lib/GRETADict.cpp; registers all GRETA streamer dictionaries. */
extern void TriggerDictionaryInitialization_GRETADict();

static bool GretaDictionaryReady()
{
  return TClass::GetClass("g2OUT") && TClass::GetClass("vector<TVector3>");
}

/** Load libGRETA and GRETA dict (g2OUT, vector<TVector3>, etc.). */
static void EnsureGretaDictionaryLoaded()
{
  if (GretaDictionaryReady()) return;

  const char *pwd = gSystem->WorkingDirectory();
  const TString libDir = TString::Format("%s/lib", pwd);
  gSystem->AddDynamicPath(libDir.Data());

  TString inc = gSystem->Getenv("ROOT_INCLUDE_PATH");
  if (inc.Length()) inc += ":";
  inc += libDir;
  gSystem->Setenv("ROOT_INCLUDE_PATH", inc.Data());

  const TString libPaths[] = {
    libDir + "/libGRETA",
    TString("lib/libGRETA"),
    TString("libGRETA"),
  };
  for (const TString &lp : libPaths) {
    if (gSystem->Load(lp.Data()) >= 0) break;
  }

  TriggerDictionaryInitialization_GRETADict();

  if (!GretaDictionaryReady()) {
    fprintf(stderr,
	    "readGreta: GRETA ROOT dictionary not loaded (need g2OUT).\n"
	    "  Run from the TrackBackToMe/GRETAAnalysis directory after:  scons readGreta\n"
	    "  Expect: %s/libGRETA.dylib (or .so) and %s/GRETADict_rdict.pcm\n",
	    libDir.Data(), libDir.Data());
    exit(1);
  }
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

  EnsureGretaDictionaryLoaded();

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

  g2OUT *g2OutBranch = nullptr;
  Long64_t coincTS = 0;

  if (!ctrl->calibrationRun) {
    gdata->Branch("g2Out", "g2OUT", &g2OutBranch);
    gdata->Branch("coincTS", &coincTS);
    printf("gdata branches: g2Out (full g2OUT), coincTS  (tracking: use trackGreta pass 2)\n");
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
    } else if (0) {
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
    } else {
      for (Int_t endNum = 1; endNum<1000; endNum++) {
	TString fName = ctrl->inputFile;
	fName += Form("%d", endNum);
	if (stat(fName.Data(), &fileStatus) == 0) {
	  runNumListEndNum.push_back(endNum);
	  bytesInFile += (int64_t)fileStatus.st_size;
	  nFiles++;
	  printf("%s\n", fName.Data());
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

  Short_t currentPercentage = 0, lastPercentage = -1;
  Int_t TSreports = 0;
  
  for (Int_t mm = 0; mm<nRuns; mm++) {
    if (multParts==0 && nRuns == 1) {
      if (!ctrl->hfc) {
	inf = fopen(ctrl->inputFile.Data(), "r");
      } else {
	ctrl->inputFile = "./GEB_HFC -p " + ctrl->inputFile;
	inf = popen(ctrl->inputFile.Data(), "r");
      }
      printf("\nOpened file - %s\n", ctrl->inputFile.Data());
    } else if (0) {
      TString fName = ctrl->inputFile;
      fName += Form("%d_%d", runNumList[mm], runNumListEndNum[mm]);
      if (!ctrl->hfc) {
	inf = fopen(fName.Data(), "r");
      } else {
	fName = "./GEB_HFC -p " + fName;
	inf = popen(fName.Data(), "r");
      }
      printf("\nOpened file - %s\n", fName.Data());
    } else {
      TString fName = ctrl->inputFile;
      fName += Form("%d", runNumListEndNum[mm]);
      if (!ctrl->hfc) {
	inf = fopen(fName.Data(), "r");
      } else {
	fName = "./GEB_HFC -p " + fName;
	inf = popen(fName.Data(), "r");
      }
      printf("\nOpened file - %s\n", fName.Data());
    }
        
    Int_t siz = 0;
    siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
    bytesRead += sizeof(struct routingHdr);
    
    long long int lastTS = 0;  long long int currTS = 0;

    Int_t lastEvtLength = 0;
    Int_t lastSeqNum = 0;
    
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
	printf(" current type: %d\n", rHeader.subtype);
	TSreports++;
      } else if (rHeader.timestamp < lastTS) {
	TSreports++;
      }
      lastTS = rHeader.timestamp;
      
      if ((rHeader.timestamp != 0) && (currTS == 0) ) {
	currTS = rHeader.timestamp;
      }
      
      const long long dt_raw = rHeader.timestamp - currTS;
      const long long dt_abs = (dt_raw >= 0) ? dt_raw : -dt_raw;
      if (DEBUG) {
	printf(" DeltaT from last TS: %lld\n", dt_raw);
      }
      
      if (dt_abs < (long long)EB_DIFF_TIME) {
	
	GetData(inf, lastEvtLength, lastSeqNum);
	//printf("Returned from GetData\n");
	bytesRead += rHeader.length;
	
      } else {
	
	if (!ctrl->calibrationRun) {
	  flushCoincidenceToGdata(gret, g2OutBranch, coincTS, gdata);
	}
	if (ctrl->calibrationRun) { gret->fillHistos(); }
	gret->Reset();
	currTS = rHeader.timestamp;
	
	GetData(inf, lastEvtLength, lastSeqNum);
	bytesRead += rHeader.length;
      }

      lastEvtLength = rHeader.length;
      lastSeqNum = rHeader.seqnum;
      
      currentPercentage = 100*bytesRead/bytesInFile;
      if (currentPercentage != lastPercentage) { progressB(currentPercentage); lastPercentage = currentPercentage; }
      siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
      bytesRead += sizeof(struct routingHdr);
      
    }

    /* Flush last coincidence (EOF / signal exit leaves built event without boundary Fill). */
    if (gret->g2Out.crystalMult() > 0 || gret->g3Out.crystalMult() > 0) {
      if (!ctrl->calibrationRun) {
	flushCoincidenceToGdata(gret, g2OutBranch, coincTS, gdata);
      }
      if (ctrl->calibrationRun) {
	gret->fillHistos();
      }
      gret->Reset();
    }

  } /* Loop over run segments */
  
  fOut->cd();
  if (ctrl->calibrationRun) { gret->gHist.writeHistos(); }
  fOut->Write();
  fOut->Close();

  for (Int_t i=0; i<121; i++) {
    if (gret->eventCnt[i] > 0) {
      printf("\n Crystal %d - total events %d", i, gret->eventCnt[i]);
    }
    if (gret->mode3TooLong[i] > 0) {
      printf("\n --> Crystal %d - %d mode3 events over frame size", i, gret->mode3TooLong[i]);
    }
    if (gret->seqNumOOO[i] > 0) {
      printf("\n --> Crystal %d - %d sequence numbers out of order", i, gret->seqNumOOO[i]);
    }
    if (gret->timesSeqNumSkipped[i] > 0) {
      printf("\n --> Crystal %d - %d times sequence numbers skipped", i, gret->timesSeqNumSkipped[i]);
    }
    if (gret->seqNumSkipped[i] > 0) {
      printf("\n --> Crystal %d - %d sequence numbers skipped", i, gret->seqNumSkipped[i]);
    }
  }

  if (TSreports > 100000) {
    printf("\033[0;31mTimestamp errors %d\033[0m\n", TSreports);
  } else if (TSreports > 10000) {
    printf("\033[0;33mTimestamp errors %d\033[0m\n", TSreports);
  } else {
    printf("\033[0;32mTimestamp errors %d\033[0m\n", TSreports);
  }
    
  std::cout << "\n --> nG2 " << gret->ng2 << std::endl;
  //long long int delta = (rHeader.timestamp-firstTS);
  //Float_t lengthInS = (Float_t)delta / 100000000.;
  //std::cout << "\n  Duration = " << Float_t(delta)/100000000. << " s" << std::endl;
  //std::cout << "\n Rate = " << (Float_t)gret->ng2/lengthInS << endl; 

  printf("\nDone and done.\n\n");
  return 0;
}

void GetData(FILE* inf, Int_t lastLength, Int_t lastSeqNum) {
  UShort_t junk[8192];

  // printf("In GetData - %d\n", rHeader.type);

  switch (rHeader.type) {
  case 3:
    {
      gret->getMode3(inf, rHeader.length, rHeader.subtype, rHeader.type);
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhTS = rHeader.timestamp;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhLength = rHeader.length;
      if (gret->lastSeqNum[rHeader.subtype] == -1) {
	gret->lastSeqNum[rHeader.subtype] = rHeader.seqnum;
      } else {
	Int_t deltaSeqNum = rHeader.seqnum - gret->lastSeqNum[rHeader.subtype];
	if (deltaSeqNum < 1) {
	  printf("Crystal %d - sequence numbers out of order.  Previous %d, now %d\n", rHeader.subtype,
		 gret->lastSeqNum[rHeader.subtype], rHeader.seqnum);
	  gret->seqNumOOO[rHeader.subtype]++;
	} else if (deltaSeqNum > 1) {
	  printf("Crystal %d - sequence numbers skipped.  Previous %d, now %d\n", rHeader.subtype,
		 gret->lastSeqNum[rHeader.subtype], rHeader.seqnum);
	  gret->seqNumSkipped[rHeader.subtype] += (deltaSeqNum-1);
	  gret->timesSeqNumSkipped[rHeader.subtype]++;
	}
	gret->lastSeqNum[rHeader.subtype] = rHeader.seqnum;
      }
      printf("SN: %d\n", gret->lastSeqNum[rHeader.subtype]);
    }
    break;
  case 4:
    {
      gret->getMode3(inf, rHeader.length, rHeader.subtype, rHeader.type);
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhTS = rHeader.timestamp;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhLength = rHeader.length;
      gret->eventCnt[rHeader.subtype]++;
      if (gret->lastSeqNum[rHeader.subtype] == -1) {
	gret->lastSeqNum[rHeader.subtype] = rHeader.seqnum;
	if (1) 	printf("\n Crystal %d - first sequence number in file %d\n", rHeader.subtype, rHeader.seqnum);
      } else {
	Int_t deltaSeqNum = rHeader.seqnum - gret->lastSeqNum[rHeader.subtype];
	if (gret->lastSeqNum[rHeader.subtype] == 65535) {
	  deltaSeqNum = rHeader.seqnum + 65536 - gret->lastSeqNum[rHeader.subtype];
	}
	if (deltaSeqNum < 1) {
	  if (1) {
	    printf("\nCrystal %d - sequence numbers out of order.  Previous %d, now %d\n", rHeader.subtype,
		   gret->lastSeqNum[rHeader.subtype], rHeader.seqnum);
	  }
	  gret->seqNumOOO[rHeader.subtype]++;
	} else if (deltaSeqNum > 1) {
	  if (1) {
	    printf("Crystal %d - sequence numbers skipped.  Previous %d, now %d\n", rHeader.subtype,
		   gret->lastSeqNum[rHeader.subtype], rHeader.seqnum);
	  }
	  gret->seqNumSkipped[rHeader.subtype] += (deltaSeqNum-1);
	  gret->timesSeqNumSkipped[rHeader.subtype]++;
	}
	gret->lastSeqNum[rHeader.subtype] = rHeader.seqnum;
      }
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
      LookForGoodData(inf, junk,lastLength, lastSeqNum);
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

void LookForGoodData(FILE *inf, UShort_t junk[], Int_t lastEvt, Int_t lastSeqNum) {
  printf("Last Event Length = %d, last sequence number = %d\n", lastEvt, lastSeqNum);

  // This happens when the last event was > 9000 (frame size).
  // It gets truncated, so we know how far back to go
  Int_t startSearch = 9000 - lastEvt;
  printf("--> starting search for good event at %d + routing header length from here.\n", startSearch);
  
  fseek(inf, -1*(sizeof(struct routingHdr)), SEEK_CUR);
  fseek(inf, startSearch, SEEK_CUR);
  Int_t success = 0;
  
  for (Int_t i=startSearch; i<lastEvt*4; i++) {
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
	printf("Recovered after bad header with data at %d (sequence number = %d)\n", i, rHeader.seqnum);
	GetData(inf, lastEvt, lastSeqNum);
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
  printf("                        -calibrationRun (fills calibration histograms, no Tree)\n");
  printf("                        -hole <X=1 to 120>  (for analyzing the files from a single Quad only if data is taken crystal-wise)\n");
  printf("\n");
  printf("  Output tree gdata has g2Out (full g2OUT) + coincTS per coincidence.\n");
  printf("  Run tracking as pass 2:  trackGreta -i <rootFile> -o tracked.root -trackingChat <chat>\n");
  printf("  Or parallel:  scripts/trackGretaParallel.sh -j N -i ... -o ... -trackingChat ...\n");
  printf("\n");
  printf("  Note - if you give a path to a run with many subfiles, just give the first part of the file name (e.g. the part\n");
  printf("         that tab-completes and it will scan ALL files.  If you specify a hole number it will only scan the \n");
  printf("         corresponding four files (if present).  If after '-f' you give a complete filename it will only scan that.\n");
}

