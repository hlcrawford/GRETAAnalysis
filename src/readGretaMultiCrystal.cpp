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

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "Globals.h"

//#include "GRETA.h"

using namespace std;

#define CFD_INT_LEN 4
#define CFD_DELAY 4
#define CFD_FRACTION 4
#define AVG_TR_LENGTH 108
#define AVG_TR_STRIDE 110
#define TR_SCALE 10000

#define EB_DIFF_TIME 100

#define NUM_CHAN 40

#define DEBUG 0

void GetData(FILE* inf);
void SkipData(FILE *inf, UShort_t junk[]);

/* Utility functions for interrupting the sort cleanly, and for terminal progress bar. */

Int_t gotSignal;
void breakhandler(int dummy) {
  cout << "Got break signal.  Aborting cleanly..." << endl;
  gotSignal = 1;
}

void progressB(int pct) {
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

  if (argc < 3) {
    printf("Usage: readGreta <InputFileWithPath> <ROOTOutputName>");
  }

  gret = new GRETA();
  gret->Initialize();

  //  gretaWaveformMsg wform;

  float ccEnergy[4];
  float segEnergy[36];
  float bl[40];
  vector<float> pz;
  float sumA = 0.; float sumB = 0.;
  int16_t trReOrder[NUM_CHAN][512]; 

  for (Int_t i=0; i<120; i++) {
    for (Int_t j=0; j<40; j++) {
      gret->gain[i][j] = 1.0;
      gret->offset[i][j] = 0.0;
    }
  }

  gret->gain[0][36] = 0.01879136; gret->offset[0][36] = -0.0385907;
  gret->gain[1][36] = 0.02006449; gret->offset[1][36] =  0.02534448;
  gret->gain[2][36] = 0.01918792; gret->offset[2][36] = -0.9469749;
  gret->gain[3][36] = 0.01979890; gret->offset[3][36] = -0.2255907;


  int mult = 0;
  float cfdTime = -1.;
  int data4net[36] = {0};
  float averageTrace[40][4096];
  int averageTraceINT[40][4096];
  float traceGain[40] = {0.0};

  gret->segMap.insert(std::pair<int, int>(29, 0));  
  gret->segMap.insert(std::pair<int, int>(23, 1));
  gret->segMap.insert(std::pair<int, int>(35, 2));   
  gret->segMap.insert(std::pair<int, int>(11, 3));
  gret->segMap.insert(std::pair<int, int>(5, 4));    
  gret->segMap.insert(std::pair<int, int>(24, 5));
  gret->segMap.insert(std::pair<int, int>(16, 6));   
  gret->segMap.insert(std::pair<int, int>(9, 7));
  gret->segMap.insert(std::pair<int, int>(3, 8));    
  gret->segMap.insert(std::pair<int, int>(33, 9));
  gret->segMap.insert(std::pair<int, int>(15, 10));  
  gret->segMap.insert(std::pair<int, int>(27, 11));
  gret->segMap.insert(std::pair<int, int>(17, 12));  
  gret->segMap.insert(std::pair<int, int>(10, 13));
  gret->segMap.insert(std::pair<int, int>(4, 14));   
  gret->segMap.insert(std::pair<int, int>(34, 15));
  gret->segMap.insert(std::pair<int, int>(28, 16));  
  gret->segMap.insert(std::pair<int, int>(22, 17));
  gret->segMap.insert(std::pair<int, int>(25, 18));  
  gret->segMap.insert(std::pair<int, int>(19, 19));
  gret->segMap.insert(std::pair<int, int>(20, 20)); 
  gret->segMap.insert(std::pair<int, int>(32, 21));
  gret->segMap.insert(std::pair<int, int>(2, 22));   
  gret->segMap.insert(std::pair<int, int>(21, 23));
  gret->segMap.insert(std::pair<int, int>(12, 24));  
  gret->segMap.insert(std::pair<int, int>(18, 25));
  gret->segMap.insert(std::pair<int, int>(6, 26));  
  gret->segMap.insert(std::pair<int, int>(0, 27));
  gret->segMap.insert(std::pair<int, int>(30, 28));  
  gret->segMap.insert(std::pair<int, int>(31, 29));
  gret->segMap.insert(std::pair<int, int>(13, 30));  
  gret->segMap.insert(std::pair<int, int>(7, 31));
  gret->segMap.insert(std::pair<int, int>(1, 32));   
  gret->segMap.insert(std::pair<int, int>(26, 33));
  gret->segMap.insert(std::pair<int, int>(14, 34));  
  gret->segMap.insert(std::pair<int, int>(8, 35));
  gret->segMap.insert(std::pair<int, int>(36, 36));  
  gret->segMap.insert(std::pair<int, int>(37, 37));
  gret->segMap.insert(std::pair<int, int>(38, 38));  
  gret->segMap.insert(std::pair<int, int>(39, 39));
  
  /* Initialize ROOT output */
  TFile *fOut = new TFile(argv[2], "RECREATE");
  TTree *data = new TTree("data", "data");
  data->Branch("g3", "g3OUT", &(gret->g3Out));
  
  unsigned char buf[65536];
  UShort_t junk[8192];

  /* Check file length for progress */
  struct stat fileStatus;
  stat(argv[1], &fileStatus);
  int64_t bytesInFile = (int64_t) fileStatus.st_size;
  cout << "Input data file size is " << (Float_t)bytesInFile/1024./1024./1024. << " GB\n\n";
  
  int64_t bytesRead = 0;

  FILE *inf;
  inf = fopen(argv[1], "r");

  Int_t siz = 0;
  siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
  bytesRead += sizeof(struct routingHdr);

  long long int lastTS = 0;  long long int currTS = 0;
  long long int deltaEvent = 0;

  Int_t TSreports = 0;

  while (siz && !gotSignal) {
    rHeader.seqnum = ntohs(rHeader.seqnum);
    rHeader.timestamp = (int64_t)ntoh64((uint64_t)rHeader.timestamp);
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
      printf("    Checksum: %i\n", rHeader.checksum);
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
      
      GetData(inf);
      //printf("Returned from GetData\n");
      bytesRead += rHeader.length;

    } else {

      data->Fill();
      gret->Reset();
      currTS = rHeader.timestamp;
      deltaEvent = (Float_t)(rHeader.timestamp - currTS);
      
      GetData(inf);
      bytesRead += rHeader.length;
    }

    progressB(100*bytesRead/bytesInFile);    
    siz = fread(&rHeader, sizeof(struct routingHdr), 1, inf);
    bytesRead += sizeof(struct routingHdr);
    
  }

  data->Write();
  fOut->Write();
  fOut->Close();
  printf("\nDone and done.\n\n");
  return 0;
}

void GetData(FILE* inf) {
  UShort_t junk[8192];

  // printf("In GetData - %d\n", rHeader.type);

  switch (rHeader.type) {
  case 3:
    {
      gret->getMode3(inf, rHeader.length, rHeader.subtype);
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSubType = rHeader.subtype;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhSequence = rHeader.seqnum;
      gret->g3Out.xtals[gret->g3Out.crystalMult()-1].rhTS = rHeader.timestamp;
    }
    break;
  default:
    {
      cout << "Routing Header type not recognized: " << rHeader.type << endl;
      SkipData(inf, junk);
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




/*
    if (rHeader.type == 3) {
      siz = fread(&wform, 1, ntohs(rHeader.length), inf);
      bytesRead += ntohs(rHeader.length);
      // tmp = (buf);
      //memmove(&wform, &buf, sizeof(struct gretaWaveformMsg));
      // wform = (gretaWaveformMsg*)buf;

      wform.tr_len = ntohs(wform.tr_len);
      wform.trig_src = ntohs(wform.trig_src);
      wform.pad = ntohs(wform.pad);
      mult = 0;
      for (int i=0; i<NUM_CHAN; i++) {
	wform.ener[i] = ntohl(wform.ener[i]);
	for (int j=0; j<512; j++) {
	  wform.tr[i][j] = ntohs(wform.tr[i][j]);	 
	}
      }

      printf("Printing trace for ch 36\n");
      for (int i=0; i<188; i++) {
	cout << i << " " << (&wform.tr[0][0])[36*188+i] << endl;
      }
      cin.get();

      /* Remap the traces now... */
      // if (wform.tr_len < 512) {
      // 	for (int ch = NUM_CHAN; ch>=0; ch--) { 
      // 	  float blsum = 0.0;
      // 	  for (int i=0; i<wform.tr_len; i++) {
      // 	    wform.tr[ch][i] = (&wform.tr[0][0])[ch*wform.tr_len+i];
      // 	    if (i<20) { 
      // 	      blsum += (float)(wform.tr[ch][i]);
      // 	    }
      // 	  }
      // 	  for (int i=wform.tr_len; i<512; i++) {
      // 	    wform.tr[ch][i] = 0.0;
      // 	  }
      // 	  bl[ch] = blsum/20.;	
      // 	}
      // }  */
  
/*
    
      if (DEBUG) {
	printf("Waveform message!\n");
	printf("  Version:  %i\n", wform.version);
	printf("  ID:           %i\n", wform.id);
	printf("  tr_len:      %d\n", wform.tr_len);
	printf("  trig_src:   %i\n", wform.trig_src);
	printf("  pad:         %i\n", wform.pad);
	printf("  timestamp: %lld\n", (int64_t)ntoh64((uint64_t)wform.timestamp));
	for (int i=0; i<40; i++) {
	  printf("  Raw energy %d: %d\n", i, wform.ener[i]);
	}
      	cin.get();
      }
	
      for (int i=0; i<NUM_CHAN; i++) {
	if (0) {	
	  //float ipz = 0;
	  //float tau =1./ 5100.;
	  float base = 0.0;
	  //sumA = 0.; sumB = 0.;
	  for (int j=0; j<30; j++) {
	    base += wform.tr[i][10+j];
	  }
	  base /= 30.;
	  //pz.clear();
	  //for (int j=50; j<wform.tr_len; j++) {
	  // pz.push_back(wform.tr[i][j] - base + ipz*tau);
	  // ipz += wform.tr[i][j] - base;
	  //}
	  //for (int j=10; j<100; j++) {
	  // sumA += pz[j];
	  //sumA += pz[j+90];//for 356-156 trace window
	  //}
	  //for (int j=0; j<100; j++) {
	  // sumB += pz[pz.size()-150-j];
	  //sumB += pz[pz.size()-10-j]; // for 356-156 trace window
	  //}
	  //if (i >= 36) {
	  // ccEnergy[i-36] = (sumB - sumA)*gain[i];
	  //} else {
	  // segEnergy[i] = (sumB - sumA)*gain[i];
	  // if (segEnergy[i] > 30.) { mult++; }
	  //}
	  
	  if (i<36) {
	    segEnergy[i] = (float)wform.ener[i]*gain[i];
	    if (segEnergy[i] > 30) { mult++; }
	  } else {
	    ccEnergy[i-36] = (float)wform.ener[i]*gain[i];
	  }
	  
	  
	  for (int j=0; j<wform.tr_len; j++) {
	    wform.tr[i][j] -= (int16_t)base;
	  }  
	}
      }

      if (0) {
	for (int i=0; i<NUM_CHAN; i++) {
	  for (int j=0; j<wform.tr_len; j++) {
	    trReOrder[i][j] = wform.tr[segMap.find(i)->second][j];
	  }
	}
	
	if (mult==1) {
	  cfdTime = 0.0;
	  for (int ch=0; ch<36; ch++) {
	    if (segEnergy[ch] >= 1330 && segEnergy[ch] <= 1334 &&
		ccEnergy[0] >= 1330 && ccEnergy[0] <= 1334) {
	      /* Time alignment! */
/*	      Int_t deriv[1024];
	      Int_t cfd[1024];
	      Int_t imax = 0, max_deriv = 0;
	      deriv[0] = 0;
	      for (int j=0; j<CFD_INT_LEN; j++) {
		deriv[0] += (wform.tr[36][j+CFD_INT_LEN] -  wform.tr[36][j]);
	      }
	      for (int j=1; j<wform.tr_len - 5 - 2*CFD_INT_LEN; j++) {
		deriv[j] = (deriv[j-1] + 
			    wform.tr[36][j+2*CFD_INT_LEN] -
			    2*wform.tr[36][j+CFD_INT_LEN] + 
			    wform.tr[36][j]);
		if (max_deriv < deriv[j]) {
		  max_deriv = deriv[j];
		  imax = j;
		}
	      }
	      for (int j=0; j<wform.tr_len - 5 - 2*CFD_INT_LEN - CFD_DELAY; j++) {
		cfd[j] = deriv[j] - deriv[j+CFD_DELAY]/CFD_FRACTION;
	      }
	      for (int j=imax + CFD_DELAY; j>0; j--) {
		if (cfd[j] <=0 && cfd[j+1] > 0) {
		  cfdTime = ((float)j) - ((float)cfd[j])/((float)(cfd[j+1]-cfd[j]));
		  break; 
		}
	      }
	      
	      if (cfdTime > 0) { 
		/* This data set is centered with tCFD = 65 */
/*		int delay = (int)(cfdTime + 0.5) - 16;
		delay -= (65-16);
		
		/* Shift traces */
/*		if (delay < 0) {
		  for (int j=0; j<40; j++) {
		    for (int k=wform.tr_len-2; k > (0-delay); k--) {
		      wform.tr[j][k] = wform.tr[j][k+delay];
		    }
		    for (int k=0; k < (0-delay); k++) {
		      wform.tr[j][k] = 0;
		    }
		  }
		} else {
		  for (int j=0; j<40; j++) {
		    for (int k=0; k<=wform.tr_len-delay-2; k++) {
		      wform.tr[j][k] = wform.tr[j][k+delay];
		    }
		    for (int k=wform.tr_len-delay-1; k < wform.tr_len; k++) {
		      wform.tr[j][k] = wform.tr[j][wform.tr_len-2];
		    }
		  }
		}
		
		data4net[ch]++;
		int whichToFill = -1;
		for (int j=0; j<36; j++) { if (segMap.find(j)->second == ch) { whichToFill = segMap.find(j)->first; } }
		for (int j=0; j<37; j++) {
		  for (int m=0; m<AVG_TR_LENGTH; m++) {
		    
		    averageTrace[whichToFill][j*(AVG_TR_STRIDE)+m] += (float)wform.tr[segMap.find(j)->second][m+65];
		    //printf("For seg %d, using trace %d\n", j, segMap.find(j)->second);
		  }
		}
	      }
	    }
	  }
	}
      }
      

      
      /* Fill tree... */
/*    data->Fill();
      
    } else {

    }
    

    // rHeader = (routingHdr*)buf;
  }

  data->Write();
  fOut->Write();
  fOut->Close();
  printf("\nDone and done.\n\n");
  
  return 0;
}
*/
