#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <signal.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "streamFunctions.h"
#include "GRETINA.h"

int gotsignal;
void breakhandler(int dummy) {
  printf("Got break signal.  Aborting cleanly...\n");
  gotsignal = 1;
}

int main(int argc, char **argv) {

  gotsignal = 0;
  signal(SIGINT, breakhandler);

  FILE *fout;
  int num, curr, newRead;

  int nReads = 10000;
  
  unsigned long long int ledCrossings;
  int ledCrossing = 0;
  
  int pileUp = 0;
  int ledCrossed = 0;
  long long int startTS = 0;
  long long int ledTS = 0;
  long long int currTS = 0;

  int withBLR = 1;
  int BLRcountdown = 0;

  if (argc < 7) {
    fprintf(stderr, "Minimum usage: Analyze -fIn <input filename> -fSet <settingsName> -n <nReads> \n");
    exit(1);
  }
  
  TString rootoutputName, inputName, fileNameSet;
  int inputFile = 0, outputFile = 0, setFile = 0;

  int i = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-fIn") == 0) {
      inputName = argv[i+1]; i++; i++;
      inputFile = 1;
    } else if (strcmp(argv[i], "-fOut") == 0) {
      rootoutputName = argv[i+1]; i++; i++;
      outputFile = 1;
    } else if (strcmp(argv[i], "-fSet") == 0) {
      fileNameSet = argv[i+1]; i++; i++;
      setFile = 1;
    } else if (strcmp(argv[i], "-n") == 0) {
      sscanf(argv[i+1], "%d", &nReads);
      if (nReads == -1) { nReads = 1000000; }
      i++; i++;
    } else {
      cout << ALERTTEXT;
      printf("Error -- unrecognized input flag: <%s>\n", argv[i]);
      cout << RESET_COLOR; fflush(stdout);
      exit(-1);
    }
  }

  if (!inputFile) {
    cout << ALERTTEXT;
    printf("Error -- missing arguments!  Try again!\n");
    cout << RESET_COLOR; fflush(stdout);
    exit(-1);
  }
  if (!setFile) {
    fileNameSet = "settings.set";
  }
  
  streamer *data = new streamer();
  curr = data->Initialize(inputName, fileNameSet);
  data->reportSettings();

  int overlapWidth = 2*(2*data->EM + data->EK);
      
  startTS = 0;
  double pzSum = 0;
  int indexStart = 0;

  for (int numberOfReads = 1; numberOfReads<=nReads; numberOfReads++) {
    
    if (numberOfReads == 1) { indexStart = 0; curr -= overlapWidth/2; startTS = 0;}
    else { indexStart = overlapWidth/2; curr += overlapWidth/2; }

    /* First basic filtering ... 
          for first read, we go from [0] to [END - 1/2 of overlap]
	  for subsequent reads, we go from [1/2 of overlap] to [END - 1/2 of overlap] */
    data->doLEDfilter(indexStart, curr, startTS);
    ledCrossings = data->getLEDcrossings(indexStart, curr, startTS);
    ledCrossing += ledCrossings;
    printf("numberOfReads = %d, ledCrossing = %d\n", numberOfReads, ledCrossing);  
    startTS += (curr - overlapWidth/2);
    

    if (numberOfReads < nReads) {
      curr = data->Reset(overlapWidth, curr+overlapWidth/2); /* Clears and gets the next data bunch... */
      if (curr <= 0) { break; }
    }

    if (gotsignal == 1) { numberOfReads = nReads+1; printf("Exiting gracefully.\n");}
  }

  printf("LED crossings observed = %d\n", ledCrossing);
  printf("MB read = %d\n", data->bytesRead/(1024*1024));

}

