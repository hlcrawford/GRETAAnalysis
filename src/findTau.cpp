#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH2.h"
#include "TF1.h"

#include "streamFunctions.h"

int gotsignal;
void breakhandler(int dummy) {
  printf("Got break signal.  Aborting cleanly...\n");
  gotsignal = 1;
}

Double_t doublePole(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return par[0]*(par[2])*TMath::Exp(-xx/par[1]) + par[0]*(1-par[2])*TMath::Exp(-xx/par[3]) + par[4];
}

Double_t singlePole(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return par[0]*TMath::Exp(-xx/par[1]) + par[2];
}

int main(int argc, char **argv) {

  FILE *fout;
  int num, curr, i, newRead;

  int nReads = 10000;
  
  unsigned long long int ledCrossings;
  int ledCrossing = 0;
  
  long long int startTS = 0;
  long long int ledTS = 0;
  long long int currTS = 0;

  int twoPoles = 0;

  Double_t xVal[110000], yVal[110000];
  for (i=0; i<110000; i++) { xVal[i] = i; }

  if (argc < 9) {
    fprintf(stderr, "Minimum usage: getTau -fIn <input filename> -fOut <rootOutputName> -fSet <settingsName> -n <nReads> \n");
    exit(1);
  }
  
  TString rootoutputName, inputName, fileNameSet;
  int inputFile = 0, outputFile = 0, setFile = 0;

  i = 1;
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
      sscanf(argv[i+1], "%d", &nReads); i++; i++;
    } else {
      std::cout << ALERTTEXT;
      printf("Error -- unrecognized input flag: <%s>\n", argv[i]);
      std::cout << RESET_COLOR; fflush(stdout);
      exit(-1);
    }
  }

  if (!inputFile || !outputFile) {
    std::cout << ALERTTEXT;
    printf("Error -- missing arguments!  Try again!\n");
    std::cout << RESET_COLOR; fflush(stdout);
    exit(-1);
  }
  if (!setFile) {
    fileNameSet = "settings.set";
  }
  
  streamer *data = new streamer();
  curr = data->Initialize(inputName, fileNameSet);
  
  int overlapWidth = 2*(2*data->EM + data->EK);

  TFile *rootOUT = new TFile(rootoutputName.Data(), "RECREATE");

  TGraph *gr;
  TH2F *tauValues;
  if (!twoPoles) {
    tauValues = new TH2F("tau", "tau", 3000, 0, 30000, 3000, 0, 3000);
  } else {
    tauValues = new TH2F("tau", "tau", 3000, 0, 30000, 3000, 0, 30000);
  }

  
  TF1 *expofit;
  if (!twoPoles) {
    expofit = new TF1("expofit", singlePole, 0, 100000, 3);
    expofit->SetParameters(100, 5000, 0);
    expofit->SetParLimits(1, 2000, 8000);
  } else {
    expofit = new TF1("expofit", doublePole, 0, 100000, 5);
    expofit->SetParameters(100, 5000, 0.9, 20000, 0);
    expofit->SetParLimits(0, 0, 100000);
    expofit->SetParLimits(1, 2000, 8000);
    expofit->SetParLimits(2, 0, 1);
    expofit->SetParLimits(3, 20, 3000000);
  }
	
  startTS = 0;
  double pzSum = 0;
  int indexStart = 0;

  for (int numberOfReads = 1; numberOfReads<=nReads; numberOfReads++) {

    if (numberOfReads == 1) { indexStart = 0; }
    else { indexStart = overlapWidth; }

    /* Actual analysis is from these functions for a given chunk of waveform */
    data->doLEDfilter(indexStart, curr, startTS);
    ledCrossings = data->getLEDcrossings(indexStart, curr, startTS);
    ledCrossing += ledCrossings;

    printf("numberOfReads = %d, ledCrossing = %d\r", numberOfReads, ledCrossing);

    for (i=0; i<data->ledOUT.size(); i++) {
      if ( (data->ledOUT[i+1] - data->ledOUT[i]) > 60000) {
	for (int j=0; j<50000; j++) { yVal[j] = data->wf[data->ledOUT[i] - startTS + 200 + j]; }
	expofit->SetParameter(0, yVal[1]);
	expofit->SetParameter(1, 5000);
	expofit->SetParameter(2, yVal[49999]);
	if (twoPoles) {
	  expofit->SetParameter(2, 0.9);
	  expofit->SetParameter(3, 20000);
	  expofit->SetParameter(4, yVal[49999]);
	}
	gr = new TGraph(50000, xVal, yVal);
	gr->Fit("expofit", "WWQ");
	if (!twoPoles && expofit->GetChisquare()/expofit->GetNDF() > 30000) {
	  expofit->SetParameters(100, 5000);
	  gr->Fit("expofit", "WWQ");
	}
	if (twoPoles && expofit->GetChisquare()/expofit->GetNDF() > 3000) {
	  expofit->SetParameters(100, 5000, 0.9, 20000, 0);
	  gr->Fit("expofit", "WWQ");
	}
	gr->SetName(Form("gr%d", i));
	gr->Write();
	if (!twoPoles) {
	  tauValues->Fill(expofit->GetParameter(1), expofit->GetChisquare()/expofit->GetNDF());
	  printf("%d:\t Tau:%0.3f,\tReducedChi2:%0.1f\n", i, expofit->GetParameter(1), expofit->GetChisquare()/expofit->GetNDF());
	} else {
	  if (expofit->GetParameter(1) < expofit->GetParameter(3)) {
	    tauValues->Fill(expofit->GetParameter(1), expofit->GetParameter(3));
	  } else {
	    tauValues->Fill(expofit->GetParameter(3), expofit->GetParameter(1));
	  }
	  printf("%d:\tTauA: %0.3f,\tTauB:%0.3f,\tFractionA: %0.3f,\tReducedChi2:%0.1f \n", i, expofit->GetParameter(1), expofit->GetParameter(3), expofit->GetParameter(2), expofit->GetChisquare()/expofit->GetNDF());
	}
      }
    }

    data->ledOUT.clear();
    
    startTS += (curr - overlapWidth);

    if (numberOfReads < nReads) {
      curr = data->Reset(overlapWidth, curr+overlapWidth/2); /* Clears and gets the next data bunch... */
      if (curr <= 0) { break; }
    }

  }

  printf("LED crossings observed = %d\n", ledCrossing);
  
  tauValues->Write();
  rootOUT->Write();
  rootOUT->Close();
}

