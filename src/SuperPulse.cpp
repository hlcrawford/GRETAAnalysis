#include "GRETA.h"

ClassImp(SuperPulse);

void SuperPulse::Initialize(controlVariables* ctrl) {
  printf("Initializing superpulse section... \n");
  CFD_INT_LEN = 4;
  CFD_DELAY = 4;
  CFD_FRACTION = 4;
  TR_SCALE = 10000.;
  AVG_TR_STRIDE = 110;
  AVG_TR_LENGTH = 108;
  
  lowE = 0; highE = 0;
  
  for (Int_t i=0; i<121; i++) {
    mult[i] = 0;
    ccE[i] = 0.;
    segE[i] = 0.;
    netSeg[i] = 0;

    for (Int_t j=0; j<40; j++) {
      segs[i][j] = 0;
      if (j<36) { data4net[i][j] = 0; }
      for (Int_t k=0; k<4096; k++) {
	averageTrace[i][j][k] = 0.0;
	averageTraceINT[i][j][k] = 0;
      }
    }
  }

  /* Check if this is a superpulse analysis, if it is, 
     get calibrations, etc. */
  if (ctrl->superPulse) {
    lowE = ctrl->spLowE;
    highE = ctrl->spHighE;
    ReadParams(ctrl->spXtalkFile, "delay1", delay1, 40);

    for (Int_t i=0; i<121; i++) {
      delay1[i][39] = 0.0f;
    }
  }
  
}

/****************************************************/

Int_t SuperPulse::ReadParams(TString filename, const char *label, Float_t x[][40], Int_t len) {

  FILE *fin;
  Int_t max_length = 120;
  char s[max_length], *tok0, *tok, *u, *s0;
  Int_t nn = 0;
  Float_t a;

  fin = fopen(filename.Data(), "r");
  if (fin == NULL) {
    printf("Could not open \"%s\".\n", filename.Data() );
  } else {
    printf("Opened \"%s\" to get the delay parameters.\n", filename.Data());
    while (fgets(s, max_length, fin) != 0) {
      tok0 = strtok(s, " ");
      if ((u = strchr(tok0, atoi(":"))) != 0) {
	*u = '\0';
	if (strcmp(label, tok0) == 0) {
	  while(1) {
	    s0 = fgets(s, max_length, fin);
	    if (s0 == 0 || s[0] == '\n') {
	      fclose(fin);
	      return nn;
	    }
	    s0 = s;
	    while ((tok = strtok(s0, " \t\n")) != 0) {
	      s0 = 0;
	      a = (Float_t) atof(tok);
	      for (Int_t xNum=0; xNum<121; xNum++) {
		x[xNum][nn++] = a;
	      }
	      if (nn >= len) {
		fclose(fin);
		return nn;
	      }
	    }
	  }
	}
      }
    }
    fclose(fin);
  }
  printf("Finished reading param file.\n");
  return 0;
}

/****************************************************/
 
void SuperPulse::MakeSuperPulses() {
  for (Int_t i=0; i<121; i++) {
    if (mult[i] == 1) { /* Segment multiplicity 1 */
      if (segE[i] >= lowE && ccE[i] >= lowE && segE[i] <= highE && ccE[i] <= highE) {
	Int_t t0 = AlignCFD(i);
	if (t0 >= 0) { /* CFD alignment was successful. */
	  data4net[i][netSeg[i]]++;
	  for (Int_t m=0; m<37; m++) { /* 36 segments + 1 CC only */
	    for (Int_t j=0; j<AVG_TR_LENGTH; j++) {
	      /* We fill one giant array with all 37 traces, with small gaps between waveforms. */
	      averageTrace[i][netSeg[i]][m*(AVG_TR_STRIDE) + j] += waves[i][m][j];
	    } /* Loop over waveform samples */
	  } /* Loop over segments */
	} /* if CFD alignment is OK */
      } /* in energy window */
    } /* if mult = 1 */
  } /* loop over MAXCRYSTALS crystals */
}

/****************************************************/

Int_t SuperPulse::AlignCFD(Int_t crystalNum) {
  Int_t i, j, k;
  Float_t t=0, tmin = 1000.0;

  /* Find the time to align to. */
  for (j=0; j<40; j++) {
    if (netSeg[crystalNum] == j || (j==36)) { /* Net segment or CC */
      t = cfdTime(crystalNum, j);
      if (t < 0.0) { return -1; }
      t -= delay1[crystalNum][j]; /* DCR - there may be a question about - sign */
      if (tmin > t) { tmin = t; }
    }
  }
  i = (Int_t)(t + 0.5) - 16; /* Starting step of the trace to extract for SP */
  if (i == 0) { return 0; }
  
  /* Shift the traces.  Order depends on if this is a left or right shift. */
  if (i < 0) {
    for (j=0; j<40; j++) {
      for (k=trLen - 2; k > (0-i); k--) {
	waves[crystalNum][j][k] = waves[crystalNum][j][k+i];
      }
      for (k=0; k<(0-i); k++) {
	waves[crystalNum][j][k] = 0;
      }
    }
  } else {
    for (j=0; j<40; j++) {
      for (k=0; k<=trLen - i - 2; k++) {
	waves[crystalNum][j][k] = waves[crystalNum][j][k+i];
      }
      for (k=trLen - i - 1; k<trLen; k++) {
	waves[crystalNum][j][k] = waves[crystalNum][j][trLen-2];
      }
    }
  }

  return 0;
}

/****************************************************/

Float_t SuperPulse::cfdTime(Int_t crystalNum, Int_t segNum) {
  /* Ge-style CFD to get time from the signal trace.  
     From DCR. */
  Int_t deriv[2048];
  Int_t cfd[2048];
  Int_t i, imax = 0, max_deriv = 0;
  
  deriv[0] = 0;
  for (i=0; i< CFD_INT_LEN; i++) {
    deriv[0] += (waves[crystalNum][segNum][i+CFD_INT_LEN] -
		 waves[crystalNum][segNum][i]);
  }
  for (i=1; i<trLen - 5 - 2*CFD_INT_LEN; i++) {
    deriv[i] = (deriv[i-1] + 
		waves[crystalNum][segNum][i+2*CFD_INT_LEN] -
		2*waves[crystalNum][segNum][i+CFD_INT_LEN] + 
		waves[crystalNum][segNum][i]);
    if (max_deriv < deriv[i]) {
      max_deriv = deriv[i];
      imax = i; 
    }
  }
  for (i=0; i<trLen - 5 - 2*CFD_INT_LEN - CFD_DELAY; i++) {
    cfd[i] = deriv[i] - deriv[i+CFD_DELAY]/CFD_FRACTION;
  }

  for (i=imax + CFD_DELAY; i>0; i--) {
    if (cfd[i] <= 0 && cfd[i+1] > 0) {
      /* Interpolate zero crossing and return time in steps. */
      return ((Float_t)i) - ((Float_t)cfd[i]) / ((Float_t)(cfd[i+1] - cfd[i]));
      break;
    }
  }
  return -1.0;
}

/****************************************************/

void SuperPulse::FinishSuperPulses() {
  /* Scale the superpulses, and get the trace gains. */
  
  for (Int_t i=0; i<121; i++) {
    
    for (Int_t j=0; j<36; j++) {
      gain[i][j] = (((Float_t)averageTrace[i][j][j*AVG_TR_STRIDE + AVG_TR_LENGTH - 1]) /
		    ((Float_t)averageTrace[i][j][36*AVG_TR_STRIDE + AVG_TR_LENGTH - 1]));
    }

    gain[i][36] = 1.0f;
    
    for (Int_t j=0; j<36; j++) {
      Float_t scaleFactor = TR_SCALE / ((Float_t)averageTrace[i][j][36*AVG_TR_STRIDE + AVG_TR_LENGTH -1]);
      for (Int_t k=0; k<37; k++) {
	if (k<36) {
	  for (Int_t m=0; m<AVG_TR_LENGTH; m++) {
	    averageTrace[i][j][k*AVG_TR_STRIDE + m] = ((Float_t)averageTrace[i][j][k*AVG_TR_STRIDE + m])*scaleFactor/gain[i][k];
	    averageTraceINT[i][j][k*AVG_TR_STRIDE + m] = (Int_t)averageTrace[i][j][k*AVG_TR_STRIDE + m];
	  }
	} else {
	  for (Int_t m=0; m<AVG_TR_LENGTH; m++) {
	    averageTrace[i][j][k*AVG_TR_STRIDE + m] = ((Float_t)averageTrace[i][j][k*AVG_TR_STRIDE + m])*scaleFactor/gain[i][k];
	    averageTraceINT[i][j][k*AVG_TR_STRIDE + m] = (Int_t)averageTrace[i][j][k*AVG_TR_STRIDE + m];
	  }
	  
	}
      }
    }
  }
}

/****************************************************/

void SuperPulse::WriteSuperPulses() {
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    if (data4net[i][0] > 0 || data4net[i][1] > 0) {
      char filenameOut[1024];
      sprintf(filenameOut, "SPCrystal%d_%d.spn", i, (Int_t)((lowE + highE)/2));
      cout << "Writing superpulses to " << filenameOut << endl;
      
      FILE *spOut = fopen(filenameOut, "wb");
      for (Int_t j=0; j<36; j++) {
	fwrite(averageTraceINT[i][j], sizeof(Int_t), 4096, spOut);
      }
      fclose(spOut);
      printf("Superpulse statistics:\n");
      for (Int_t k=0; k<36; k++) {
	printf("%d -- %d\n", k, data4net[i][k]);
      }
    }
  }
}

    
							    

