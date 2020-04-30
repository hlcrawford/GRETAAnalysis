#include "streamFunctions.h"

int streamer::Initialize(TString inputFileName, TString setFileName) {
   
  tau = 0.; DV = 0; LEDThreshold = 0; BLRinhibitLength = 0; BLRretrigger = 0; 
  EK = 0;  EM = 0;  LEDK = 0; 
  useBLR = 0;  usePO = 0;
  POtime = 0;

  invertWF = 0;

  READSIZE = 1000*1000*10; // 1 second at a time for the original streamer

  getSettings(setFileName);

  if (tau == 0) { tau = 5000.; }
  if (EK == 0) { EK = 50; }
  if (EM == 0) { EM = 400; }
  if (LEDK == 0) { LEDK = 5; }
  if (DV == 0) { DV = 8; }
  if (LEDThreshold == 0) { LEDThreshold = 50; }
  if (BLRinhibitLength == 0) { BLRinhibitLength = EM*2 + EK + 300; }
  if (BLRretrigger == 0) { BLRretrigger = -1; }

  fin = fopen(inputFileName.Data(), "r");
  if (fin == 0) {
    fprintf(stderr, "Cannot open the input file %s\n", inputFileName.Data());
    exit(1);
  }
  
  bytesRead = 0;
  int sizeArrays = 0;	
  
  if (FWB) { 
    fread(&hdr, sizeof(struct routingHdr), 1, fin);
    bytesRead += sizeof(struct routingHdr);
    READSIZE = ntohs(hdr.length);
    
    printf("First Routing Header: \n");
    printf("    Version:  0x%x\n", hdr.version);
    printf("    Flags:    %i\n", hdr.flags);
    printf("    Type:     %i\n", hdr.type);
    printf("    SubType:  %i\n", hdr.subtype);
    printf("    Length:   %i\n", ntohs(hdr.length));
    printf("    SeqNum:   %i\n", ntohs(hdr.seqnum));
    
    READSIZE /= sizeof(short int);
    sizeArrays = 101*READSIZE;
    
  } else {
    sizeArrays = 2*READSIZE;
  }
  
  wf = (short int*)calloc(sizeArrays, sizeof(short int));
  ledBuf = (int*)calloc(sizeArrays, sizeof(int));
  trapBuf = (int*)calloc(sizeArrays, sizeof(int));
  pzBuf = (double*)calloc(sizeArrays, sizeof(double));
  subBL1 = (double*)calloc(sizeArrays, sizeof(double));
  subBL2 = (double*)calloc(sizeArrays, sizeof(double));
  BL1 = (double*)calloc(sizeArrays, sizeof(double));
  BL2 = (double*)calloc(sizeArrays, sizeof(double));
  pzBLBuf = (double*)calloc(sizeArrays, sizeof(double));
  tempBuf = (double*)calloc(sizeArrays, sizeof(double));
  BL1[0] = 0;  BL2[0] = 0; subBL1[0] = 0; subBL2[0] = 0;

  int num = fread(wf, 1, READSIZE*sizeof(short int), fin);
  bytesRead += num;

  if (FWB) { for (int j=0; j<READSIZE; j++) { wf[j] = ntohs(wf[j]); } }
  if (invertWF) { for (int j=0; j<2*READSIZE; j++) { wf[j] = -wf[j]; } }
  
  if (FWB) {
    READSIZE *= 101; // Temporary hack until I clean up code, make sure algorithms run fully
  }

  return (num / sizeof(short int));
}

/**************************************************************/

int streamer::Reset(int overlap, int endIndex) {
  
  int index = 0;

  for (int i=0; i<overlap; i++) {
    index = endIndex-overlap+i;
    if (invertWF) wf[i] = -1*wf[index];
    ledBuf[i] = ledBuf[index];
    trapBuf[i] = trapBuf[index];
    pzBuf[i] = pzBuf[index];
    subBL1[i] = subBL1[index];
    subBL2[i] = subBL2[index];
    BL1[i] = BL1[index];
    BL2[i] = BL2[index];
    pzBLBuf[i] = pzBLBuf[index];
    tempBuf[i] = 0.0;
  }
  
  for (int i=overlap; i<endIndex; i++) {
    wf[i] = 0;
    ledBuf[i] = 0;
    trapBuf[i] = 0;
    pzBuf[i] = 0.0;
    subBL1[i] = 0.0;
    subBL2[i] = 0.0;
    BL1[i] = 0.0;
    BL2[i] = 0.0;
    pzBLBuf[i] = 0.0;
    tempBuf[i] = 0.0;
  }

  ledOUT.clear();

int num = 0, total = 0;
int count = 0;
  
  if (FWB) { 
    while (num >= 0 && count < 100) {
        fread(&hdr, sizeof(struct routingHdr), 1, fin);
    	bytesRead += sizeof(struct routingHdr);
    	READSIZE = ntohs(hdr.length);
    
	if (0) {
	  printf("Routing Header: \n");
	  printf("    Version:  0x%x\n", hdr.version);
	  printf("    Flags:    %i\n", hdr.flags);
	  printf("    Type:     %i\n", hdr.type);
	  printf("    SubType:  %i\n", hdr.subtype);
	  printf("    Length:   %i\n", ntohs(hdr.length));
	  printf("    SeqNum:   %i\n", ntohs(hdr.seqnum));
	}
    
    READSIZE /= sizeof(short int);

    num = fread(wf+overlap+count*READSIZE, 1, READSIZE*sizeof(short int), fin);
    total += num;
    bytesRead += READSIZE*sizeof(short int);
    count++;
  }
} 

if (!FWB) {
   num = fread(wf+overlap, 1, READSIZE*sizeof(short int), fin);
  bytesRead += num;
}

  if (FWB) { for (int j=overlap; j<READSIZE*count+overlap; j++) { wf[j] = ntohs(wf[j]); } }
  if (invertWF) { for (int j=0; j<2*(READSIZE*count+overlap); j++) { wf[j] = -wf[j]; } }

  if (FWB) {
    READSIZE *= 101; // Temporary hack until I clean up code, make sure algorithms run fully
  }
    
  if (FWB) {
    return (total/sizeof(short int));
  } else {
    return (num/sizeof(short int));
  }
}

/**************************************************************/

int streamer::getSettings(TString setFile) {
  FILE *fin;
  char *p, *line, str[512] = {'0'}, str1[512] = {'0'};
  int nn = 0, nret = 0;
  int d1 = 0;
  float f1 = 0.0;
  
  if ((fin = fopen(setFile.Data(), "r")) == NULL) {
    std::cout << ALERTTEXT;
    printf("getSettings(): Could not open settings file \"%s\"\n", setFile.Data());
    std::cout << RESET_COLOR;
    return -1;
  } else {
    printf("getSettings(): Opened file \"%s\"\n", setFile.Data());
    fflush(stdout);
  }

  line = fgets(str, 512, fin);
  while (line != NULL) {
    if (str[0] == 35 || str[0] == 59 || str[0] == 10) {
      /* Do nothing -- empty or comment lines */
    } else if ((p = strstr(str, "Tau")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &f1);
      setTau(f1);
    } else if ((p = strstr(str, "Integration time")) != NULL) {
      nret = sscanf(str, "%s %s %d", str1, str1, &d1);
      setIntTime(d1);
    } else if ((p = strstr(str, "Gap (flat-top) time")) != NULL) {
      nret = sscanf(str, "%s %s %s %d", str1, str1, str1, &d1);
      setGapTime(d1);
    } else if ((p = strstr(str, "LED gap time")) != NULL) {
      nret = sscanf(str, "%s %s %s %d", str1, str1, str1, &d1);
      setLEDGapTime(d1);
    } else if ((p = strstr(str, "LED threshold")) != NULL) {
      nret = sscanf(str, "%s %s %d", str1, str1, &d1);
      setLEDThresh(d1);
    } else if ((p = strstr(str, "BLR time constant")) != NULL) {
      nret = sscanf(str, "%s %s %s %d", str1, str1, str1, &d1);
      setBLRValue(d1);
    } else if ((p = strstr(str, "BLR inhibit length")) != NULL) {
      nret = sscanf(str, "%s %s %s %d", str1, str1, str1, &d1);
      setBLRInhibit(d1);
    } else if ((p = strstr(str, "BLR retrigger")) != NULL) {
      nret = sscanf(str, "%s %s %d", str1, str1, &d1);
      setBLRTrigger(d1);
    } else if ((p = strstr(str, "Do BLR?")) != NULL) {
      nret = sscanf(str, "%s %s %d", str1, str1, &d1);
      setBLR(d1);
    } else if ((p = strstr(str, "Do fixed time pickoff for energy?")) != NULL) {
      nret = sscanf(str, "%s %s %s %s %s %s %d", str1, str1, str1, str1, str1, str1, &d1);
      setEnergyPO(d1);
    } else if ((p = strstr(str, "Fixed pickoff time")) != NULL) {
      nret = sscanf(str, "%s %s %s %d", str1, str1, str1, &d1);
      setPOTime(d1);
    } else if ((p = strstr(str, "Invert?")) != NULL) {
      nret = sscanf(str, "%s %d", str1, &d1);
      invertWF = d1;
    } else {
      /* Unknown command */
      std::cout << ALERTTEXT;
      printf("getSettings(): Unknown command \"%s\"\n", str);
      printf("               ----> Aborting!!\n");
      std::cout << RESET_COLOR;  fflush(stdout);
      return -1;
    }

    nn++; line = fgets(str, 512, fin);
  }

  fclose(fin);
  printf("getSettings(): Settings file closed.\n");
  fflush(stdout);
  return (0);

}

/**************************************************************/

void streamer::reportSettings() {
  std::cout << DCYAN;
  printf(" Analysis settings:\n");
  printf("   Tau (in 10ns clicks) = %0.2f\n", getTau());
  printf("   LED gap time = %d\n", getLEDGapTime());
  printf("   LED threshold = %d\n", getLEDThresh());
  printf("   Trapezoid = %d / %d / %d\n", getIntTime(), getGapTime(), getIntTime());
  printf("   Fixed pickoff = "); 
  if (getEnergyPO()) { printf("Yes, at %d\n", getPOTime()); }
  if (!getEnergyPO()) { printf("No.\n"); }
  printf("   Baseline restoration = ");
  if (!getBLR()) { printf("No.\n"); }
  if (getBLR()) { 
    printf("Yes...\n");
    printf("     Time constant = %d\n", getBLRValue());
    printf("     Inhibit length = %d\n", getBLRInhibit());
    printf("     Retrigger = %d\n", getBLRTrigger());
  }
  printf("\n\n");
  std::cout << RESET_COLOR;  fflush(stdout);
  
}

/**************************************************************/

int streamer::calculateLEDlevel(int index, int thresh) {
  short int d1, d2, d3, d4, d5, dd, dd1, dd2, dd3;
  index-=2;
  d1 = wf[index-2] - wf[index-LEDK-2];
  d2 = wf[index-1] - wf[index-LEDK-1];
  d3 = wf[index] - wf[index-LEDK];
  d4 = wf[index+1] - wf[index-LEDK+1];
  d5 = wf[index+2] - wf[index-LEDK+2];
  dd1 = (d1 + 2*d2 + d3)/4;
  dd2 = (d2 + 2*d3 + d4)/4;
  dd3 = (d3 + 2*d4 + d5)/4;
  dd = (dd1 + 2*dd2 + dd3)/4;
  return dd;
  
  /* Alternative - just average on each side and subtract */
  Float_t d1b = 0, d2b = 0;
  for (int i=-1; i<2; i++) {
    d1b += (Float_t)wf[index+i];
    d2b += (Float_t)wf[index-LEDK+i];
  }
  d1b /= 3.;
  d2b /= 3.;
  //printf("%d: d1b = %0.1f, d2b = %0.1f, LED = %d\n", index, d1b, d2b, (int)(d1b-d2b));
  return (int)(d1b-d2b);
}

/**************************************************************/

void streamer::doLEDfilter(int startIndex, int endIndex, long long int startTS) {

  for (int i=startIndex; i<endIndex; i++) {
    /* LED-like filter */
    if ( (i+startTS) > (LEDK + 8)) {
      ledBuf[i] = calculateLEDlevel(i, LEDThreshold);
    } else { ledBuf[i] = 0; }
  }
}

/**************************************************************/

int streamer::getLEDcrossings(int startIndex, int endIndex, long long int startTS) {
  int noiseCountdown = 0;
  for (int i=startIndex; i<endIndex; i++) {
    if (noiseCountdown>0) { noiseCountdown--; }
    if ((ledBuf[i] >= LEDThreshold) && (ledBuf[i-1] < LEDThreshold) && (noiseCountdown==0) && i>0) {
      ledOUT.push_back((long long int)i+startTS);
      noiseCountdown = LEDK;
    }
  }
  return ledOUT.size();
}


/**************************************************************/

double streamer::baseline(int index) {
  int sum = 0;
  for (int i=0; i<10; i++) {
    sum += wf[index-i];
  }
  return (double)(sum)/10.;
}

/**************************************************************/

void streamer::doTrapezoid(int startIndex, int endIndex, long long int startTS, int first) {

  /* Trapezoid filter straight out of Jordanov and Knoll, NIM A 345, 337 (1994). */
  /* Initialize for now assuming no pulse at the start of the data... */

  for (int i=startIndex; i<endIndex; i++) {
    if (i == 0 && first == 1) { trapBuf[i] = 0; }
    else if (i <= (2*EM+EK) && first == 1) { trapBuf[i] = 0; }
    else if ( (i > (2*EM+EK) && first == 1) || first != 1) {
      int scratch =  ( (wf[i]) - (wf[i-EM]) - (wf[i-EM-EK]) + (wf[i-2*EM-EK]) );
      trapBuf[i] = trapBuf[i-1] + scratch;
    }
  }
    
  /* This shouldn't be allowed, but I'm going to give a shot.
     Find the DC offset and just remove it.  */
  long long int values = 0;
  int nValue = 0;
  long long int offset = 0;
  int iLED = 0;
  if (ledOUT.size() > 0) {
    for (int j=startIndex; j<endIndex; j++) {
      if (j < ledOUT[iLED]) {
	values += trapBuf[j];
	nValue++;
      } else if (j > (ledOUT[iLED]-startTS) && j<(ledOUT[iLED]-startTS+500*(2*EK+EM))) {
	/* Nothing. */
      } else if (j==(ledOUT[iLED]-startTS+500*(2*EK+EM)) || (j==startIndex && j>(ledOUT[iLED]-startTS)) ) {
	iLED++;
      }
      if (iLED >= ledOUT.size()) { j = endIndex; }
    }
  }
  if (nValue != 0) {
    offset = values/nValue;
  } else { values = 0; }
  
  for (int i=startIndex; i<endIndex; i++) {
    if (i == 0 && first == 1) { trapBuf[i] = 0; }
    else if (i <= (2*EM+EK) && first == 1) { trapBuf[i] = 0; }
    else if ( (i > (2*EM+EK) && first == 1) || first != 1) {
      trapBuf[i] -= offset;
    }
  }
  
}

/**************************************************************/

double streamer::doPolezeroBasic(int startIndex, int endIndex, double sum, int first) {

  /* Initialization is a problem here probably */

  for (int i=startIndex; i<endIndex; i++) {
    if (i < (2*EM+EK) && first == 1) { pzBuf[i] = 0.; }
    else if (i == (2*EM+EK) && first == 1) { pzBuf[i] = 0.; }
    else if ((i > (2*EM+EK) && first == 1) || (first != 1)) {
      sum += (double)trapBuf[i];
      pzBuf[i] = trapBuf[i] + sum*(1/tau); // tau is in clock ticks
    }
  }

  return sum;
  
}

/**************************************************************/

double streamer::twoPolePolezero(int startIndex, int endIndex, double sum, int first) {

  for (int i=startIndex; i<endIndex; i++) {
    if (i < (2*EM+EK) && first == 1) { pzBuf[i] = 0.; }
    else if (i == (2*EM+EK) && first == 1) { pzBuf[i] = 0.; }
    else if ((i > (2*EM+EK) && first == 1) || (first != 1)) {
      sum += (double)trapBuf[i];
      pzBuf[i] = trapBuf[i] + sum*0.9874*(1/tau) + sum*0.0126*(1/318.215);
    }
  }
  return sum;

}

/**************************************************************/

double streamer::doLocalPZandEnergy(int startIndex, int endIndex, int LEDIndex, double tau) {

  int sumL = 0;

  double BL = 0.0;
  int nBL = 0;
  
  for (int i=startIndex; i<LEDIndex - 30; i++) {
    BL += (double)wf[i];
    nBL++;
  }
  BL /= (double)nBL;
  
  for (int i=startIndex; i<endIndex; i++) {
    if (i < (2*EM+EK)) { tempBuf[i] = 0.; }
    else if (i == (2*EM+EK)) { tempBuf[i] = 0.; }
    else if (i > (2*EM+EK)) {
      sumL += (double)wf[i] - BL;
      tempBuf[i] = wf[i-1] - BL + sumL*(1/tau);
    }
  }
  double pre = 0.0, post = 0.0;
  for (int i=0; i<EM; i++) {
    pre += tempBuf[i+startIndex];
    post += tempBuf[i+LEDIndex+100];
    //pre += tempBuf[LEDIndex-2*EM-EK-100+i];
    //post += tempBuf[LEDIndex-100-EM+i];
  }
  // return (((post-pre)/(double)EM) + 1000.0);
  return ((post-pre)/32.);
  
}

/**************************************************************/

int streamer::doBaselineRestorationCC(int startIndex, int endIndex, long long int startTS, int first, int countdown) {

  int iLED =0;
  int BLRinhibit = 0;
  //  int countdown = 0;
  int numRetriggers = 0;
  int j=0;
  long int ledTS = 0;

  double div = 1/pow(2., DV);

  for (int i=startIndex+1; i<endIndex; i++) {

    if (ledOUT.size() > 0) {
      if (i > ledOUT[iLED]-startTS && i == startIndex) { iLED++; }
      
      /* We see an LED -- inhibit for duration inhibitLength */
      if (i == (ledOUT[iLED]-startTS-100)) {
	/* Check if this is a re-trigger */
	if (countdown !=0) {
	  numRetriggers++;
	  if ((numRetriggers <= BLRretrigger) ||(BLRretrigger = -1)) {
	    BLRinhibit = 1;  
	    countdown = BLRinhibitLength;
	  } else { /* Don't retrigger anymore, just cross our fingers this
		      doesn't become a mess. */
	    BLRinhibit = 0;
	    countdown = 0;
	    numRetriggers = 0;
	  }
	  iLED++;
	} else {
	  BLRinhibit = 1;
	  countdown = BLRinhibitLength;
	  iLED++;
	}
      } else {
	if (countdown > 0) { 
	  countdown--;
	  BLRinhibit = 1;
	} else {
	  BLRinhibit = 0;
	  numRetriggers = 0;
	}
      }
    }	

    subBL1[i] = (pzBuf[i-1]-BL1[i-1])*div;
    subBL2[i] = (BL1[i-1]-BL2[i-1])*div;
    
    if (!BLRinhibit) {
      BL1[i] = BL1[i-1] + subBL1[i-1];
      BL2[i] = BL2[i-1] + subBL2[i-1];
    } else {
      BL1[i] = BL1[i-1];
      BL2[i] = BL2[i-1];
    }
  
    pzBLBuf[i] = pzBuf[i] - BL2[i];
  }

  return countdown;
}

/**************************************************************/

void streamer::doBaselineRestorationM2(int startIndex, int endIndex, long long int startTS, int first) {

  /* This corresponds to my understanding of the MatLab function "baserestoration11_v2.m", 
     but limited to 'self' inhibition only (none of this segment vs. CC inhibiting stuff). */

  int iLED =0;
  int BLRinhibit = 0, BLRinhibitP = 0;
  int countdown = 0;
  int numRetriggers = 0;

  double div = 1/pow(2., DV);
  
  for (int i=startIndex; i<endIndex; i++) {

    if (ledOUT.size() > 0) {
      if (i > ledOUT[iLED]-startTS && i == startIndex) { iLED++; }
      
      /* We see an LED -- inhibit for duration inhibitLength */
      if (i == (ledOUT[iLED]-startTS)) {
	/* Check if this is a re-trigger */
	if (countdown !=0) {
	  numRetriggers++;
	  if ((numRetriggers <= BLRretrigger) ||(BLRretrigger = -1)) {
	    BLRinhibit = 1;  
	    countdown = BLRinhibitLength;
	  } else { /* Don't retrigger anymore, just cross our fingers this
		      doesn't become a mess. */
	    BLRinhibit = 0;
	    countdown = 0;
	    numRetriggers = 0;
	  }
	  iLED++;
	} else {
	  BLRinhibit = 1;
	  countdown = BLRinhibitLength;
	  iLED++;
	}
      } else {
	if (countdown > 0) { 
	  countdown--;
	  BLRinhibit = 1;
	} else {
	  BLRinhibit = 0;
	  numRetriggers = 0;
	}
      }
    }	

    if (!BLRinhibit) {
      BL1[i] = BL1[i-1]*(1-div) + pzBuf[i];
      subBL1[i] = BL1[i]*(div);
      BL2[i] = BL2[i-1]*(1-div) + subBL1[i];
      subBL2[i] = BL2[i]*(div);
    } else {
      BL1[i] = BL1[i-1];
      subBL1[i] = subBL1[i-1];
      BL2[i] = BL2[i-1];
      subBL2[i] = subBL2[i-1];
    }
    
    BLRinhibitP = BLRinhibit;

    pzBLBuf[i] = pzBuf[i] - subBL2[i];
    
  }

}

/**************************************************************/

std::vector<double> streamer::doEnergyPeakFind(double *in, int startIndex, int endIndex, long long int startTS, int *pileUp) {

  std::vector<double> energies;
  int piledUp = 0;
  double peak = 0;

  for (int i=0; i<ledOUT.size(); i++) {
    /* Check for pileup condition... */
    if (i!=0 && i!=ledOUT.size()-1) {
      if ( ((ledOUT[i]-ledOUT[i-1])< EM+EK) ||
	   ((ledOUT[i+1]-ledOUT[i])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    } else if (i==0) { /* Edge case of first LED... */
      if ( ((ledOUT[i+1]-ledOUT[i])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    } else if (i==ledOUT.size()-1) { /* Edge case of last LED... */
      if ( ((ledOUT[i]-ledOUT[i-1])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    }
    peak = 0;
    for (int j=0; j<EM+EK; j++) {
      if (ledOUT[i]-startTS+j < READSIZE) {
	if (in[ledOUT[i]-startTS+j] > peak) {
	  peak = in[ledOUT[i]-startTS+j];
	}
      }
    }
    energies.push_back(peak/32.);
    energies.push_back(piledUp);
  }

  return energies;
}

std::vector<double> streamer::doPeakSensing(int startIndex, int endIndex, long long int startTS, int *pileUp) {
  std::vector<double> energies;
  int piledUp = 0;
  double peak = 0;
  
  for (int i=0; i<ledOUT.size(); i++) {
    peak = 0;
    for (int j=0; j<3000; j++) {
      if (ledOUT[i]-startTS+j < READSIZE) {
	if (wf[ledOUT[i]-startTS+j] > peak) {
	  peak = (double)wf[ledOUT[i]-startTS+j];
	}
      }
    }
    energies.push_back(peak - wf[ledOUT[i]-startTS-200]);
    energies.push_back(piledUp);
  }

  return energies;

}

std::vector<double> streamer::doPeakIntegrate(int startIndex, int endIndex, long long int startTS, int *pileUp) {
  std::vector<double> energies;
  int piledUp = 0;
  double peak = 0;
  int nInt = 0;
  
  for (int i=0; i<ledOUT.size(); i++) {
    peak = 0;  nInt = 0;
    for (int j=0; j<3000; j++) {
      if (ledOUT[i]-startTS+j < READSIZE && wf[ledOUT[i]-startTS+j] > wf[ledOUT[i]-startTS]) {
	peak += (double)wf[ledOUT[i]-startTS+j];
	nInt++;
      }
    }
    energies.push_back(peak - wf[ledOUT[i]-startTS-200]*nInt);
    energies.push_back(piledUp);
  }

  return energies;

}

std::vector<double> streamer::doEnergyFixedPickOff(double *in, int startIndex, int endIndex, long long int startTS, int *pileUp) {

  std::vector<double> energies;
  int piledUp = 0;
  
  for (int i=0; i<ledOUT.size(); i++) {
    /* Check for pileup condition... */
    if (i!=0 && i!=ledOUT.size()-1) {
      if ( ((ledOUT[i]-ledOUT[i-1])< EM+EK) ||
	   ((ledOUT[i+1]-ledOUT[i])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    } else if (i==0) { /* Edge case of first LED... */
      if ( ((ledOUT[i+1]-ledOUT[i])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    } else if (i==ledOUT.size()-1) { /* Edge case of last LED... */
      if ( ((ledOUT[i]-ledOUT[i-1])< EM+EK) ) {
	(*pileUp)++; piledUp = 1;
      } else { piledUp = 0; }
    }
    /* Push back the energies and the pileup, in that order */
    if ((ledOUT[i]-startTS+POtime) < READSIZE) {
      energies.push_back((in[ledOUT[i]-startTS + POtime])/32.);
    } else { energies.push_back(-1); }
    energies.push_back(piledUp);
  }
  
  return energies;
}
