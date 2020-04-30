#include "GRETINA.h"

ClassImp(globalHeader);
ClassImp(mode3DataPacket);

/**************************************************************/
/* g3CrystalEvent Class Functions *****************************/
/**************************************************************/

/*! Puts channels from a mode3 crystal event in numerical order
    within the vector chn<g3ChannelEvent>

    \return No return value -- directly alters class structures
*/

void g3CrystalEvent::OrderChannels() {
  g3ChannelEvent temp;
  Bool_t finished = 0;
  while (!finished) {
    finished = 1;
    for (UInt_t ui=0; ui<chn.size()-1; ui++) {
      if (chn[ui].chnNum() > chn[ui].chnNum()) {
	temp = chn[ui];
	chn[ui] = chn[ui+1];
	chn[ui+1] = temp;
	finished = 0;
      }
    }
  }
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the calibrated CC1 energy
*/

Float_t g3CrystalEvent::cc1() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the calibrated CC2 energy
*/

Float_t g3CrystalEvent::cc2() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the calibrated CC3 energy
*/

Float_t g3CrystalEvent::cc3() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the calibrated CC4 energy
*/

Float_t g3CrystalEvent::cc4() { 
 for (UInt_t ui=0; ui<chn.size(); ui++) {
   if (chn[ui].chnNum() == 39) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC1 energy
*/

Float_t g3CrystalEvent::cc1Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC2 energy
*/

Float_t g3CrystalEvent::cc2Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC3 energy
*/

Float_t g3CrystalEvent::cc3Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC4 energy
*/

Float_t g3CrystalEvent::cc4Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 39) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the waveform calculated CC1 energy
*/

Float_t g3CrystalEvent::cc1Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eCalc; }
  }
  return -1.0; 
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the waveform calculated CC2 energy
*/

Float_t g3CrystalEvent::cc2Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the waveform calculated CC3 energy
*/

Float_t g3CrystalEvent::cc3Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the waveform calculated CC4 energy
*/

Float_t g3CrystalEvent::cc4Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 39) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Calculates a raw (NOT Dino corrected) segment sum for the crystal, 
    with no threshold on the segments included

    \return Returns the float value of the raw segment sum for
            the crystal 
*/

Float_t g3CrystalEvent::segSumRaw() {
  return segSumRawThresh(0.0);
}

/**************************************************************/

/*! Calculates a raw (NOT Dino corrected) segment sum for the crystal,
    including segments that are above a threshold value

    \param thresh The float energy threshold; any segments below
           this value are not included in the sum
    \return Returns the float value of the raw segment sum for the crystal
*/

Float_t g3CrystalEvent::segSumRawThresh(Float_t thresh) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      segSum += chn[ui].eCal;
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumRawTimed(Float_t minT, Float_t maxT) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36) {
      if (chn[ui].calcTime > minT &&
	  chn[ui].calcTime > maxT) {
	segSum += chn[ui].eCal;
      }
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumRawThreshTimed(Float_t thresh, Float_t minT, Float_t maxT) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      if (chn[ui].calcTime > minT &&
	  chn[ui].calcTime > maxT) {
	segSum += chn[ui].eCal;
      }
    }
  }
  return segSum;
}

/**************************************************************/

/*! Determines the number of segments hit in the crystal above a given
    energy threshold

    \param thresh The float energy threshold -- segments with energy
           depositions above threshold are counted as hit
    \return Returns integer value of the number of segments hit above 
            the energy threshold
*/

Int_t g3CrystalEvent::segsHit(Float_t thresh) {
  Int_t hit = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      hit++;
    }
  }
  return hit;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumCalc() {
  Float_t segSum = 0.;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36) {
      segSum += chn[ui].eCalc;
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumCalcThresh(Float_t thresh) {
  Float_t segSum = 0.;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36 && chn[ui].eCalc < thresh) {
      segSum += chn[ui].eCalc;
    }
  }
  return segSum;
}

/**************************************************************/

Int_t g3CrystalEvent::segsHitCalc(Float_t thresh) {
  Int_t hit = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > thresh && chn[ui].segNum < 36) {
      hit++;
    }
  }
  return hit;
}

/**************************************************************/

Float_t g3CrystalEvent::maxSegE() {
  Float_t maxE = 0.;  Int_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCal;
    }
  }
  return maxE;
}

/**************************************************************/

Float_t g3CrystalEvent::maxSegECalc() {
  Float_t maxE = 0.;  Int_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCalc;
    }
  }
  return maxE;
}

/**************************************************************/

Int_t g3CrystalEvent::maxSegNum() {
  Float_t maxE = 0.;  Int_t max = -1;
 
  for (Int_t u=0; u<chn.size(); u++) {
    if (chn[u].eCal > maxE && chn[u].segNum < 36) {
      max = u;
      maxE = chn[u].eCal;
    }
  }

  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::secondSegNum() {
  Float_t maxE = 0.;  Int_t max = -1;
  
  for (Int_t u=0; u<chn.size(); u++) {
    if (chn[u].eCal > maxE && chn[u].segNum < 36 && chn[u].eCal != maxSegE()) {
      max = u;
      maxE = chn[u].eCal;
    }
  }

  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::maxSegNumCalc() {
  Float_t maxE = 0.;  UInt_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCalc;
    }
  }
  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::deepRingHit(Float_t thresh) {
  Int_t deep = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 6 && chn[ui].eCal > thresh && deep < 1) {
      deep = 1; 
    } else if (chn[ui].segNum < 12 && chn[ui].segNum > 5 && 
	       chn[ui].eCal > thresh && deep < 2) {
      deep = 2;
    } else if (chn[ui].segNum < 18 && chn[ui].segNum > 11 && 
	       chn[ui].eCal > thresh && deep < 3) {
      deep = 3;
    } else if (chn[ui].segNum < 24 && chn[ui].segNum > 17 && 
	       chn[ui].eCal > thresh && deep < 4) {
      deep = 4;
    } else if (chn[ui].segNum < 30 && chn[ui].segNum > 23 && 
	       chn[ui].eCal > thresh && deep < 5) {
      deep = 5;
    } else if (chn[ui].segNum < 36 && chn[ui].segNum > 29 && 
	       chn[ui].eCal > thresh && deep < 1) {
      deep = 6;
    }
  }
  return deep;
}

/**************************************************************/

long long int g3CrystalEvent::LEDLow() {
  long long int ledLow = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (ledLow == -1) { ledLow = chn[ui].timestamp; }
    if (ledLow > chn[ui].timestamp) {
      ledLow = chn[ui].timestamp;
    }
  }
  return ledLow;
}

/**************************************************************/

long long int g3CrystalEvent::LEDHigh() {
  long long int ledHigh = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (ledHigh < chn[ui].timestamp) {
      ledHigh = chn[ui].timestamp;
    }
  }
  return ledHigh;
}

/**************************************************************/

UInt_t g3CrystalEvent::crystalBuild() { return chn.size(); }

/**************************************************************/

/*! Function calculates the difference between the high and low 
    LED timestamps for channel events within the mode3 crystal event,
    essentially the 'width' of the event for the crystal

    \return Returns the difference between the first and last LED
            timestamps as a long long int
*/

long long int g3CrystalEvent::LEDRange() {
  return (LEDHigh() - LEDLow());
}

/**************************************************************/

void g3CrystalEvent::mode3ToFile(FILE *outFile) {

  mode3DataPacket *dp;
  dp = (mode3DataPacket*)malloc(sizeof(dp->aahdr) + 
				sizeof(dp->hdr) + 
				sizeof(dp->waveform));
  memset(dp->waveform, 1, MAX_TRACE_LENGTH * sizeof(UShort_t));

  globalHeader *gHeader;
  gHeader = (globalHeader*)malloc(sizeof(gHeader));

  unsigned char outBuf[32*1024];
  unsigned char *tmp = (outBuf);

  dp->aahdr[0] = 0xAAAA;
  dp->aahdr[1] = 0xAAAA;
  
  for (Int_t i=0; i<chn.size(); i++) {

    gHeader->type = RAW;
    gHeader->timestamp = chn[i].timestamp;
    Int_t evtLength = (sizeof(dp->aahdr) + sizeof(dp->hdr) + chn[i].tracelength()*sizeof(UShort_t));
    gHeader->length = evtLength;

    dp->hdr[0] = chn[i].hdr0;
    dp->hdr[1] = chn[i].hdr1;
    dp->hdr[2] = ((chn[i].timestamp & 0xffff0000) >> 16);
    dp->hdr[3] = ((chn[i].timestamp & 0xffff));
    dp->hdr[5] = ((chn[i].timestamp & 0xffff00000000) >> 32);
    dp->hdr[6] = ((chn[i].CFDtimestamp & 0xffff));
    dp->hdr[8] = ((chn[i].CFDtimestamp & 0xffff00000000) >> 32);
    dp->hdr[9] = ((chn[i].CFDtimestamp & 0xffff0000) >> 16);

    Bool_t sign = 0;
    Int_t tmpIntEnergy = (Int_t)(32*chn[i].eRaw);
    UInt_t tmpEnergy = 0;
    if (chn[i].chanID()%10 == 9) { /* CC */
      if (chn[i].eRaw < 0) {
	tmpEnergy = (tmpIntEnergy + (UInt_t)(0x01000000));
	sign = 1;
      } else {
	tmpEnergy = (UInt_t)(tmpIntEnergy);
      }
    } else {
      if (chn[i].eRaw < 0) {
	tmpEnergy = (UInt_t)(-tmpIntEnergy);
      } else {
	tmpEnergy = (UInt_t)((0x01000000)-tmpIntEnergy);
	sign = 1;
      }
    }
    dp->hdr[4] = (tmpEnergy & 0xffff);
    dp->hdr[7] = ((tmpEnergy & 0x00ff0000) >> 16);
    if (sign) { dp->hdr[7] = (dp->hdr[7] & 0x0100); }
    dp->hdr[7] += (chn[i].pileUp() << 15);

    for (Int_t j=0; j<chn[i].tracelength()+1; j=j+2) {
      if (chn[i].wf.raw[j] < 0) {
    	dp->waveform[j+1] = (UShort_t)(chn[i].wf.raw[j] + std::numeric_limits<unsigned int>::max());
      } else {
    	dp->waveform[j+1] = (UShort_t)(chn[i].wf.raw[j]);
      }
      if (chn[i].wf.raw[j+1] < 0) {
    	dp->waveform[j] = (UShort_t)(chn[i].wf.raw[j+1] + std::numeric_limits<unsigned int>::max());
      } else {
    	dp->waveform[j] = (UShort_t)(chn[i].wf.raw[j+1]);
      }
    }

    memmove(tmp, &dp->aahdr[0], evtLength);
    for (Int_t j=0; j<evtLength; j=j+2) {
      swap(*(outBuf + j), *(outBuf + j + 1));
    }
    
    fwrite(gHeader, sizeof(struct globalHeader), 1, outFile);
    fwrite(tmp, evtLength, 1, outFile);
    
  }

}

/**************************************************************/
/* g3OUT Class Functions **************************************/
/**************************************************************/

/*! Resets the mode3 output class -- clears channel vector

    \return No return value -- clears data structures in class
            directly
*/

void g3OUT::Reset() {
  UInt_t crystals = crystalMult();
  for (UInt_t ui=0; ui<crystals; ui++) {
    xtals[ui].chn.clear();
  }
  xtals.clear();
}

/**************************************************************/

/*! Extracts the number of crystals hit within a mode3 event, namely
    the size of the xtals vector

    \return Returns the unsigned integer value of the size of the
            xtals array
*/

UInt_t g3OUT::crystalMult() { return xtals.size(); }

/**************************************************************/

Float_t g3OUT::calorimeterE() {
  Float_t sum = 0.;
  for (UInt_t ui=0; ui<crystalMult(); ui++) {
    sum += xtals[ui].cc1();
  }
  return sum;
}
