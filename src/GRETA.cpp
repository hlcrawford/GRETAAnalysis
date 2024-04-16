#include "GRETA.h"

#define DEBUG 0
#define KEEP_WAVEFORMS 0

ClassImp(gretaWaveformMsg);

/* Function for 64-bit numbers endian-ness stuff */
uint64_t GRETA::ntoh64(uint64_t input) {
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

uint64_t GRETA::hton64(uint64_t input) {
  return (ntoh64(input));
}

void g3CrystalEvent::Clear() {
  rhSequence = 0;
  for (Int_t i=0; i<4; i++) { ccE[i] = 0; }
  for (Int_t i=0; i<36; i++) { segE[i] = 0; }
}

void g3OUT::Reset() {
  UInt_t crystals = crystalMult();
  xtals.clear();
}

UInt_t g3OUT::crystalMult() { return xtals.size(); }

Float_t g3OUT::calorimeterE() {
  Float_t sum = 0.;
  for (UInt_t ui = 0; ui<crystalMult(); ui++) {
    sum += xtals[ui].ccE[0];
  }
  return sum;
}

void GRETA::Initialize() { ; }

void GRETA::Reset() { 
  g3X.Clear();
  g3Out.Reset();
}

Int_t GRETA::getMode3(FILE *inf, Int_t evtLength, Int_t subType) {
  gretaWaveformMsg wform;
  Int_t siz = 0;
  
  siz = fread(&wform, 1, evtLength, inf);
  
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
    std::cin.get(); 
  } 
  
  g3X.Clear();
  g3X.version = wform.version;
  g3X.id = wform.id;
  g3X.trLen = ntohs(wform.tr_len);
  g3X.trigSrc = ntohs(wform.trig_src);
  g3X.pad = ntohs(wform.pad);
  g3X.timestamp = ntoh64((uint64_t)wform.timestamp);
  g3X.pileup = ntoh64((uint64_t)wform.pileup);
  for (Int_t i=0; i<40; i++) {
    g3X.ener[i] = ntohl(wform.ener[i]);
    if (KEEP_WAVEFORMS) {
      for (Int_t j=0; j<512; j++) {
	wform.tr[i][j] = ntohs(wform.tr[i][j]);	 
	g3X.tr[i][j] = wform.tr[i][j];
      }
    }
  }
  g3X.t0 = ntohs(wform.t0);
  g3X.subt0 = ntohs(wform.sub_t0);
  g3X.tLEDCore = ntohs(wform.t_led_core);
  g3X.tCFDCore = ntohs(wform.t_cfd_core);
  g3X.tLEDFirst = ntohs(wform.t_led_first);
  g3X.ccE[0] = (g3X.ener[segMap.find(36)->second]/32.)*gain[g3X.id][36] + offset[g3X.id][36];
  g3X.ccE[1] = (g3X.ener[segMap.find(37)->second]/32.)*gain[g3X.id][37] + offset[g3X.id][37];
  g3X.ccE[2] = (g3X.ener[segMap.find(38)->second]/32.)*gain[g3X.id][38] + offset[g3X.id][38];
  g3X.ccE[3] = (g3X.ener[segMap.find(39)->second]/32.)*gain[g3X.id][39] + offset[g3X.id][39];
  for (Int_t i=0; i<36; i++) {
    g3X.segE[i] = g3X.ener[segMap.find(i)->second]*gain[g3X.id][i] + offset[g3X.id][i];
  }
  
  //printf("Made it to the end of getMode3, just the vector push_back now.\n");

  analyzeMode3(&g3X);
  g3Out.xtals.push_back(g3X);
  return 0;
  
}

Int_t GRETA::analyzeMode3(g3CrystalEvent *g3) {

  return 0;

}
