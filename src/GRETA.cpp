#include "GRETA.h"

#define DEBUG 0
#define KEEP_WAVEFORMS 0

ClassImp(gretaWaveformMsg);

/*********************************************************************************/

/* This next section of code is basically taken from the dptc_unpack files,
   because it didn't work when I tried to just include the headers etc.
   So I reproduced the whole thing to make the compressed wavefunctions work. */

#define DPTC_IMPL_UNPACK  dptc_unpack16
#define DPTC_OUTPUT_TYPE  uint16_t

#include <unistd.h>

#define DPTC_INT_DEBUG(...)  do { } while (0)

/* If the read pointer has not passed the end of the input buffer,
 * read the next word.  This is protect against buffor overflow on
 * corrupt data that.
 *
 * Then shift the data bits to the correct position in the buffer, and
 * mix them in.
 *
 * If the buffer actually had less than 32 bits available, increment
 * the read pointer.
 */

#define DPTC_INT_FETCH_WORD(src, src_end, src_val,	\
			    buffer, avail)		\
  do {							\
    src_val = 0;					\
    if (((ssize_t) (src - src_end)) < 0)		\
      src_val = ((uint64_t) *src);			\
    buffer |= src_val << avail;				\
    /* The following is equivalent to: */		\
    /* if (avail < 32) { avail += 32; src++; }*/	\
    src += ((~avail) & 0x20) >> 5;			\
    avail |= 0x20;					\
  } while (0)


int dptc_unpack16(uint32_t *compr, size_t ncompr,
		     uint16_t *output, size_t ndata,
		     int bits )
{
  /* *** Constants based on the number of bits. *** */

  /* Mask for the data bits. */
  uint32_t mask = ((((uint32_t) 1) << bits)-1);
  /* Signbit for the data. */
  uint32_t signbit   = mask ^ (mask >> 1);

  /* Number of bits required to store the number of bits used in a chunk. */
  /* Only allowing 'bits' values up to 16, we can do log2 with a lookup. */
  uint32_t storebits = (bits-3-1 < 16) ?
    ((0x4444444433332211ll >> ((bits-3-1)*4)) & 0xf) : 5;
  /* Mask for the same. */
  uint32_t storebitsmask = ((((uint32_t) 1) << storebits)-1);

  /* *** Data source pointers. *** */ 

  /* Pointer to next word to use. */
  uint32_t *src = compr;
  /* Pointer to one after last word to use.  (Will never be read). */
  uint32_t *src_end = compr + ncompr;

  /* *** Local buffer of next bits to use. *** */

  /* Compressed data.  This buffer always contain at least 32 bits of data. */
  uint64_t buffer = 0;
  /* Number of bits currently available. */
  uint32_t avail = 0;
  /* Value of next word in data source. */
  uint64_t src_val;

  /* *** Chunk handling. *** */

  /* Number of bits stored for this chunk. */
  uint32_t numbits;
  /* Temporary value to adjusted wrapped values. */
  uint32_t numbits_wrap;
  /* Number of bits stored in previous chunk. */
  uint32_t prevnumbits;
  /* Loop variables. */
  size_t i, j;

  /* *** Compression handling. *** */
  
  /* Mask for data with number of bits of this chunk. */
  uint32_t storedmask;
  /* Base value (most negative possible) for data within this chunk. */
  uint32_t storedbase;
  /* Actual value stored (recovered after mask and base add). */
  uint32_t stored;

  /* Difference value stored might have had its sign flipped. */
  uint32_t diff_flipped;
  /* Keep track if value is stored with sign change? */
  int flipsign = 0;

  /* Difference value after sign flip correction. */
  uint32_t diff;
  /* Keep track of last three differences. */
  uint32_t d1 = 0, d2 = 0, d3 = 0;

  /* Previous value.  (To undo difference operation.) */
  uint32_t prev = 0;

  if (bits > 8 * sizeof (DPTC_OUTPUT_TYPE))
    {
      DPTC_INT_DEBUG("Cannot handle more bits (%d) than output type (%d).\n",
		     bits, 8 * sizeof (DPTC_OUTPUT_TYPE));
      return 1;
    }

  DPTC_INT_DEBUG("-----------------------------------------------------\n");

  DPTC_INT_FETCH_WORD(src, src_end, src_val, buffer, avail);
  
  if (ndata)
    {
      prev = (buffer & mask) ^ signbit;
      d1 = prev;
      flipsign = (prev & signbit);

      *(output++) =
	prev;

      DPTC_INT_DEBUG(" : %04x\n", prev);

      buffer >>= bits;
      avail -= bits;
    }
  
  prevnumbits = bits;

  for (i = 1; i < ndata; )
    {
      DPTC_INT_FETCH_WORD(src, src_end, src_val, buffer, avail);

      size_t numstore = DPTC_CHUNK_SIZE;

      if (numstore > ndata - i)
	numstore = ndata - i;

      /* Get the number of bits. */

      numbits = buffer & 0x03;

      buffer >>= 2;
      avail -= 2;

      /* This branch is taken very seldomly.  It is faster to keep it
       * as a branch than to evaluate both paths and assign results
       * using bitmasks.
       */
      if (numbits == 0)
	{
	  numbits = (buffer & storebitsmask);
	  DPTC_INT_DEBUG("long bits: +%d+2 [%d]\n", numbits, storebits);
	  numbits = prevnumbits + numbits + 2;
	  buffer >>= storebits;
	  avail -= storebits;
	}
      else
	{
	  DPTC_INT_DEBUG("short bits: +%d-2\n", numbits);
	  numbits = prevnumbits + numbits + bits - 2;
	}
      DPTC_INT_DEBUG("temp bits: %d\n", numbits);
      /* This version is generally faster than the modulo operation. */
      numbits_wrap = numbits - bits;
      if (numbits > bits)
	numbits = numbits_wrap;
      numbits_wrap = numbits - bits;
      if (numbits > bits)
	numbits = numbits_wrap;

      storedmask = (1 << numbits) - 1;
      storedbase = -((1 << numbits) >> 1);

      DPTC_INT_DEBUG("numbits: %d (mask 0x%04x, base 0x%04x)\n",
		     numbits, storedmask, storedbase);

      for (j = 0; j < numstore; j++)
	{
	  /* Fetch more bits if needed. */
	  DPTC_INT_FETCH_WORD(src, src_end, src_val, buffer, avail);	  

	  /* Get data. */

	  stored = buffer & storedmask;

	  DPTC_INT_DEBUG("%04x+b=%04x", stored, (stored + storedbase) & mask);
	  
	  stored += storedbase;

	  buffer >>= numbits;
	  avail -= numbits;

	  diff_flipped = stored;

#if DPTC_ENABLE_FLIPSIGN
	  /* 0- instead of - to avoid Visual2019 error. */
	  diff = flipsign ? (0-diff_flipped) : diff_flipped;

	  if (diff)
	    flipsign = (diff & signbit);
#else
	  diff = diff_flipped;
#endif

	  DPTC_INT_DEBUG(" unflipped %04x [%d] ",
			 diff & mask, flipsign ? 1 : 0);

	  /* Do the three last differences have the same sign?
	   * A zero value breaks the sequence.
	   * As we never use the full range of an 32-bit int,
	   * we can look for all three values being postive by the
	   * combinded sign bit of their negative representations.
	   */

	  /* It has been tested to first run the unpacking, store the
	   * intermediate differences, and then do the final predictor
	   * calculations.  That is slower.
	   */

#if DPTC_ENABLE_PREDICTOR
	  if (((d1 & d2 & d3) | ((-d1) & (-d2) & (-d3))) & signbit)
	    {
	      DPTC_INT_DEBUG("[%d %d %d] ", d1, d2, d3);
	      diff += d1;
	    }

	  d3 = d2;
	  d2 = d1;
	  d1 = diff;
#endif

	  prev = (prev + diff) & mask;

	  *(output++) =
	    prev;

	  DPTC_INT_DEBUG("%04x : %04x\n", diff & mask, prev & mask);

	  /* There is a speed advantage of explicitly unrolling this
	   * loop.  We gain by only fetching into the buffer for every
	   * two values.  Possible as long as we only support 16-bit
	   * values.
	   */
	  if (++j >= numstore)
	    break;

	  if (sizeof (DPTC_OUTPUT_TYPE) > 2)
	    {
	      /* We might consume more than 16 bits per cycle,
	       * fetch more bits if needed.
	       */
	      DPTC_INT_FETCH_WORD(src, src_end, src_val, buffer, avail);
	    }

	  /* Get data. */

	  stored = buffer & storedmask;

	  DPTC_INT_DEBUG("%04x+b=%04x", stored, (stored + storedbase) & mask);
	  
	  stored += storedbase;

	  buffer >>= numbits;
	  avail -= numbits;

	  diff_flipped = stored;

#if DPTC_ENABLE_FLIPSIGN
	  diff = flipsign ? (0-diff_flipped) : diff_flipped;

	  if (diff)
	    flipsign = (diff & signbit);
#else
	  diff = diff_flipped;
#endif

	  DPTC_INT_DEBUG(" unflipped %04x [%d] ",
			 diff & mask, flipsign ? 1 : 0);

#if DPTC_ENABLE_PREDICTOR
	  if (((d1 & d2 & d3) | ((-d1) & (-d2) & (-d3))) & signbit)
	    {
	      DPTC_INT_DEBUG("[%d %d %d] ", d1, d2, d3);
	      diff += d1;
	    }

	  d3 = d2;
	  d2 = d1;
	  d1 = diff;
#endif
      
#if DPTC_IMPL_ADAPTIVE_DOWNSAMPLING
	  l_shift = (step_plus_1 >> 1);
	  r_shift = (1 >> step_plus_1);

	  prev = (((prev << l_shift) + diff) >> r_shift) & mask;
	  DPTC_INT_DEBUG("[r:%d l:%d]\n", r_shift, l_shift);
#else
	  prev = (prev + diff) & mask;
#endif

	  *(output++) =
	    prev;

	  DPTC_INT_DEBUG("%04x : %04x\n", diff & mask, prev & mask);
	}

      i += numstore;

      prevnumbits = numbits;
    }

  DPTC_INT_DEBUG("%zd \n", src - compr);

  if (buffer)
    {
      DPTC_INT_DEBUG("non-zero bits left in buffer!: 0x%016llx\n", buffer);
      return 1;
    }

  if (((ssize_t) (src - src_end)) < 0)
    {
      DPTC_INT_DEBUG("words left in input: %zd\n", src_end - src);
      return 1;
    }

  return 0;
}

#undef DPTC_INT_DEBUG

/*********************************************************************************/

/* Now back to GRETA-specific code */


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

void GRETA::readGRETACalibration(TString filename) {
  FILE *fp;
  Int_t i1, i2, nn;  Float_t f1, f2;
  char *st, str[128];

  /* Open file */
  fp = fopen(filename.Data(), "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", filename.Data());
    exit(1);
  }
  printf("\"%s\" open...", filename.Data());

  nn = 0;
  st = fgets(str, 64, fp);
  while (st != NULL) {
    if (str[0] == 35 || str[0] == 59 || str[0] == 10) {
      /* Comment or blank line, do nothing. */
    } else {
      sscanf(str, "%i %i %f %f\n", &i1, &i2, &f1, &f2);
      if(i1 >= 0 && i1 <= MAXCRYSTALS && i2 >= 0 && i2 < NUM_CHAN) {
	gain[i1][i2] = f1;
	offset[i1][i2] = f2;
      }
      nn++;
    }

    st = fgets(str, 64, fp);
  }
  printf("Read %i calibrations.\n", nn);

  fclose(fp);
}

void g3CrystalEvent::Clear() {
  rhSequence = 0;
  ccE.clear();
  segE.clear();
  segsHit.clear();
  tr.clear();
}

TGraph* g3CrystalEvent::GetSegTrace(Int_t seg) {
  vector<int> xVal;
  vector<int> yVal;
  for (Int_t i=0; i<trLen; i++) {
    xVal.push_back(i);
    yVal.push_back((Int_t)(tr[seg][i]));
  }
  TGraph *trSeg = new TGraph(trLen, xVal.data(), yVal.data());
  return trSeg;
}

TGraph* g3CrystalEvent::GetFullTrace() {
  vector<int> xVal;
  vector<int> yVal;
  for (Int_t seg=0; seg<40; seg++) {
    for (Int_t i=0; i<trLen; i++) {
      xVal.push_back(seg*trLen + i);
      yVal.push_back((Int_t)(tr[seg][i]));
    }
  }
  TGraph *trAll = new TGraph(trLen*40, xVal.data(), yVal.data());
  return trAll;
}

void g3OUT::Reset() {
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

gHistos::gHistos() {
  for (Int_t i=0; i<121; i++) {
    for (Int_t j=0; j<4; j++) {
      eRawCC[i][j] = NULL;
    }
    for (Int_t j=0; j<36; j++) {
      eRawSeg[i][j] = NULL;
    }
  }
}

gHistos::~gHistos() { ; }

void gHistos::writeHistos() {
  for (Int_t i=0; i<121; i++) {
    for (Int_t j=0; j<4; j++) {
      if (eRawCC[i][j]!=NULL) { eRawCC[i][j]->Write(); eRawCC[i][j]->Delete(); }
    }
    for (Int_t j=0; j<36; j++) {
      if (eRawSeg[i][j]!=NULL) { eRawSeg[i][j]->Write(); eRawSeg[i][j]->Delete(); }
    }
  }    
}

void GRETA::fillHistos() {
  for (UInt_t xtal=0; xtal<g3Out.xtals.size(); xtal++) {
    Int_t id = g3Out.xtals[xtal].rhSubType;
    for (Int_t i=0; i<4; i++) {
      if (!gHist.eRawCC[id][i]) {
	gHist.eRawCC[id][i] = new TH1F(Form("eRawCC_%d_%d", id, i), Form("eRawCC_%d_%d", id, i), 10000, 0, 3000/gain[id][i+36]);
      } else { ; }
      gHist.eRawCC[id][i]->Fill(g3Out.xtals[xtal].ener[i+36]/128.);
    }
    for (Int_t i=0; i<36; i++) {
      if (!gHist.eRawSeg[id][i]) {
	gHist.eRawSeg[id][i] = new TH1F(Form("eRawSeg_%d_%d", id, i), Form("eRawSeg_%d_%d", id, i), 10000, 0, 3000/gain[id][i+36]);
      } else { ; }
      gHist.eRawSeg[id][i]->Fill(g3Out.xtals[xtal].ener[i]/128.);
    }					    
  }
}
      
void GRETA::Initialize(controlVariables *ctrl) {
  ng2 = 0;

  for (Int_t i=0; i<121; i++) {
    lastSeqNum[i] = -1;
    seqNumOOO[i] = 0;
    seqNumSkipped[i] = 0;
    timesSeqNumSkipped[i] = 0;
    eventCnt[i] = 0;
    mode3TooLong[i] = 0;
  }
  
}

void GRETA::Reset() { 
  g3X.Clear();
  g3Out.Reset();
}

Int_t GRETA::getMode3(FILE *inf, Int_t evtLength, Int_t subType, Int_t type, controlVariables *ctrl) {
  gretaWaveformMsg wform;
  Int_t siz = 0;

  if (DEBUG) {
    printf("Entering getMode3: %d %d %d\n", evtLength, subType, type);
  }

  if (evtLength > 9000) {
    if (DEBUG) {
      printf("GetMode3: Long event length %d - over frame size.\n", evtLength);
      printf("Will skipping waveforms, under assumption it is not all here.\n");
    }
    mode3TooLong[subType]++;
  }
    
  siz = fread(&wform, 1, sizeof(struct gretaWaveformMsg), inf);
  
  if (DEBUG) {
    printf("Waveform message!\n"); 
    printf("  Version:  %i\n", wform.version); 
    printf("  ID:           %i\n", wform.id);
    printf("  tr_len:      %d\n", ntohs(wform.tr_len)); 
    printf("  trig_src:   %i\n", ntohs(wform.trig_src)); 
    printf("  pad:         %i\n", ntohs(wform.pad));  
    printf("  timestamp: %lld\n", (int64_t)ntoh64((uint64_t)wform.timestamp)); 
    for (int i=0; i<40; i++) { 
      printf("  Raw energy %d: %d\n", i, ntohl(wform.ener[i]));
    }  
    std::cin.get(); 
  } 

  g3X.Clear();
  /* Minimal unpack for superpulse analysis...*/
  g3X.id = subType;
  g3X.trLen = ntohs(wform.tr_len);
  for (Int_t i=0; i<40; i++) {
    g3X.ener[i] = ntohl(wform.ener[i]);
  }
  if (ctrl->superPulse)  {  sp.trLen = g3X.trLen; }
  for (Int_t i=36; i<40; i++) {
    g3X.ccE.push_back((g3X.ener[i]/128.)*gain[subType][i] + offset[subType][i]);
  }
  for (Int_t i=0; i<36; i++) {
    g3X.segE.push_back((g3X.ener[i]/128.)*gain[subType][i] + offset[subType][i]);
  }

  if (evtLength <= 9000) {
    vector<int16_t> itr;
    if (ctrl->keepWaveforms || ctrl->superPulse) {
      if (type==3) { // Uncompressed waveforms, type 3, need to propagate this forward
	int16_t trace[5000];
	if (DEBUG) {
	  printf("Getting uncompressed waveforms.\n");
	}
	for (Int_t i=0; i<40; i++) {
	  if (DEBUG) { printf("Segment %d\n", i); }
	  siz = fread(trace, 2, g3X.trLen, inf);
	  for (Int_t j=0; j<g3X.trLen; j++) {
	    itr.push_back(ntohs(trace[j]));
	    if (DEBUG) {
	      printf("  %d  %d\n", j, itr.back());
	    }
	  }
	  g3X.tr.push_back(itr);
	  itr.clear();
	}
      } else if (type==4) { // Compressed waveforms, type 4
	/* Get total compressed length of waveforms */
	Int_t sumCompressedLengths = 0;
	uint16_t sumDecompressionErrs = 0;
	for (Int_t i=0; i<40; i++) {
	  sumCompressedLengths += ntohs(wform.wfLength[i]);
	}
	if (sumCompressedLengths*4 != evtLength - sizeof(struct gretaWaveformMsg)) {
	  printf("Compressed waveform lengths do not match expected, %u != %lu\n", sumCompressedLengths*4, evtLength-sizeof(struct gretaWaveformMsg));
	}
	uint32_t wfData[16400];
	fread(wfData, 4, sumCompressedLengths, inf);
	for (Int_t i=0; i<sumCompressedLengths; i++) {
	  wfData[i] = ntohl(wfData[i]);
	}
	uint32_t *temp32 = (wfData);
	for (Int_t ch=0; ch<40; ch++) {
	  int16_t trTemp[5000] = {0};
	  unsigned char* tmp = (unsigned char*)(trTemp);
	  sumDecompressionErrs += dptc_unpack16(temp32, ntohs(wform.wfLength[ch]), (uint16_t*)( &trTemp ), g3X.trLen, 16);
	  for (Int_t l=0; l<g3X.trLen; l++) {
	    itr.push_back(trTemp[l]);
	  }
	  g3X.tr.push_back(itr);
	  itr.clear();
	  temp32 += ntohs(wform.wfLength[ch]);      
	}
      }
    } else {
      fseek(inf, (evtLength-sizeof(struct gretaWaveformMsg)), SEEK_CUR);
    }
  } else {
    fseek(inf, (9000-sizeof(struct gretaWaveformMsg)), SEEK_CUR);
  }

  if (ctrl->superPulse) {
    for (Int_t j=0; j<sp.trLen; j++) {
      for (Int_t seg=0; seg<36; seg++) {
	sp.waves[subType][seg][j] = g3X.tr[seg][j]; 
      }
      sp.waves[subType][36][j] = g3X.tr[36][j];
    }
  }
  
  g3Out.xtals.push_back(g3X);
  return 0;
  
}

Int_t GRETA::checkSP() {
  for (Int_t i=0; i<121; i++) {
    sp.netSeg[i] = -1;
    sp.mult[i] = 0;
    sp.ccE[i] = -1;
    sp.segE[i] = -1;
    
  }

  Int_t xtalNum = -1;
  
  for (UInt_t xtal = 0; xtal<g3Out.xtals.size(); xtal++) {
    xtalNum = g3Out.xtals[xtal].id;
    /* Loop over signals */
    for (Int_t i=0; i<37; i++) {
      /* Adjust trace baseline offsets */
      Int_t s=0;
      for (Int_t b=0; b<25; b++) {
	s += sp.waves[xtalNum][i][b];
      }
      if (s >= 0) { s = (s+7)/25; }
      else { s = (s-7)/25; }
      for (Int_t b=0; b<sp.trLen; b++) {
	sp.waves[xtalNum][i][b] -= s;
      }
    }

     /* Check for net energy. */
    for (Int_t i=0; i<36; i++) {
      Int_t avg = 0;
      for (Int_t b=sp.trLen-10; b<sp.trLen; b++) {
	avg += TMath::Abs(sp.waves[xtalNum][i][b]);
      }
      avg /= 10;
      
      if (avg > 20 && g3Out.xtals[xtal].segE[i]) { /* Net */
	if (g3Out.xtals[xtal].segE[i] > sp.segE[xtalNum]) {
	  sp.segE[xtalNum] = g3Out.xtals[xtal].segE[i];
	  sp.netSeg[xtalNum] = i;
	}
	sp.segs[xtalNum][sp.mult[xtalNum]] = i;
	sp.mult[xtalNum]++;
      } 
    }

    sp.ccE[xtalNum] = g3Out.xtals[xtal].ccE[0]; 
  }

  return 0;
}

