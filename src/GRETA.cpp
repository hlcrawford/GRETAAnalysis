#include "GRETA.h"

#define DEBUG 0
#define KEEP_WAVEFORMS 0

ClassImp(rotationMatrix);
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

/**************************************************************/
/* rotationMatrix class functions *****************************/
/**************************************************************/

/*! Loads the rotation matrix that maps from crystal coordinates 
    to world coordinates based on the hole position of a detector.
    This function reads the text format file, NOT the binary file.

    \param file A string value, the filename of the text format 
           rotation matrix (typically crmat.dat)
    \return Returns 0 if successful, and 1 if the file couldn't be 
            opened 
*/

Int_t rotationMatrix::ReadMatrix(TString file) {
  FILE *fp;

  Float_t f1, f2, f3, f4;
  Int_t pos, xtal;
  Int_t nn = 0;
  char *st, str[256];
  
  fp = fopen(file.Data(), "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", file.Data());
    exit(1);
  } else {
    printf("\"%s\" open...", file.Data());
  }
  
  nn = 0;
  st = fgets(str, 256, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %i", &pos, &xtal);
      for (Int_t i=0; i<4; i++) {
	st = fgets(str, 256, fp);
	sscanf(str, "%f %f %f %f", &f1, &f2, &f3, &f4);
	crmat[pos-1][xtal][i][0] = f1;
	crmat[pos-1][xtal][i][1] = f2;
	crmat[pos-1][xtal][i][2] = f3;
	crmat[pos-1][xtal][i][3] = f4;
      }
      nn++;
    }
    /* Attempt to read the next line */
    st = fgets(str, 256, fp);
  }
  
  printf("Read %i rotation matrix coefficients.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

/**************************************************************/

/*! Calculates the world position for a point in crystal coordinates
    based on the rotation matrix provided and the hole number (position)
    of the crystal

    \param crystalID The integer crystal identification - the hole 
           number comes from crystalID/4, and the crystal number from
	   crystalID%4
    \param xyz TVector3 value of the position to be transformed into 
           world coordinates space
    \return Returns the TVector3 corresponding to the point in world
            coordinate space
*/

TVector3 rotationMatrix::crys2Lab(Int_t crystalID, TVector3 xyz) {

  Int_t detectorPosition = ((crystalID)/4);
  Int_t crystalNumber = (crystalID%4);
  if (crystalNumber==0) { crystalNumber = 4; detectorPosition-= 1; }

  crystalNumber -= 1; // Crystal numbers need to be 0..3

  TVector3 xyzLab;
  xyzLab.SetX((Double_t)((crmat[detectorPosition][crystalNumber][0][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][0][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][0][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][0][3]) ));
  xyzLab.SetY((Double_t)((crmat[detectorPosition][crystalNumber][1][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][1][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][1][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][1][3]) ));
  xyzLab.SetZ((Double_t)((crmat[detectorPosition][crystalNumber][2][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][2][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][2][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][2][3]) ));
  
  return xyzLab;
}

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

void g2IntPts::Reset() {
  e = 0;
}

void g2CrystalEvent::Reset() {
  intPts.clear();
  for (Int_t i=0; i<40; i++) {
    tr[i].clear();
  }
}

TVector3 g2CrystalEvent::maxIntPtXYZLab() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intPts.size(); i++) {
    if (intPts[i].e > maxE) {
      maxE = intPts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intPts[max].xyzLab; }
  else { return TVector3(0., 0., 0.); }
}

Double_t g2CrystalEvent::maxIntPtX() {
  Float_t maxE = 0.0;
  Int_t maxI = -1;
  if (intPts.size() > 0) {
    for (Int_t i=0; i<intPts.size(); i++) {
      if (intPts[i].e > maxE) {
	maxE = intPts[i].e;
	maxI = i;
      }
    }
    return intPts[maxI].xyzCrystal.X();
  } else {
    return -10000.;
  }
}

Double_t g2CrystalEvent::maxIntPtY() {
  Float_t maxE = 0.0;
  Int_t maxI = -1;
  if (numIntPts() > 0) {
    for (Int_t i=0; i<intPts.size(); i++) {
      if (intPts[i].e > maxE) {
	maxE = intPts[i].e;
	maxI = i;
      }
    }
    return intPts[maxI].xyzCrystal.Y();
  } else {
    return -10000.;
  }
}

Double_t g2CrystalEvent::maxIntPtR() {
  if (numIntPts() > 0) {
    return TMath::Sqrt(maxIntPtY()*maxIntPtY() + maxIntPtX()*maxIntPtX());
  }
  else {
    return -1;
  }
}

Double_t g2CrystalEvent::maxIntPtZ() {
  Float_t maxE = 0.0;
  Int_t maxI = -1;
  if (numIntPts() > 0) {
    for (Int_t i=0; i<intPts.size(); i++) {
      if (intPts[i].e > maxE) {
	maxE = intPts[i].e;
	maxI = i;
      }
    }
    return intPts[maxI].xyzCrystal.Z();
  } else {
    return -10000.;
  }
}

Float_t g2CrystalEvent::gTheta() {
  if (numIntPts() > 0) {
    return maxIntPtXYZLab().Theta();
  } else { return 0.0; }
}

Float_t g2CrystalEvent::gPhi() {
  if (numIntPts() > 0) {
    if (maxIntPtXYZLab().Phi() < 0) {
      return (maxIntPtXYZLab().Phi() + TMath::TwoPi());
    } else {
      return maxIntPtXYZLab().Phi();
    }
  } else { return 0.0; }
}

UInt_t g2CrystalEvent::numIntPts() { return intPts.size(); }



void g2OUT::Reset() {
  xtals.clear();
}

UInt_t g2OUT::crystalMult() { return xtals.size(); } 

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

  ReadSegmentCenters("gretaCalibrations/segmentCenters.dat");

  if (ctrl->calibrationRun) {
    doingACalibration = 1;
  } else {
    doingACalibration = 0;
  }

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

  g2X.Clear();
  g2Out.Reset();
}

void GRETA::ReadSegmentCenters(TString filename) {
  FILE *fp;
  Int_t i1, i2, nn;  Float_t f1, f2, f3;
  char *st, str[128];
  
  /* Open file */
  fp = fopen(filename.Data(), "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", filename.Data());
    exit(1);
  }
  printf("\"%s\" open...", filename.Data());
  
  /* Read values */
  nn = 0;
  st = fgets(str, 64, fp);
  while (st!=NULL) {
    if (str[0] == 35 || str[0] == 59 || str[0] == 10) {
      /* Comment or blank line.  Do nothing. */
    } else {
      sscanf(str, "%i %i %f %f %f", &i1, &i2, &f1, &f2, &f3);
      if (i1>=0 && i1<=1 && i2>=0 && i2<=35) {
        segCenter[i1][0][i2] = f1;
        segCenter[i1][1][i2] = f2;
        segCenter[i1][2][i2] = f3;
      }
      nn++;
    }
    
    /* Attempt to read next line. */
    st = fgets(str, 64, fp);
  }
  printf("Read %i segment positions.\n", nn-1);
  
  fclose(fp);
}
  

Int_t GRETA::getMode3(FILE *inf, Int_t evtLength, Int_t subType, Int_t type) {
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
  g3X.version = wform.version;
  g3X.id = wform.id;
  g3X.trLen = ntohs(wform.tr_len);
  // g3X.trigSrc = ntohs(wform.trig_src);
  // g3X.pad = ntohs(wform.pad);
  g3X.timestamp = ntoh64((uint64_t)wform.timestamp);
  g3X.pileup = ntoh64((uint64_t)wform.pileup);
  for (Int_t i=0; i<40; i++) {
    g3X.ener[i] = ntohl(wform.ener[i]);
  }
  g3X.hist_corr[0][0] = ntohs(wform.hist_corr[0][0]);
  g3X.hist_corr[0][1] = ntohs(wform.hist_corr[0][1]);
  g3X.hist_corr[1][0] = ntohs(wform.hist_corr[1][0]);
  g3X.hist_corr[1][1] = ntohs(wform.hist_corr[1][1]);
  g3X.t0 = ntohs(wform.t0);
  g3X.subt0 = ntohs(wform.sub_t0);
  g3X.tLEDCore = ntohs(wform.t_led_core);
  g3X.tCFDCore = ntohs(wform.t_cfd_core);
  g3X.tLEDFirst = ntohs(wform.t_led_first);

  // dT[subType]->Fill(g3X.hist_corr[0][1]);
  // dT[subType]->Fill(g3X.hist_corr[1][1]);
  for (Int_t i=36; i<40; i++) {
    g3X.ccE.push_back((g3X.ener[i]/128.)*gain[subType][i] + offset[subType][i]);
  }
  for (Int_t i=0; i<36; i++) {
    g3X.segE.push_back((g3X.ener[i]/128.)*gain[subType][i] + offset[subType][i]);
  }

  if (evtLength <= 9000) {
    vector<int16_t> itr;
    if (KEEP_WAVEFORMS) {
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

  analyzeMode3(&g3X);
  g3Out.xtals.push_back(g3X);
  return 0;
  
}

Int_t GRETA::analyzeMode3(g3CrystalEvent *g3) {

  for (Int_t i=0; i<36; i++) {
    if (g3->segE[i] > 5) {
      g3->segsHit.push_back(i);
    }
  }

  g3->noNeighbors = 1;

 
  //  if (g3->segsHit.size() > 1) {
  //  for (Int_t i=0; i<g3->segsHit.size(); i++) {
  //  for (Int_t j=i+1; j<g3->segsHit.size(); j++) {
  //	if (abs(j-i) == 6) {
  //	  g3->noNeighbors = 0;
  //	}
  //	if (abs(j-i) == 1) {
  //	  if ((i%5==0 && j%6 == 0) || (i%6==0 && j%5 == 0)) {
  //	    g3->noNeighbors = 1;
  //	  } else {
  //	    g3->noNeighbors = 0;
  //	  }
  //	}
  //  }
  //  }
  // }
  
  return 0;

}

Int_t GRETA::getMode2(FILE *inf, Int_t evtLength, Int_t subType) {

  gretaIntPtMsg ipmsg;
  Int_t siz = 0;

  Int_t bytesRead = 0;
  
  siz = fread(&ipmsg, 1, sizeof(gretaIntPtMsg), inf);
  bytesRead = siz;
  // printf("bytesRead %d\n", bytesRead);
  
  if (DEBUG) {
    printf("Decomp message!\n"); 
    printf("  Version:  %i\n", ipmsg.version); 
    printf("  ID:           %i\n", ipmsg.id);
    printf("  timestamp: %lld\n", ipmsg.timestamp);
    printf("  numFits: %d\n", ipmsg.num_fits);
    printf("  event length: %d\n", evtLength);
    printf("Length remaining for intpts: %lu\n", evtLength - sizeof(gretaIntPtMsg));
  } 

  vector<intPtFit> intPtFits;
  vector<intPt> intpts;
  intPtFit ipfit;
  intPt ip;
  
  for (Int_t i=0; i<ipmsg.num_fits; i++) {
    siz = fread(&ipfit, 1, sizeof(intPtFit), inf);
    intPtFits.push_back(ipfit);
    bytesRead += siz;

    for (Int_t j=0; j<(ipfit.num_interactions); j++) {
      siz = fread(&ip, 1, sizeof(intPt), inf);
      intpts.push_back(ip);
      bytesRead += siz;
    }
  }
  
  g2X.Reset();
  g2X.id = ipmsg.id;
  // g2X.trigSrc = ntohs(ipmsg.trig_src);
  g2X.timestamp = ipmsg.timestamp;
  // g2X.pileup = ntoh64((uint64_t)ipmsg.pileup);
  for (Int_t i=0; i<40; i++) {
    g2X.ener[i] = (Float_t)ipmsg.ener[i]/1000.;
  }
  if (DEBUG) {
    printf("Energies: \n");
    for (Int_t i=0; i<40; i++) {
      cout << i << " " << g2X.ener[i] << endl;
    }
  }
  
  //  g2X.t0 = ntohs(wform.t0);
  //  g2X.subt0 = ntohs(wform.sub_t0);
  //  g2X.tLEDCore = ntohs(wform.t_led_core);
  //  g2X.tCFDCore = ntohs(wform.t_cfd_core);
  //  g2X.tLEDFirst = ntohs(wform.t_led_first);
  g2X.num_fits = ipmsg.num_fits;
  
  Int_t ipIndex = 0;
  for (Int_t i=0; i<ipmsg.num_fits; i++) {
    g2X.errorCode = ipfit.errorCode;
    for (Int_t j=0; j<ipfit.num_interactions; j++) {
      g2Xip.Clear();
      g2Xip.ir = (intpts[ipIndex].ir);
      g2Xip.ip = (intpts[ipIndex].ip);
      g2Xip.iz = (intpts[ipIndex].iz);
      g2Xip.seg = (intpts[ipIndex].seg);
      g2Xip.xyzCrystal = TVector3((intpts[ipIndex].x), (intpts[ipIndex].y), (intpts[ipIndex].z));
      g2Xip.xyzLab = rot.crys2Lab(subType, g2Xip.xyzCrystal);
      g2Xip.xyzLabSeg = rot.crys2Lab(subType, TVector3(segCenter[subType%2][0][g2Xip.seg], segCenter[subType%2][1][g2Xip.seg], segCenter[subType%2][2][g2Xip.seg]));
      g2Xip.e = (intpts[ipIndex].e);
      g2X.intPts.push_back(g2Xip);
      ipIndex++;
      if (DEBUG) {
	printf("ip%d -- %f %f %f %f\n", ipIndex-1, g2Xip.xyzCrystal.X(), g2Xip.xyzCrystal.Y(), g2Xip.xyzCrystal.Z(), g2Xip.e);
      }
    }
  }

  /* Mode 2+3 appends the waveform to the Mode2 basically */

  // printf("Length used %d\n", bytesRead);
  if (evtLength - bytesRead != 0) {
    // printf("I still have data, probably waveforms...\n");

    uint16_t wfLength[40];
    siz = fread(wfLength, 1, 2*40, inf);
    bytesRead += siz;

    /* Get total compressed length of waveforms */
    Int_t sumCompressedLengths = 0;
    uint16_t sumDecompressionErrs = 0;
    for (Int_t i=0; i<40; i++) {
      sumCompressedLengths += ntohs(wfLength[i]);
    }
    if (sumCompressedLengths*4 != evtLength - bytesRead) {
      printf("Compressed waveform lengths do not match expected, %u != %u\n", sumCompressedLengths*4, evtLength - bytesRead);
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
      sumDecompressionErrs += dptc_unpack16(temp32, ntohs(wfLength[ch]), (uint16_t*)( &trTemp ), 192, 16);
      for (Int_t l=0; l<192; l++) {
    	g2X.tr[ch].push_back(trTemp[l]);
      }
      temp32 += ntohs(wfLength[ch]);      
    }
  }
  
  g2Out.xtals.push_back(g2X);
  return 0;
  
}
