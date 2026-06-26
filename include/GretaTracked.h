#ifndef __GRETA_TRACKED_H
#define __GRETA_TRACKED_H

#include "Rtypes.h"
#include "TObject.h"

#include <vector>

/* Upper bound when filling from trackEvent (clusters + pair production). */
#define GRETA_TRACK_MAX_GAMMAS 200

/**
 * One tracked gamma ray (lab frame). Units match TRACKED_GAMMA_RAY / writeTrack_addtrack:
 * positions in mm; eSum and interaction energies in keV.
 */
class GretaTrackedGamma : public TObject {
 public:
  Float_t eSum;
  Int_t ndet;
  Float_t fom;
  Short_t tracked;
  Long64_t timestamp;

  Float_t x0_mm;
  Float_t y0_mm;
  Float_t z0_mm;
  Float_t e0_keV;

  Float_t x1_mm;
  Float_t y1_mm;
  Float_t z1_mm;
  Float_t e1_keV;

  Short_t fhcrID;

  Float_t doppler;
  
 public:
  GretaTrackedGamma();
  ~GretaTrackedGamma() { ; }
  void Clear();

 private:
  ClassDef(GretaTrackedGamma, 2);
};

/**
 * Tracked gamma rays for one built coincidence event.
 */
class GretaTrackedEvent : public TObject {
 public:
  /** Return code from trackEvent (0 = success). Tracking disabled / skipped uses -1. */
  Int_t trackReturnCode;
  std::vector<GretaTrackedGamma> gammas;

 public:
  GretaTrackedEvent();
  ~GretaTrackedEvent() { ; }
  void Reset();
  Int_t nGamma() const { return static_cast<Int_t>(gammas.size()); }

 private:
  ClassDef(GretaTrackedEvent, 2);
};

#endif
