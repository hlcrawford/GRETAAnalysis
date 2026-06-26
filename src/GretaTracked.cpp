#include "GretaTracked.h"

ClassImp(GretaTrackedGamma);
ClassImp(GretaTrackedEvent);

GretaTrackedGamma::GretaTrackedGamma() { Clear(); }

void GretaTrackedGamma::Clear() {
  eSum = 0.f;
  ndet = 0;
  fom = 0.f;
  tracked = 0;
  timestamp = 0;
  x0_mm = y0_mm = z0_mm = 0.f;
  e0_keV = 0.f;
  x1_mm = y1_mm = z1_mm = 0.f;
  e1_keV = 0.f;
  fhcrID = 0;
  doppler = 0;
}

GretaTrackedEvent::GretaTrackedEvent() { Reset(); }

void GretaTrackedEvent::Reset() {
  trackReturnCode = -1;
  gammas.clear();
}
