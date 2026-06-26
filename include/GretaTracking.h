#ifndef __GRETA_TRACKING_H
#define __GRETA_TRACKING_H

#include "TString.h"

class GRETA;
class g2OUT;
class GretaTrackedEvent;
struct GretaTrackingWorkspace;

/** Initialize CTK tracking (chat file). Call once after ROOT startup. Exits on chat errors (same as trackMain). */
void GretaTrackingInit(const TString &chatFile);

/** Run clustering + tracking on the current g2Out and fill gret->tracked. */
void GretaTrackingRun(GRETA *gret);

/**
 * Thread-safe entry point: use the given workspace (set active for trackEvent / pairprod).
 * g2 must be the pre-Doppler coincidence (same as GretaTrackingRun in readGreta).
 */
void GretaTrackingRunOnG2(const g2OUT &g2, GretaTrackingWorkspace &ws, GretaTrackedEvent &out);

#endif
