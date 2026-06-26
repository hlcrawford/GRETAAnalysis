#ifndef __GRETA_TRACKING_DEBUG_H
#define __GRETA_TRACKING_DEBUG_H

/**
 * Tracking debug prints for readGreta / GretaTracking.
 *
 * Enable at build time:
 *   scons TRACK_DEBUG=1 readGreta
 *
 * When disabled, all GretaTrackDbg* calls compile to no-ops.
 *
 * Verbose per-event trace prints only when mode-2 has > 10 crystals
 * (see GretaTrackDbgSetMode2CrystalMult).
 */

#ifdef GRETA_TRACK_DEBUG

#include "GretaTracked.h"
#include "ctk.h"

class g2OUT;
class GretaTrackedEvent;

/** Enable verbose stderr trace for this event if nMode2Crystals > 10. */
void GretaTrackDbgSetMode2CrystalMult(unsigned nMode2Crystals);
bool GretaTrackDbgIsVerbose();

void GretaTrackDbgSetEventLabel(unsigned long long label);
void GretaTrackDbgPrintG2Input(const g2OUT &g2Out);
void GretaTrackDbgPrintPackedInput(int nCrystals, GEBDATA *gd, PAYLOAD *pl);
void GretaTrackDbgPrintShellHits(const char *stage, const SHELLHIT *shellhit);
void GretaTrackDbgPrintClusters(const char *stage, int nClusters, CLUSTER_INTPTS clstr[]);
void GretaTrackDbgPrintPairProduction(int npp);
void GretaTrackDbgPrintTrackedOutput(int trackStatus, const GretaTrackedEvent &tracked);

#else

#define GretaTrackDbgSetMode2CrystalMult(n) ((void)(n))
#define GretaTrackDbgIsVerbose() (false)
#define GretaTrackDbgSetEventLabel(label) ((void)(label))
#define GretaTrackDbgPrintG2Input(g2Out) ((void)(g2Out))
#define GretaTrackDbgPrintPackedInput(n, gd, pl) ((void)(n), (void)(gd), (void)(pl))
#define GretaTrackDbgPrintShellHits(stage, sh) ((void)(stage), (void)(sh))
#define GretaTrackDbgPrintClusters(stage, n, cl) ((void)(stage), (void)(n), (void)(cl))
#define GretaTrackDbgPrintPairProduction(n) ((void)(n))
#define GretaTrackDbgPrintTrackedOutput(st, tr) ((void)(st), (void)(tr))

#endif

/* Run-level counters (always compiled; summary printed when tracking is enabled). */
class g2OUT;
class GretaTrackedEvent;

void GretaTrackStatsReset();
void GretaTrackStatsRecordCoincidence(const g2OUT &g2Out, int nPackedForTrack,
				      int trackReturnCode, const GretaTrackedEvent &tracked);
void GretaTrackStatsPrintSummary(FILE *out = stderr);

#endif
