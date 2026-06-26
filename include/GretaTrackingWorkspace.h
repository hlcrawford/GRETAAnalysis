#ifndef __GRETA_TRACKING_WORKSPACE_H
#define __GRETA_TRACKING_WORKSPACE_H

#include <sys/times.h>

#include "ctk.h"

/**
 * Per-event mutable state for CTK tracking (replaces file-scope npp / static buffers in GretaTracking.cpp).
 * Set the active workspace with GretaTrackingSetActiveWorkspace() before trackEvent.
 */
struct GretaTrackingWorkspace {
  struct tms timesThen;
  STAT ctkStat;
  SHELLHIT shellhit;
  CLUSTER_INTPTS clstr[MAXCLUSTERHITS];
  int nclusters;

  int nPairProd;
  float ppE[MAX_NDET];
  float ppX[MAX_NDET];
  float ppY[MAX_NDET];
  float ppZ[MAX_NDET];
  float ppFom[MAX_NDET];
  long long int ppTs[MAX_NDET];

  GEBDATA *gd;
  PAYLOAD *pl;

  GretaTrackingWorkspace();
  ~GretaTrackingWorkspace();
  GretaTrackingWorkspace(const GretaTrackingWorkspace &) = delete;
  GretaTrackingWorkspace &operator=(const GretaTrackingWorkspace &) = delete;

  void EnsureTrackBuffers();
  void ResetPairProd();  /* clears nPairProd */
};

/** Active workspace for tracking TU (thread_local; default workspace for single-threaded readGreta). */
GretaTrackingWorkspace *GretaTrackingActiveWorkspace();
GretaTrackingWorkspace *GretaTrackingSetActiveWorkspace(GretaTrackingWorkspace *ws);

/** Default workspace used by GretaTrackingInit / GretaTrackingRun(GRETA*). */
GretaTrackingWorkspace &GretaTrackingDefaultWorkspace();

#endif
