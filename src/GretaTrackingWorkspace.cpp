#include "GretaTrackingWorkspace.h"

#include <cstdlib>
#include <cstring>

namespace {

GretaTrackingWorkspace g_default_workspace;
thread_local GretaTrackingWorkspace *g_active_workspace = &g_default_workspace;

}  // namespace

GretaTrackingWorkspace::GretaTrackingWorkspace()
    : nclusters(0), nPairProd(0), gd(nullptr), pl(nullptr) {
  std::memset(&timesThen, 0, sizeof(timesThen));
  std::memset(&ctkStat, 0, sizeof(ctkStat));
  std::memset(&shellhit, 0, sizeof(shellhit));
  std::memset(clstr, 0, sizeof(clstr));
  std::memset(ppE, 0, sizeof(ppE));
  std::memset(ppX, 0, sizeof(ppX));
  std::memset(ppY, 0, sizeof(ppY));
  std::memset(ppZ, 0, sizeof(ppZ));
  std::memset(ppFom, 0, sizeof(ppFom));
  std::memset(ppTs, 0, sizeof(ppTs));
}

GretaTrackingWorkspace::~GretaTrackingWorkspace() {
  if (gd) {
    std::free(gd);
    gd = nullptr;
  }
  if (pl) {
    std::free(pl);
    pl = nullptr;
  }
}

void GretaTrackingWorkspace::EnsureTrackBuffers() {
  if (!gd) {
    gd = (GEBDATA *) std::calloc((size_t) MAXTRACK, sizeof(GEBDATA));
    pl = (PAYLOAD *) std::calloc((size_t) MAXTRACK, sizeof(PAYLOAD));
  }
}

void GretaTrackingWorkspace::ResetPairProd() { nPairProd = 0; }

GretaTrackingWorkspace *GretaTrackingActiveWorkspace() { return g_active_workspace; }

GretaTrackingWorkspace *GretaTrackingSetActiveWorkspace(GretaTrackingWorkspace *ws) {
  GretaTrackingWorkspace *prev = g_active_workspace;
  g_active_workspace = ws ? ws : &g_default_workspace;
  return prev;
}

GretaTrackingWorkspace &GretaTrackingDefaultWorkspace() { return g_default_workspace; }
