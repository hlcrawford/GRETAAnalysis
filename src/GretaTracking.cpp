/**
 * Bridge GRETA g2OUT → trackMain TRACK_STRUCT / trackEvent, then fill GretaTrackedEvent.
 */

#include <sys/times.h>

#include <algorithm>
#include <cstdlib>

#include "GRETA.h"
#include "GretaTracking.h"
#include "GretaTracked.h"
#include "GretaTrackingDebug.h"
#include "GretaTrackingWorkspace.h"

#include "ctk.h"
#include "gdecomp.h"

extern TRACKINGPAR Pars;

namespace {

bool g_tracking_ready = false;

/** Pack one crystal into CRYS_INTPTS (lab coordinates from g2, cm; energies keV). Rotation in CTK is disabled — see GretaTrackingInit. */
void packCrystal(const g2CrystalEvent &cx, CRYS_INTPTS *cr, GEBDATA *gd) {
  std::fill_n(reinterpret_cast<unsigned char *>(cr), sizeof(CRYS_INTPTS), 0);
  gd->type = GEB_TYPE_DECOMP;
  gd->timestamp = (long long) cx.timestamp;
  gd->length = (int) sizeof(CRYS_INTPTS);

  cr->crystal_id = (int)(unsigned char) cx.id;
  cr->timestamp = (long long) cx.timestamp;
  int ni = (int) cx.intPts.size();
  if (ni > MAX_INTPTS) ni = MAX_INTPTS;
  cr->num = ni;

  cr->tot_e = cx.ener[36] * 1000.f;

  cr->t0 = (float) cx.t0;
  cr->cfd = (float) cx.t_cfd_core;
  if (ni > 0) {
    cr->chisq = cx.intPts[0].chisq;
    cr->norm_chisq = cx.intPts[0].chisq;
  }
  cr->pad = 0;

  for (int j = 0; j < ni; j++) {
    const g2IntPts &ip = cx.intPts[j];
    cr->intpts[j].x = ip.xyzLab.X();
    cr->intpts[j].y = ip.xyzLab.Y();
    cr->intpts[j].z = ip.xyzLab.Z();
    cr->intpts[j].e = ip.e * 1000.f;
    cr->intpts[j].seg = ip.seg;
    int sg = (int) ip.seg;
    if (sg >= 0 && sg < NUM_CHAN)
      cr->intpts[j].seg_ener = cx.ener[sg] * 1000.f;
    else
      cr->intpts[j].seg_ener = 0.f;
  }
}

void fillTrackedOutput(GretaTrackedEvent &out, int trackStatus, CLUSTER_INTPTS Clstr[MAXCLUSTERHITS],
                      int *nClusters, GretaTrackingWorkspace &ws) {
  out.Reset();
  out.trackReturnCode = trackStatus;
  if (trackStatus != 0) return;

  int wo[MAXCLUSTERHITS];
  long long mTS[MAXCLUSTERHITS];
  int ng = 0;
  int nCl = *nClusters;
  if (nCl > MAXCLUSTERHITS) nCl = MAXCLUSTERHITS;

  for (int iCluster = 0; iCluster < nCl; iCluster++) {
    wo[iCluster] = 0;
    mTS[iCluster] = (long long int) 0;
    int nmTS = 0;

    if (Clstr[iCluster].valid) {
      if (Clstr[iCluster].tracked && Clstr[iCluster].ndet > 0 && Clstr[iCluster].ndet < MAX_NDET) {
        mTS[iCluster] = 0;
        for (int j = 0; j < Clstr[iCluster].ndet; j++) {
          mTS[iCluster] += Clstr[iCluster].intpts[j].timestamp;
          nmTS++;
        }
        if (nmTS > 1) mTS[iCluster] /= (long long int) nmTS;

        if (nmTS > 0 && mTS[iCluster] > (long long int) 0 && Clstr[iCluster].valid == 1 &&
            Clstr[iCluster].tracked == 1) {
          wo[iCluster] = 1;
          ng++;
        }
      }
    }
  }

  if (ng == 0 && ws.nPairProd <= 0) return;

  out.gammas.reserve(static_cast<size_t>(ng + ws.nPairProd));

  for (int iCluster = 0; iCluster < nCl; iCluster++) {
    if (!wo[iCluster]) continue;
    if (static_cast<int>(out.gammas.size()) >= GRETA_TRACK_MAX_GAMMAS) break;

    GretaTrackedGamma g;
    const Float_t e_eV = Clstr[iCluster].esum * 1000.f;
    g.eSum = e_eV / 1000.f;
    g.ndet = Clstr[iCluster].ndet;
    g.fom = Clstr[iCluster].fom;
    g.tracked = (Short_t) Clstr[iCluster].tracked;

    bool have0 = false, have1 = false;
    for (int j = 0; j < Clstr[iCluster].ndet; j++) {
      if (Clstr[iCluster].intpts[j].order == 0) {
        g.x0_mm = Clstr[iCluster].intpts[j].xx * 10.f;
        g.y0_mm = Clstr[iCluster].intpts[j].yy * 10.f;
        g.z0_mm = Clstr[iCluster].intpts[j].zz * 10.f;
        g.e0_keV = Clstr[iCluster].intpts[j].edet * 1000.f;
        g.timestamp = Clstr[iCluster].intpts[j].timestamp;
        g.fhcrID = (Short_t) Clstr[iCluster].intpts[j].detno;
        have0 = true;
      }
      if (Clstr[iCluster].intpts[j].order == 1) {
        g.x1_mm = Clstr[iCluster].intpts[j].xx * 10.f;
        g.y1_mm = Clstr[iCluster].intpts[j].yy * 10.f;
        g.z1_mm = Clstr[iCluster].intpts[j].zz * 10.f;
        g.e1_keV = Clstr[iCluster].intpts[j].edet * 1000.f;
        have1 = true;
      }
    }
    if (!have0) {
      g.timestamp = mTS[iCluster];
    }
    if (!have1) {
      g.x1_mm = g.y1_mm = g.z1_mm = g.e1_keV = 0.f;
    }
    out.gammas.push_back(g);
  }

  for (int i = 0; i < ws.nPairProd; i++) {
    if (static_cast<int>(out.gammas.size()) >= GRETA_TRACK_MAX_GAMMAS) break;

    GretaTrackedGamma g;
    const Float_t e_eV = ws.ppE[i] * 1000.f;
    g.eSum = e_eV / 1000.f;
    g.ndet = 1;
    g.fom = ws.ppFom[i];
    g.tracked = 1;
    g.x0_mm = ws.ppX[i] * 10.f;
    g.y0_mm = ws.ppY[i] * 10.f;
    g.z0_mm = ws.ppZ[i] * 10.f;
    g.e0_keV = ws.ppE[i] * 1000.f;
    g.x1_mm = g.y1_mm = g.z1_mm = g.e1_keV = 0.f;
    g.timestamp = ws.ppTs[i];
    g.fhcrID = 0;
    out.gammas.push_back(g);
  }
}

}  // namespace

void GretaTrackingInit(const TString &chatFile) {
  if (chatFile.IsNull() || chatFile.Length() == 0) return;
  GretaTrackingWorkspace &ws = GretaTrackingDefaultWorkspace();
  setupTrack(&ws.timesThen, &ws.ctkStat, &ws.shellhit);
  readChatFile(const_cast<char *>(chatFile.Data()));
  Pars.nocrystaltoworldrot = 1;
  g_tracking_ready = true;
}

void GretaTrackingRunOnG2(const g2OUT &g2, GretaTrackingWorkspace &ws, GretaTrackedEvent &out) {
  out.Reset();

  if (!g_tracking_ready) {
    out.trackReturnCode = -1;
    return;
  }

  GretaTrackingWorkspace *prev = GretaTrackingSetActiveWorkspace(&ws);

  const UInt_t nx = static_cast<UInt_t>(g2.xtals.size());

  unsigned long long evtLabel = 0;
  if (nx > 0) evtLabel = (unsigned long long) g2.xtals[0].timestamp;
  GretaTrackDbgSetEventLabel(evtLabel);
  GretaTrackDbgPrintG2Input(g2);

  ws.EnsureTrackBuffers();
  TRACK_STRUCT trackStruct;
  trackStruct.gd = ws.gd;
  trackStruct.payload = ws.pl;

  int k = 0;
  for (UInt_t i = 0; i < nx; i++) {
    const g2CrystalEvent &cx = g2.xtals[i];
    if (cx.intPts.size() == 0) continue;
    CRYS_INTPTS *cr = (CRYS_INTPTS *) (void *) &ws.pl[k].p[0];
    packCrystal(cx, cr, &ws.gd[k]);
    k++;
  }

  if (k == 0) {
    out.trackReturnCode = -1;
    GretaTrackingSetActiveWorkspace(prev);
    return;
  }

  trackStruct.n = k;

  GretaTrackDbgPrintPackedInput(k, ws.gd, ws.pl);

  float target_pos[3];
  target_pos[0] = -Pars.target_x;
  target_pos[1] = -Pars.target_y;
  target_pos[2] = -Pars.target_z;

  int st = trackEvent(target_pos, &trackStruct, &ws.ctkStat, &ws.shellhit, ws.clstr, &ws.nclusters);

  GretaTrackDbgPrintClusters("after trackEvent", ws.nclusters, ws.clstr);
  GretaTrackDbgPrintPairProduction(ws.nPairProd);

  fillTrackedOutput(out, st, ws.clstr, &ws.nclusters, ws);
  GretaTrackDbgPrintTrackedOutput(st, out);

  GretaTrackingSetActiveWorkspace(prev);
}

void GretaTrackingRun(GRETA *gret) {
  GretaTrackingWorkspace &ws = GretaTrackingDefaultWorkspace();
  GretaTrackingRunOnG2(gret->g2Out, ws, gret->tracked);
}
