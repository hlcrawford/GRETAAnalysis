/**
 * Tracking debug helpers: optional per-event trace (GRETA_TRACK_DEBUG) and
 * always-on run summary counters (GretaTrackStats*).
 */

#include "GretaTrackingDebug.h"

#include <cstdio>

#include "GRETA.h"
#include "GretaTracked.h"

namespace {

struct RunStats {
  unsigned long long nCoincidences = 0;
  unsigned long long nG2Crystals = 0;
  unsigned long long nG2CrystalsWithIPs = 0;
  unsigned long long nCrystalsPackedForTrack = 0;
  unsigned long long nTrackedGammas = 0;
  unsigned long long nCoincidencesZeroGammas = 0;
  unsigned long long nTrackReturnNonZero = 0;
};

RunStats s_stats;

}  // namespace

void GretaTrackStatsReset() { s_stats = RunStats{}; }

void GretaTrackStatsRecordCoincidence(const g2OUT &g2Out, int nPackedForTrack,
				      int trackReturnCode, const GretaTrackedEvent &tracked)
{
  const unsigned long long nx = static_cast<unsigned long long>(g2Out.xtals.size());
  unsigned long long withIP = 0;
  for (unsigned long long i = 0; i < nx; i++) {
    if (g2Out.xtals[static_cast<size_t>(i)].intPts.size() > 0) withIP++;
  }

  s_stats.nCoincidences++;
  s_stats.nG2Crystals += nx;
  s_stats.nG2CrystalsWithIPs += withIP;
  s_stats.nCrystalsPackedForTrack += static_cast<unsigned long long>(nPackedForTrack > 0 ? nPackedForTrack : 0);
  if (trackReturnCode != 0) s_stats.nTrackReturnNonZero++;

  const Int_t ng = tracked.nGamma();
  if (ng > 0) {
    s_stats.nTrackedGammas += static_cast<unsigned long long>(ng);
  } else {
    s_stats.nCoincidencesZeroGammas++;
  }
}

void GretaTrackStatsPrintSummary(FILE *out)
{
  if (!out) out = stderr;

  const double invCoinc = (s_stats.nCoincidences > 0)
			      ? 1.0 / static_cast<double>(s_stats.nCoincidences)
			      : 0.0;
  const double gammasPerCoinc = static_cast<double>(s_stats.nTrackedGammas) * invCoinc;
  const double crystalsPerCoinc = static_cast<double>(s_stats.nG2Crystals) * invCoinc;
  const double ratioGammasToCrystals =
    (s_stats.nG2Crystals > 0)
	? static_cast<double>(s_stats.nTrackedGammas) / static_cast<double>(s_stats.nG2Crystals)
	: 0.0;
  const double ratioGammasToCrystalsWithIP =
    (s_stats.nG2CrystalsWithIPs > 0)
	? static_cast<double>(s_stats.nTrackedGammas) / static_cast<double>(s_stats.nG2CrystalsWithIPs)
	: 0.0;

  fprintf(out,
	  "\n"
	  "======== GRETA tracking run summary (coincidence aggregates) ========\n"
	  "  coincidences processed:              %llu\n"
	  "  mode-2 crystals (g2, all):         %llu  (%.2f per coincidence)\n"
	  "  mode-2 crystals with intPts:       %llu  (%.2f per coincidence)\n"
	  "  crystals packed for trackEvent:    %llu\n"
	  "  tracked gammas written (tracked branch): %llu  (%.2f per coincidence)\n"
	  "  coincidences with zero gammas:     %llu  (%.1f%%)\n"
	  "  coincidences trackReturnCode!=0:   %llu  (%.1f%%)\n"
	  "  tracked gammas / mode-2 crystals:  %.3f\n"
	  "  tracked gammas / crystals w/ IPs:  %.3f\n"
	  "====================================================================\n"
	  "  Compare spectrum totals: sum(ener36) ~ mode-2 crystals; sum(eSum) ~ tracked gammas.\n"
	  "====================================================================\n\n",
	  s_stats.nCoincidences,
	  s_stats.nG2Crystals, crystalsPerCoinc,
	  s_stats.nG2CrystalsWithIPs, static_cast<double>(s_stats.nG2CrystalsWithIPs) * invCoinc,
	  s_stats.nCrystalsPackedForTrack,
	  s_stats.nTrackedGammas, gammasPerCoinc,
	  s_stats.nCoincidencesZeroGammas,
	  100.0 * static_cast<double>(s_stats.nCoincidencesZeroGammas) * invCoinc,
	  s_stats.nTrackReturnNonZero,
	  100.0 * static_cast<double>(s_stats.nTrackReturnNonZero) * invCoinc,
	  ratioGammasToCrystals,
	  ratioGammasToCrystalsWithIP);
}

#ifdef GRETA_TRACK_DEBUG

#include "gdecomp.h"

extern TRACKINGPAR Pars;
extern int npp;
extern float ppe[MAX_NDET];
extern float ppx[MAX_NDET];
extern float ppy[MAX_NDET];
extern float ppz[MAX_NDET];
extern float ppfom[MAX_NDET];
extern long long int ppTS[MAX_NDET];

namespace {

unsigned long long s_evtLabel = 0;
unsigned s_evtSeq = 0;
bool s_verbose = false;

void banner(const char *title)
{
  fprintf(stderr,
	  "\n"
	  "================================================================================\n"
	  " GRETA TRACK DEBUG  event_seq=%u  label=%llu\n"
	  " %s\n"
	  "================================================================================\n",
	  s_evtSeq, (unsigned long long)s_evtLabel, title);
}

void rule() { fprintf(stderr, "--------------------------------------------------------------------------------\n"); }

}  // namespace

void GretaTrackDbgSetMode2CrystalMult(unsigned nMode2Crystals)
{
  s_verbose = (nMode2Crystals > 15);
}

bool GretaTrackDbgIsVerbose() { return s_verbose; }

void GretaTrackDbgSetEventLabel(unsigned long long label)
{
  s_evtLabel = label;
  s_evtSeq++;
}

void GretaTrackDbgPrintG2Input(const g2OUT &g2Out)
{
  if (!s_verbose) return;
  banner("Mode-2 input (g2OUT) — per crystal before packing");
  const UInt_t nx = (UInt_t)g2Out.xtals.size();
  fprintf(stderr, "  crystals in coincidence: %u\n", nx);
  for (UInt_t i = 0; i < nx; i++) {
    const g2CrystalEvent &cx = g2Out.xtals[i];
    fprintf(stderr,
	    "  [%u] id=%u rhSubType=%u ts=%lld num_fits=%d nIntPts=%u\n"
	    "       ener36=%.3f ener37=%.3f ener38=%.3f (keV)  t0=%d t_cfd_core=%d\n",
	    i, (unsigned)cx.id, (unsigned)cx.rhSubType, (long long)cx.timestamp,
	    (int)cx.num_fits, (unsigned)cx.intPts.size(),
	    cx.ener[36], cx.ener[37], cx.ener[38],
	    (int)cx.t0, (int)cx.t_cfd_core);
    for (UInt_t j = 0; j < cx.intPts.size(); j++) {
      const g2IntPts &ip = cx.intPts[j];
      fprintf(stderr,
	      "       intPt[%u] seg=%d e=%.4f MeV  lab(cm)=(%.2f, %.2f, %.2f)  crystal=(%.2f, %.2f, %.2f)\n",
	      j, (int)ip.seg, ip.e,
	      ip.xyzLab.X(), ip.xyzLab.Y(), ip.xyzLab.Z(),
	      ip.xyzCrystal.X(), ip.xyzCrystal.Y(), ip.xyzCrystal.Z());
    }
  }
  rule();
}

void GretaTrackDbgPrintPackedInput(int nCrystals, GEBDATA *gd, PAYLOAD *pl)
{
  if (!s_verbose) return;
  banner("Packed TRACK_STRUCT input (sent to trackEvent)");
  fprintf(stderr, "  crystals with hits sent to trackEvent: %d\n", nCrystals);
  fprintf(stderr, "  target (mm, from Pars): (%.2f, %.2f, %.2f)  [negated in GretaTrackingRun]\n",
	  -Pars.target_x, -Pars.target_y, -Pars.target_z);
  for (int i = 0; i < nCrystals; i++) {
    const CRYS_INTPTS *cr = (const CRYS_INTPTS *)(void *)&pl[i].p[0];
    fprintf(stderr,
	    "  [%d] geb type=%d crystal_id=%d ts=%lld tot_e=%.2f keV num=%d\n",
	    i, gd[i].type, cr->crystal_id, (long long)gd[i].timestamp,
	    cr->tot_e, cr->num);
    for (int j = 0; j < cr->num; j++) {
      fprintf(stderr,
	      "       ip[%d] seg=%d e=%.2f keV  lab(cm)=(%.2f, %.2f, %.2f)  seg_ener=%.2f keV\n",
	      j, cr->intpts[j].seg, cr->intpts[j].e,
	      cr->intpts[j].x, cr->intpts[j].y, cr->intpts[j].z,
	      cr->intpts[j].seg_ener);
    }
  }
  rule();
}

void GretaTrackDbgPrintShellHits(const char *stage, const SHELLHIT *shellhit)
{
  if (!s_verbose || !shellhit) return;
  fprintf(stderr, "\n--- Shell hits / clustering [%s] ---\n", stage);
  fprintf(stderr, "  nhit=%d  NumClusters=%d\n", shellhit->nhit, shellhit->NumClusters);
  for (int j = 0; j < shellhit->nhit; j++) {
    fprintf(stderr,
	    "  hit[%3d] crystal_id=%d detno=%d cluster=%d  pos(cm)=(%.2f, %.2f, %.2f)  edet=%.4f MeV  esum=%.4f MeV  ts=%lld\n",
	    j, shellhit->crystal_id[j], shellhit->detno[j], shellhit->ClusterNumber[j],
	    shellhit->XX[j], shellhit->YY[j], shellhit->ZZ[j],
	    shellhit->edet[j], shellhit->esum[j], (long long)shellhit->timestamp[j]);
  }
  rule();
}

void GretaTrackDbgPrintClusters(const char *stage, int nClusters, CLUSTER_INTPTS clstr[])
{
  if (!s_verbose) return;
  fprintf(stderr, "\n--- Clusters [%s]  nClusters=%d ---\n", stage, nClusters);
  for (int ic = 0; ic < nClusters; ic++) {
    const CLUSTER_INTPTS &c = clstr[ic];
    fprintf(stderr,
	    "  cluster[%2d] valid=%d tracked=%d ndet=%d esum=%.4f MeV fom=%.4f trackno=%d bestPerm=%d\n",
	    ic, c.valid, c.tracked, c.ndet, c.esum, c.fom, c.trackno, c.bestPermutation);
    for (int j = 0; j < c.ndet; j++) {
      const CL_INTPTS &h = c.intpts[j];
      fprintf(stderr,
	      "       hit[%d] order=%d detno=%d  pos(cm)=(%.2f, %.2f, %.2f)  edet=%.4f MeV  ts=%lld\n",
	      j, h.order, h.detno, h.xx, h.yy, h.zz, h.edet, (long long)h.timestamp);
    }
  }
  rule();
}

void GretaTrackDbgPrintPairProduction(int npp)
{
  if (!s_verbose || npp <= 0) return;
  fprintf(stderr, "\n--- Pair-production candidates (npp=%d) ---\n", npp);
  for (int i = 0; i < npp; i++) {
    fprintf(stderr,
	    "  pp[%d] e=%.4f MeV  pos(cm)=(%.2f, %.2f, %.2f)  fom=%.4f  ts=%lld\n",
	    i, ppe[i], ppx[i], ppy[i], ppz[i], ppfom[i], (long long)ppTS[i]);
  }
  rule();
}

void GretaTrackDbgPrintTrackedOutput(int trackStatus, const GretaTrackedEvent &tracked)
{
  if (!s_verbose) return;
  banner("Final tracked output (GretaTrackedEvent)");
  fprintf(stderr, "  trackReturnCode=%d  nGamma=%d\n", trackStatus, tracked.nGamma());
  for (Int_t i = 0; i < tracked.nGamma(); i++) {
    const GretaTrackedGamma &g = tracked.gammas[static_cast<size_t>(i)];
    fprintf(stderr,
	    "  gamma[%d] eSum=%.4f keV  fom=%.4f  ndet=%d  tracked=%d  fhcrID=%d  ts=%lld\n"
	    "           x0(mm)=(%.2f, %.2f, %.2f)  e0=%.2f keV\n"
	    "           x1(mm)=(%.2f, %.2f, %.2f)  e1=%.2f keV  doppler=%.5f\n",
	    (int)i, g.eSum, g.fom, g.ndet, (int)g.tracked, (int)g.fhcrID, (long long)g.timestamp,
	    g.x0_mm, g.y0_mm, g.z0_mm, g.e0_keV,
	    g.x1_mm, g.y1_mm, g.z1_mm, g.e1_keV, g.doppler);
  }
  fprintf(stderr, "================================================================================\n\n");
}

#endif  /* GRETA_TRACK_DEBUG */
