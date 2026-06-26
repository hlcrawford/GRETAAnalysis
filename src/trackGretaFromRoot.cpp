/**
 * Pass 2: read mode-2 coincidences (g2Out branch) from gdata written by readGreta.
 *
 * Example:
 *   readGreta -f data -rootFile out.root
 *   trackGreta -i out.root -o tracked.root -trackingChat tracking/greta_readGreta.chat
 */

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "TClass.h"
#include "TFile.h"
#include "TROOT.h"
#include "TVirtualStreamerInfo.h"
#include "TSystem.h"
#include "TTree.h"

#include "GRETA.h"
#include "GretaTracking.h"
#include "GretaTracked.h"
#include "GretaTrackingWorkspace.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern void TriggerDictionaryInitialization_GRETADict();

namespace {

struct TrackGretaRow {
  g2OUT g2;
  Long64_t coincTS;
  GretaTrackedEvent tracked;
  bool haveG2;
};

bool GretaDictionaryReady()
{
  return TClass::GetClass("g2OUT") && TClass::GetClass("GretaTrackedEvent");
}

void EnsureGretaDictionaryLoaded()
{
  if (GretaDictionaryReady()) return;

  const char *pwd = gSystem->WorkingDirectory();
  const TString libDir = TString::Format("%s/lib", pwd);
  gSystem->AddDynamicPath(libDir.Data());

  TString inc = gSystem->Getenv("ROOT_INCLUDE_PATH");
  if (inc.Length()) inc += ":";
  inc += libDir;
  gSystem->Setenv("ROOT_INCLUDE_PATH", inc.Data());

  const TString libPaths[] = {
    libDir + "/libGRETA",
    TString("lib/libGRETA"),
    TString("libGRETA"),
  };
  for (const TString &lp : libPaths) {
    if (gSystem->Load(lp.Data()) >= 0) break;
  }

  TriggerDictionaryInitialization_GRETADict();

  if (!GretaDictionaryReady()) {
    fprintf(stderr,
	    "trackGreta: GRETA ROOT dictionary not loaded (need g2OUT, GretaTrackedEvent).\n"
	    "  Run from TrackBackToMe/GRETAAnalysis after:  scons trackGreta\n");
    exit(1);
  }
}

/** Force TStreamerInfo to compile before parallel TTree reads (avoids BranchElement errors with -j>1). */
void CompileClassStreamer(const char *className)
{
  TClass *cl = TClass::GetClass(className, true, true);
  if (!cl) {
    fprintf(stderr, "trackGreta: warning: dictionary missing class %s\n", className);
    return;
  }
  TVirtualStreamerInfo *si = cl->GetStreamerInfo();
  if (si) si->Compile();
}

void CompileGretaStreamers()
{
  static const char *const kClasses[] = {
    "TVector3",
    "vector<TVector3>",
    "g2IntPts",
    "vector<g2IntPts>",
    "g2CrystalEvent",
    "vector<g2CrystalEvent>",
    "g2OUT",
    "GretaTrackedGamma",
    "vector<GretaTrackedGamma>",
    "GretaTrackedEvent",
    nullptr,
  };
  for (int i = 0; kClasses[i]; i++) CompileClassStreamer(kClasses[i]);
}

void BindInputBranches(TTree *tin, g2OUT *&g2, Long64_t &coincTS, bool &haveCoincTS);

/** One real read on the main thread so branch streamers match the on-disk layout. */
void WarmupInputTree(TTree *tin, Long64_t entry)
{
  if (!tin || entry < 0 || entry >= tin->GetEntries()) return;

  g2OUT *g2 = nullptr;
  Long64_t coincTS = 0;
  bool haveCoincTS = false;
  BindInputBranches(tin, g2, coincTS, haveCoincTS);
  tin->GetEntry(entry);
  if (g2) (void) g2->xtals.size();
}

void PrepareRootForParallel(const TString &inFile, const TString &treeName, Long64_t warmupEntry)
{
  ROOT::EnableThreadSafety();
  CompileGretaStreamers();

  TFile fin(inFile.Data(), "READ");
  if (fin.IsZombie()) return;
  TTree *tin = dynamic_cast<TTree *>(fin.Get(treeName.Data()));
  if (tin) WarmupInputTree(tin, warmupEntry);
  fin.Close();
}

void PrintUsage(const char *prog)
{
  fprintf(stderr,
	  "Usage: %s -i <mode2.root> -o <tracked.root> -trackingChat <chatfile>\n"
	  "       [-tree gdata] [-entryMin N] [-entryMax M] [-j N]\n"
	  "\n"
	  "  -j N   OpenMP worker threads (default 1). Prefer scripts/trackGretaParallel.sh for -j>1.\n"
	  "  -queryEntries  print TTree entry count on stdout and exit (needs -i [-tree])\n"
	  "  Input must contain g2Out (g2OUT) on gdata from readGreta\n",
	  prog);
}

int ResolveThreadCount(int jFlag)
{
  if (jFlag == 0) {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  }
  if (jFlag < 0) return 1;
  return jFlag;
}

bool OpenTree(TFile &fin, const char *treeName, TTree *&tin, bool needCoincTS)
{
  tin = dynamic_cast<TTree *>(fin.Get(treeName));
  if (!tin) {
    fprintf(stderr, "trackGreta: tree \"%s\" not found\n", treeName);
    return false;
  }
  if (!tin->GetBranch("g2Out") && !tin->GetBranch("g2")) {
    fprintf(stderr, "trackGreta: branch \"g2Out\" (g2OUT) missing on tree %s\n", treeName);
    return false;
  }
  (void) needCoincTS;
  return true;
}

const char *G2InputBranchName(TTree *tin)
{
  if (tin->GetBranch("g2Out")) return "g2Out";
  if (tin->GetBranch("g2")) return "g2";
  return nullptr;
}

void BindInputBranches(TTree *tin, g2OUT *&g2, Long64_t &coincTS, bool &haveCoincTS)
{
  g2 = nullptr;
  coincTS = 0;
  haveCoincTS = (tin->GetBranch("coincTS") != nullptr);
  const char *g2Branch = G2InputBranchName(tin);
  tin->SetBranchAddress(g2Branch, &g2);
  if (haveCoincTS) tin->SetBranchAddress("coincTS", &coincTS);
}

/**
 * Read one entry and run tracking.
 * parallelIO: serialize TTree::GetEntry (ROOT streamers); tracking still runs in parallel.
 */
bool ProcessOneEntry(TTree *tin, Long64_t entry, g2OUT *&g2, Long64_t &coincTS, bool haveCoincTS,
                     GretaTrackingWorkspace &ws, TrackGretaRow &row, bool parallelIO)
{
  auto readEntry = [&]() {
    tin->GetEntry(entry);
    row.coincTS = haveCoincTS ? coincTS : 0;
    row.haveG2 = (g2 != nullptr);
    if (row.haveG2) row.g2 = *g2;
  };

  if (parallelIO) {
#pragma omp critical(greta_tree_io)
    readEntry();
  } else {
    readEntry();
  }

  if (!row.haveG2) {
    row.tracked.Reset();
    row.tracked.trackReturnCode = -1;
    return false;
  }
  GretaTrackingRunOnG2(row.g2, ws, row.tracked);
  return true;
}

void WriteOutput(TTree *tout, const std::vector<TrackGretaRow> &rows)
{
  CompileGretaStreamers();

  g2OUT *g2Ptr = nullptr;
  Long64_t coincTS = 0;
  GretaTrackedEvent trackedObj;
  GretaTrackedEvent *tracked = &trackedObj;

  tout->Branch("g2Out", "g2OUT", &g2Ptr);
  tout->Branch("coincTS", &coincTS);
  tout->Branch("tracked", "GretaTrackedEvent", &tracked);

  for (const TrackGretaRow &row : rows) {
    if (!row.haveG2) continue;
    trackedObj = row.tracked;
    coincTS = row.coincTS;
    g2OUT g2ForBranch = row.g2;
    g2Ptr = &g2ForBranch;
    tout->Fill();
    g2Ptr = nullptr;
  }
}

int RunSerial(const TString &inFile, const TString &outFile, const TString &treeName, Long64_t iStart,
              Long64_t iStop)
{
  TFile fin(inFile.Data(), "READ");
  if (fin.IsZombie()) {
    fprintf(stderr, "trackGreta: cannot open input %s\n", inFile.Data());
    return 1;
  }

  TTree *tin = nullptr;
  if (!OpenTree(fin, treeName.Data(), tin, true)) return 1;

  g2OUT *g2 = nullptr;
  Long64_t coincTS = 0;
  bool haveCoincTS = false;
  BindInputBranches(tin, g2, coincTS, haveCoincTS);

  const Long64_t nProc = iStop - iStart;
  std::vector<TrackGretaRow> rows(static_cast<size_t>(nProc));

  GretaTrackingWorkspace ws;
  for (Long64_t i = 0; i < nProc; i++) {
    ProcessOneEntry(tin, iStart + i, g2, coincTS, haveCoincTS, ws, rows[static_cast<size_t>(i)], false);

    if (nProc > 20 && (i + 1) % (nProc / 20 + 1) == 0) {
      const int pct = (int) (100.0 * (i + 1) / nProc);
      fprintf(stderr, "\r  %d%%", pct);
      fflush(stderr);
    }
  }
  if (nProc > 20) fprintf(stderr, "\n");

  fin.Close();

  TFile fout(outFile.Data(), "RECREATE");
  TTree *tout = new TTree(treeName.Data(), "CTK tracking from mode-2 g2OUT");
  WriteOutput(tout, rows);
  fout.cd();
  fout.Write();
  fout.Close();

  fprintf(stderr, "trackGreta: wrote %lld entries to %s\n", (long long) rows.size(), outFile.Data());
  return 0;
}

int RunParallel(const TString &inFile, const TString &outFile, const TString &treeName, Long64_t iStart,
                Long64_t iStop, int nThreads)
{
  const Long64_t nProc = iStop - iStart;
  std::vector<TrackGretaRow> rows(static_cast<size_t>(nProc));

  fprintf(stderr,
	  "trackGreta: OpenMP parallel mode, %d threads (I/O serialized, tracking parallel)\n", nThreads);

  TFile fin(inFile.Data(), "READ");
  if (fin.IsZombie()) {
    fprintf(stderr, "trackGreta: cannot open input %s\n", inFile.Data());
    return 1;
  }
  TTree *tin = nullptr;
  if (!OpenTree(fin, treeName.Data(), tin, true)) return 1;

  g2OUT *g2 = nullptr;
  Long64_t coincTS = 0;
  bool haveCoincTS = false;
  BindInputBranches(tin, g2, coincTS, haveCoincTS);

  std::atomic<Long64_t> entriesDone{0};
  std::atomic<int> nullG2Count{0};

#ifdef _OPENMP
  omp_set_num_threads(nThreads);

#pragma omp parallel
  {
    GretaTrackingWorkspace ws;
    GretaTrackingSetActiveWorkspace(&ws);

#pragma omp for schedule(dynamic, 32)
    for (Long64_t i = 0; i < nProc; i++) {
      TrackGretaRow row;
      ProcessOneEntry(tin, iStart + i, g2, coincTS, haveCoincTS, ws, row, true);
      if (!row.haveG2) nullG2Count++;
      rows[static_cast<size_t>(i)] = row;

      const Long64_t done = ++entriesDone;
      if (nProc > 20 && done % (nProc / 20 + 1) == 0) {
#pragma omp critical(trackgreta_progress)
        {
          const int pct = (int) (100.0 * done / nProc);
          fprintf(stderr, "\r  %d%%", pct);
          fflush(stderr);
        }
      }
    }
  }
  fin.Close();
#else
  (void) tin;
  (void) g2;
  (void) coincTS;
  (void) haveCoincTS;
  (void) nThreads;
  fprintf(stderr, "trackGreta: OpenMP not available at build time; use -j 1\n");
  return 1;
#endif

  if (nProc > 20) fprintf(stderr, "\n");
  if (nullG2Count.load() > 0)
    fprintf(stderr, "trackGreta: warning: %d entries had null g2 branch\n", nullG2Count.load());

  TFile fout(outFile.Data(), "RECREATE");
  TTree *tout = new TTree(treeName.Data(), "CTK tracking from mode-2 g2OUT");
  WriteOutput(tout, rows);
  fout.cd();
  fout.Write();
  fout.Close();

  fprintf(stderr, "trackGreta: wrote %lld entries to %s\n", (long long) rows.size(), outFile.Data());
  return 0;
}

}  // namespace

int main(int argc, char *argv[])
{
  TString inFile;
  TString outFile;
  TString chatFile;
  TString treeName = "gdata";
  Long64_t entryMin = 0;
  Long64_t entryMax = -1;
  int jFlag = 1;
  bool queryEntries = false;

  for (int i = 1; i < argc; i++) {
    TString flag(argv[i]);
    if (flag == "-i" && i + 1 < argc) {
      inFile = argv[++i];
    } else if (flag == "-o" && i + 1 < argc) {
      outFile = argv[++i];
    } else if (flag == "-trackingChat" && i + 1 < argc) {
      chatFile = argv[++i];
    } else if (flag == "-tree" && i + 1 < argc) {
      treeName = argv[++i];
    } else if (flag == "-entryMin" && i + 1 < argc) {
      entryMin = atoll(argv[++i]);
    } else if (flag == "-entryMax" && i + 1 < argc) {
      entryMax = atoll(argv[++i]);
    } else if (flag == "-j" && i + 1 < argc) {
      jFlag = atoi(argv[++i]);
    } else if (flag == "-queryEntries") {
      queryEntries = true;
    } else if (flag == "-h" || flag == "--help") {
      PrintUsage(argv[0]);
      return 0;
    } else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      PrintUsage(argv[0]);
      return 1;
    }
  }

  if (inFile.IsNull()) {
    PrintUsage(argv[0]);
    return 1;
  }

  if (queryEntries) {
    EnsureGretaDictionaryLoaded();
    TFile fin(inFile.Data(), "READ");
    if (fin.IsZombie()) {
      fprintf(stderr, "trackGreta: cannot open input %s\n", inFile.Data());
      return 1;
    }
    TTree *tin = dynamic_cast<TTree *>(fin.Get(treeName.Data()));
    if (!tin) {
      fprintf(stderr, "trackGreta: tree \"%s\" not found in %s\n", treeName.Data(), inFile.Data());
      return 1;
    }
    printf("%lld\n", (long long) tin->GetEntries());
    return 0;
  }

  if (outFile.IsNull() || chatFile.IsNull()) {
    PrintUsage(argv[0]);
    return 1;
  }

  int nThreads = ResolveThreadCount(jFlag);
#ifndef _OPENMP
  if (jFlag != 1 && jFlag != 0) {
    fprintf(stderr, "trackGreta: warning: built without OpenMP; ignoring -j %d (using 1 thread)\n", jFlag);
    nThreads = 1;
  }
#endif

  EnsureGretaDictionaryLoaded();
  GretaTrackingInit(chatFile);

  TFile finProbe(inFile.Data(), "READ");
  if (finProbe.IsZombie()) {
    fprintf(stderr, "trackGreta: cannot open input %s\n", inFile.Data());
    return 1;
  }
  TTree *tinProbe = dynamic_cast<TTree *>(finProbe.Get(treeName.Data()));
  if (!tinProbe) {
    fprintf(stderr, "trackGreta: tree \"%s\" not found in %s\n", treeName.Data(), inFile.Data());
    return 1;
  }
  const Long64_t nIn = tinProbe->GetEntries();
  finProbe.Close();

  Long64_t iStart = entryMin;
  Long64_t iStop = (entryMax >= 0 && entryMax < nIn) ? entryMax + 1 : nIn;
  if (iStart < 0) iStart = 0;
  if (iStart >= nIn) {
    fprintf(stderr, "trackGreta: entryMin %lld >= entries %lld\n", (long long) iStart, (long long) nIn);
    return 1;
  }
  if (iStop > nIn) iStop = nIn;

  const Long64_t nProc = iStop - iStart;
  fprintf(stderr, "trackGreta: %lld entries (%lld .. %lld) from %s", (long long) nProc,
	  (long long) iStart, (long long) iStop - 1, inFile.Data());
  if (nThreads > 1) fprintf(stderr, ", %d threads", nThreads);
  fprintf(stderr, "\n");

  if (nThreads > 1) PrepareRootForParallel(inFile, treeName, iStart);

  if (nThreads <= 1) return RunSerial(inFile, outFile, treeName, iStart, iStop);
  return RunParallel(inFile, outFile, treeName, iStart, iStop, nThreads);
}
