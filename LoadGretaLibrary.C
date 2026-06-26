/**
 * Load libGRETA and register ROOT streamer dictionaries (for interactive ROOT).
 *
 * Usage from GRETAAnalysis / TrackBackToMe project directory:
 *   .x LoadGretaLibrary.C
 */
void LoadGretaLibrary(const char *libdir = "lib")
{
  gSystem->AddDynamicPath(libdir);

  TString inc = gSystem->Getenv("ROOT_INCLUDE_PATH");
  if (inc.Length()) inc += ":";
  inc += libdir;
  gSystem->Setenv("ROOT_INCLUDE_PATH", inc.Data());

  /* Prefer native shared library (avoid stale Linux .so on macOS). */
  TString libPath;
#if defined(__APPLE__)
  libPath = TString::Format("%s/libGRETA.dylib", libdir);
#elif defined(__linux__)
  libPath = TString::Format("%s/libGRETA.so", libdir);
#else
  libPath = TString::Format("%s/libGRETA", libdir);
#endif

  if (gSystem->Load(libPath.Data()) < 0) {
    /* Fallback for platforms that build without extension in Load(). */
    if (gSystem->Load(TString::Format("%s/libGRETA", libdir).Data()) < 0) {
      ::Error("LoadGretaLibrary", "failed to load GRETA library from %s", libdir);
      return;
    }
  }

  if (!TClass::GetClass("g2OUT") || !TClass::GetClass("GretaTrackedEvent")) {
    ::Warning("LoadGretaLibrary",
	      "library loaded but g2OUT/GretaTrackedEvent dictionary missing; "
	      "rebuild with: scons readGreta");
    return;
  }
  printf("LoadGretaLibrary: dictionaries OK (g2OUT, GretaTrackedEvent, vector<TVector3>)\n");
}
