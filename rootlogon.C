#include "TString.h"

TFile *fIn;

TTree *openFileTEB(TString fname)
{
  fIn = new TFile(fname.Data());
  TTree *teb = nullptr;
  fIn->GetObject("teb", teb);
  return teb;
}

void rootlogon()
{
  /* Full dictionary setup (see LoadGretaLibrary.C). */
  gROOT->ProcessLine(".x LoadGretaLibrary.C");

  const char *searchP = "./";
  const char *foundlib = gSystem->Which(searchP, "traceUtilities_cpp.so", kFileExists);
  if (foundlib) {
    gSystem->Load("traceUtilities_cpp.so");
    cout << "Loaded traceUtilities." << endl;
  }

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
}
