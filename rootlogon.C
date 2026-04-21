#include "TString.h"

TFile *fIn;

TTree* openFileTEB(TString fname) {
  fIn = new TFile(fname.Data());
  TTree *teb;
  fIn->GetObject("teb", teb);
  return teb;
}

void rootlogon () {
  TString gretaLib = "lib/libGRETA.so";
  const char *searchP = "./";
  const char* foundlib;
  foundlib = gSystem->Which(searchP, gretaLib, EAccessMode::kFileExists);
  if (foundlib) {
    gSystem->Load("lib/libGRETA.so");
    cout << "Loaded: libGRETA.so";
    cout << endl;
  } else {
    gretaLib = "lib/libGRETA.dll";
    foundlib = gSystem->Which(searchP, gretaLib, EAccessMode::kFileExists);
    if (foundlib) {
      gSystem->Load("lib/libGRETA.dll");
      cout << "Loaded: libGRETA.dll" << endl;
    } else {
      gretaLib = "lib/libGRETA.dylib";
      foundlib = gSystem->Which(searchP, gretaLib, EAccessMode::kFileExists);
      if (foundlib) {
	gSystem->Load("lib/libGRETA.dylib");
	cout << "Loaded: libGRETA.dylib" << endl;
      } else { 
	cout << "Could not load shared libraries!" << endl;
      }
    }
  }

  foundlib = gSystem->Which(searchP, "traceUtilities_cpp.so", EAccessMode::kFileExists);
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
