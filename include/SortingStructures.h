#ifndef __SORTINGSTRUCTURES_H
#define __SORTINGSTRUCTURES_H

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"
#include <vector>
#include <stdint.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

using std::cout; using std::endl;
using std::cin; using std::cerr;

/* Analysis control variables...   */

class controlVariables : public TObject {
  
 public:
  Int_t specifyCalibration;
  TString calibrationFile;
  TString inputFile;
  TString rootFile;
  Int_t config;
  Int_t holeNum;
  Int_t calibrationRun;
  
 public:
  controlVariables();
  ~controlVariables() { ; }
  void Initialize();
  Int_t InterpretCommandLine(int argc, char *argv[]);  
  
  private:  
    ClassDef(controlVariables, 1);
};

#endif
