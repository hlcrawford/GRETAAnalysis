/**********************************************************************/
/* File: SortingStructures.C                                          */
/* Description: Functions for GRETINA/Aux control variables and       */
/*              counting variables in analysis                        */
/* Author: H. Crawford                                                */
/* Date: January 2013                                                 */
/**********************************************************************/

#include <stdlib.h>

#include "SortingStructures.h"

ClassImp(controlVariables);

/****************************************************/
/* Control variable functions...                    */
/****************************************************/

controlVariables::controlVariables() { 
  specifyCalibration = 0;
  calibrationFile = "";
  inputFile = "";
  rootFile = "";
  holeNum = -1;
  calibrationRun = 0;
  hfc = 0;
}

/****************************************************/

void controlVariables::Initialize() {  
  specifyCalibration = 0;
  calibrationFile = "";
  inputFile= "";
  rootFile = "";
  holeNum = -1;
  calibrationRun = 0;
  hfc = 0;
}

/****************************************************/

Int_t controlVariables::InterpretCommandLine(int argc, char *argv[]) {
  Int_t i = 1;
  while (i < argc) {
    TString flag(argv[i]);
    if (flag == "-f") {
      inputFile = argv[i + 1];
      printf("Input file: %s\n", inputFile.Data());
      i += 2;
    } else if (flag == "-rootFile") {
      rootFile = argv[i + 1];
      i += 2;
      printf("ROOT file will be saved as %s\n", rootFile.Data());
    } else if (flag == "-readCal") {
      specifyCalibration = 1;
      calibrationFile = argv[i + 1];
      i += 2;
    } else if (flag == "-hole") {
      holeNum = atoi(argv[i + 1]);
      i += 2;
    } else if (flag == "-calibrationRun") {
      calibrationRun = 1;
      i++;
    } else if (flag == "-withHFC") {
      hfc = 1;
      i++;
    } else  {
      cout << "Error -- unrecognized input flag: " << flag << endl;
      return -1;
    }
  }
  return 1;
}
