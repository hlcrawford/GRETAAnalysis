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
}

/****************************************************/

void controlVariables::Initialize() {  
  specifyCalibration = 0;
  calibrationFile = "";
  inputFile= "";
  rootFile = "";
  holeNum = -1;
  calibrationRun = 0;
}

/****************************************************/

Int_t controlVariables::InterpretCommandLine(int argc, char *argv[]) {
  Int_t i=1;
  while (i < argc) {
    if (strcmp(argv[i], "-f") == 0) {
      inputFile = argv[i+1];
      printf("Input file: %s\n", inputFile.Data());
      i+=2;
    } else if (strcmp(argv[i], "-rootFile") == 0) {
      rootFile = argv[i+1];
      i+=2;
      printf("ROOT file will be saved as %s\n", rootFile.Data());
    } else if (strcmp(argv[i], "-readCal") == 0) {
      specifyCalibration = 1; calibrationFile = argv[i+1];
      i+=2;
    } else if (strcmp(argv[i], "-hole") == 0) {
      holeNum = atoi(argv[i+1]);
      i+=2;
    } else if (strcmp(argv[i], "-calibrationRun") == 0) {
      calibrationRun = 1;
      i+=2;    
    } else {
      cout << "Error -- unrecognized input flag: " << argv[i] << endl;
      return -1;
    }
  }
  return 1;
}
