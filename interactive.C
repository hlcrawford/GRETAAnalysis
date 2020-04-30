#include <iomanip>

Int_t variable;
Int_t definefitregion;
Int_t beginfitregion;
Int_t fitregionlo;
Double_t fitregionlo_use;
Int_t fitregionhi;
Double_t fitregionhi_use;
Double_t fitresults[8];
Double_t fitresults_err[8];
Double_t fitintegral;
Double_t fitintegral_err;

TCanvas *c9;

Int_t last_x;
Int_t last_y;

Int_t peakmark[500];

/***************************************************************/

Double_t back(Double_t *x, Double_t *par) {
  /* Linear background - par[0] = offset, par[1] = slope */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[0] + xx*par[1] < 0) { 
    f = 1000000000000;
  } else {
    f = (par[0] + xx*par[1]);
  }
  return f;
}

/***************************************************************/

Double_t plainGaus(Double_t *x, Double_t *par) {
  /* Gaussian - par[2] = amplitude, par[3] = amplitude scaling,
     par[4] = centroid, par[5] = sigma */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[2]*((100-par[3])/100))*
	 TMath::Exp(-0.5*((xx-par[4])/(par[5]))*((xx-par[4])/(par[5]))));
  } else {
    f = 1000000000000;
  } 
  return f;
}

/***************************************************************/

Double_t skewGaus(Double_t *x, Double_t *par) {
  /* Skewed Gaussian peak -- a la Radford gf3 */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[2]*par[3]/100)*
	 TMath::Exp((xx-par[4])/par[6])*
	 TMath::Erfc(((xx-par[4])/(TMath::Sqrt(2)*par[5])) +
		     ((par[5])/(TMath::Sqrt(2)*par[6]))));
  } else {
    f = 1000000000000;
  } 
  return f;
}

/***************************************************************/

Double_t stepFnc(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[7])*
	 TMath::Erfc((xx-par[4])/(TMath::Sqrt(2)*par[5])));
    if (f < 0) {
      f =  1000000000000;
    }
  } else {
    f = 1000000000000;
  }
  return f;
}

/***************************************************************/

Double_t skewFit(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return f;
}

/***************************************************************/

Double_t skewLimitP(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return TMath::Sqrt(f);
}

Double_t skewLimitM(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return -1*TMath::Sqrt(f);
}

/***************************************************************/

Double_t gausFit(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return (back(x, &par[0]) + plainGaus(x, &par[0]) +
	  stepFnc(x, &par[0]));
}

/***************************************************************/

Double_t gausLimitP(Double_t *x, Double_t *par) {
  return TMath::Sqrt((back(x, &par[0]) + plainGaus(x, &par[0]) +
		      stepFnc(x, &par[0])));
}


Double_t gausLimitM(Double_t *x, Double_t *par) {
  return -1*TMath::Sqrt((back(x, &par[0]) + plainGaus(x, &par[0]) +
			 stepFnc(x, &par[0])));
}

/***************************************************************/

Double_t peakOnly(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return (plainGaus(x, &par[0]) + skewGaus(x, &par[0]));
}

/***************************************************************/

Double_t dofit(Double_t ilo=-1, Double_t ihi=-1, const char *name="",char *write="n", Int_t skewYes = 1){

  Double_t lo = ilo; Double_t hi = ihi;
  
  /* To select fit region, mouse should be moved from left to right
     If the mouse is moved from right to left the fitting range will
     be inverted.  Fix this here. */
  if(lo > hi) {
    Double_t swap = lo;
    lo = hi;  hi = swap;
  }
  
  /* Find the histogram to be fit. */
  TH1F *h = (TH1F*)gROOT->FindObject(name);
  TH1F *h1 = (TH1F*)h->Clone();  

  /* Retrieve the histogram name and use it as the base for file output. */
  string histoname = h->GetName();
  histoname = histoname + ".txt";
  
  /* Retrieve the width of a bin, use bin 1 just because. */
  /** NOTE: invalid for variable width bins. **/
  Double_t binwidth = h->GetBinWidth(1);
  
  /* Set up the linear fit. */
  TF1 *backlin = new TF1("backlin",back,lo,hi, 2);
  backlin->SetLineColor(2);
  
  /* Take range passed into fitting routine and convert into bin numbers. */
  Int_t blo = Int_t (lo / binwidth) + 1;
  Int_t bhi = Int_t (hi / binwidth) + 1;
  
  /* Assume centroid of the gaussian is the bin with the maximum counts
     that is in range.... */
  Double_t maximum = 0;
  Double_t binmax = 0;
  for(Int_t i=blo; i<bhi; i++) {
    Double_t val = h->GetBinContent(i);
    if(val > maximum) { maximum = val;  binmax = i; }
  }

  printf("Starting centroid/height = %f/%f\n", maximum, binmax);
  
  /* Convert bin value back to an xaxis value. */
  Int_t guess_centroid = binmax * binwidth;
  
  /* Create Gaussian for display and intergration purposes. */
  TF1 *gausD = new TF1("gausD", plainGaus, lo, hi, 8);
  gausD->SetLineColor(4);

  TF1 *skewgausD;
  if (skewYes) {
    skewgausD = new TF1("skewgausD", skewGaus, lo, hi, 8);
    skewgausD->SetLineColor(5);
  }
  
  TF1 *stepD = new TF1("stepD", stepFnc, lo, hi, 8);
  stepD->SetLineColor(6);
  
  /* Fit the combined spectrum with a gaussian and linear curves. */
  TF1 *userfit, *peaks, *limitP, *limitM;
  if (skewYes) {
    userfit = new TF1("userfit", skewFit, lo, hi, 8);
    peaks = new TF1("peaks", peakOnly, lo, hi, 8);
    limitP = new TF1("limitP", skewLimitP, lo, hi, 8);
    limitM = new TF1("limitM", skewLimitM, lo, hi, 8);
  } else {
    userfit = new TF1("userfit", gausFit, lo, hi, 8);
    peaks = new TF1("peaks", plainGaus, lo, hi, 8);
    limitP = new TF1("limitP", gausLimitP, lo, hi, 8);
    limitM = new TF1("limitM", gausLimitM, lo, hi, 8);
  }
  
  /* Create arrays for storing values. */
  Int_t npar = 8;
  Double_t params[8] = {0};
  Double_t params_err[8];
    
  /* Set the total fit function parameters to start the fit. */
  userfit->SetParameters(params);
 
  /* Background offset and slope */
  userfit->SetParLimits(0, -1., 10000); 
  userfit->SetParameter(0, h->GetBinContent(blo));
  userfit->SetParLimits(1, -100., 1.);
  userfit->SetParameter(1, 0.0);
    
  /* Total height, H */
  userfit->SetParLimits(2, 0.5*maximum, 2*maximum);
  userfit->SetParameter(2, maximum);

  /* Ratio of heights, R */
  userfit->SetParLimits(3, 0, 100);
  userfit->SetParameter(3, 50.);
  if (!skewYes) { userfit->FixParameter(3, 0); }

  /* Peak Position */
  userfit->SetParLimits(4, 0.99*binmax, 1.01*binmax);
  userfit->SetParameter(4, binmax);

  /* Peak width (sigma) */
  userfit->SetParLimits(5, binmax*0.0001, binmax*0.02);
  userfit->SetParameter(5, binmax*0.001);
  
  /* Beta for skewness */
  userfit->SetParLimits(6, 0.00001, 12000);
  userfit->SetParameter(6, 50.);
  if (!skewYes) { userfit->FixParameter(6, 0); }
  
  /* Step height */
  userfit->SetParLimits(7, 0, h->GetBinContent(blo)*1.1);
  userfit->SetParameter(7, h->GetBinContent(blo));
  
  //userfit->GetParameters(&params[0]);

  /* Perform the total fit within a range, all weights set to 1 
     with loglike method. */
  h->Fit("userfit","WWLMRQ");
  
  /* Retrieve the fit parameters and associated errors. */
  userfit->GetParameters(&params[0]);
  params_err[0] = userfit->GetParError(0);
  params_err[1] = userfit->GetParError(1);
  params_err[2] = userfit->GetParError(2);
  params_err[3] = userfit->GetParError(3);
  params_err[4] = userfit->GetParError(4);
  params_err[5] = userfit->GetParError(5);
  params_err[6] = userfit->GetParError(6);
  params_err[7] = userfit->GetParError(7);
  cout << "Background: offset " << params[0] << " +/- " << params_err[0] << endl;
  cout << "Background: slope " << params[1] << " +/- " << params_err[1] << endl;
  cout << "Peak height: " << params[2] << " +/- " << params_err[2] << endl;
  cout << "Ratio of heights (R): " << params[3] << " +/- " << params_err[3] << endl;
  cout << "Peak Centroid: " << params[4] << " +/- " << params_err[4] << endl;
  cout << "Peak Sigma: " << params[5] << " +/- " << params_err[5] << endl;
  cout << "Peak Skewness (beta): " << params[6] << " +/- " << params_err[6] << endl;
  cout << "Step height: " << params[7] << " +/- " << params_err[7] << endl;
  
  peaks->SetParameters(&params[0]);

  Double_t maxValue = peaks->GetMaximum();
  Double_t fwhm, fwtm;
  Double_t valueLow = peaks->GetX(maxValue*0.5, lo, params[4]);
  Double_t valueHigh = userfit->GetX(maxValue*0.5, params[4], hi);
  fwhm = valueHigh - valueLow;
  cout << "FWHM: " << fwhm << endl;
  valueLow = peaks->GetX(maxValue*0.1, lo, params[4]);
  valueHigh = peaks->GetX(maxValue*0.1, params[4], hi);
  fwtm = valueHigh - valueLow;
  cout << "FW(1/10)M: " << fwtm << endl;
  cout << " Ratio: " << fwtm/fwhm << endl;

  /* Retrieve linear parameters. */
  backlin->SetParameters(&params[0]);
  gausD->SetParameters(&params[0]);
  if (skewYes) { skewgausD->SetParameters(&params[0]); }
  stepD->SetParameters(&params[0]);
  limitP->SetParameters(&params[0]);
  limitM->SetParameters(&params[0]);

  TF1 *fPO = new TF1("fPO", peakOnly, lo, hi, 8);
  fPO->SetParameters(&params[0]);
  maxValue = fPO->GetMaximum(lo, hi);
  Double_t maxValueX = fPO->GetMaximumX(lo, hi);
  Double_t low = fPO->GetX(maxValue/2., lo, maxValueX);
  Double_t high = fPO->GetX(maxValue/2., maxValueX, hi);
  low = fPO->GetX(maxValue/10., lo, maxValueX);
  high = fPO->GetX(maxValue/10., maxValueX, hi);
  
  cout << std::setw(15) << params[2] << std::setw(15) << params[3] << std::setw(15) << params[4] << std::setw(15)
       << params[5] << std::setw(15) << params[6] <<  std::setw(15) << params[7] << endl;
 
  Int_t np = 10000;
  Double_t *xD = new Double_t[np];
  Double_t *wD = new Double_t[np];
  Double_t *xS = new Double_t[np];
  Double_t *wS = new Double_t[np];
  gausD->CalcGaussLegendreSamplingPoints(np,xD,wD,1e-15);
  if (skewYes) { skewgausD->CalcGaussLegendreSamplingPoints(np,xS,wS,1e-15); }
    
  /* Draw the functions. */
  backlin->Draw("same");
  gausD->Draw("same");
  if (skewYes) { skewgausD->Draw("same"); }
  stepD->Draw("same");
  h1->Add(userfit, -1);
  for (Int_t i=0; i<h1->GetNbinsX(); i++) {
    if (h1->GetBinCenter(i) < ilo || h1->GetBinCenter(i) > ihi) {
      h1->SetBinContent(i,0);
    }
  }

  c9->cd(2);
  h1->Draw();

  limitP->SetLineStyle(2);
  limitP->SetLineWidth(1);
  limitP->SetLineColor(kGray+2);
  limitM->SetLineStyle(2);
  limitM->SetLineWidth(1);
  limitM->SetLineColor(kGray+2);

  limitP->Draw("Csame");
  limitM->Draw("Csame");
  
  /* If function dofit is invoked from command line and option write
     is invoked write out the results of the fit. */
  if(write=="w") {
    ofstream efile(histoname.c_str(),ios_base::app); //append to end of file
    
    if(!efile){
      return 1;
    } else {
      efile << setw(12) << params[0] << " " << setw(12) << params_err[0] << " " 
	    << setw(12) << params[1] << " " << setw(12) << params_err[1] << " "
	    << setw(12) << params[2] << " " << setw(12) << params_err[2] << " "
	    << setw(12) << params[3] << " " << setw(12) << params_err[3] << " "
	    << setw(12) << params[4] << " " << setw(12) << params_err[4] << " "	   
	    << setw(12) << params[5] << " " << setw(12) << params_err[5] << " "
	    << setw(12) << params[6] << " " << setw(12) << params_err[6] << " "
	    << setw(12) << params[7] << " " << setw(12) << params_err[7] << " "
	    << setw(12) << (gausD->Integral(lo,hi)/binwidth + skewgausD->Integral(lo,hi)/binwidth);
    }
    
    efile.close();
  }
  
  /* Put fit parameters in global variables to enable later access. */
  for(Int_t j=0;j<8;j++){
    fitresults[j] = params[j];
    fitresults_err[j] = params_err[j];
  }
  return 0;
}

Double_t doIntegral(Double_t lo=-1, Double_t hi=-1, const char *name="",
		    char *write="n"){
    
    /* To select fit region, mouse should be moved from left to right
       If the mouse is moved from right to left the fitting range will
       be inverted.  Fix this here. */
    if (lo > hi) { Double_t swap = lo;  lo = hi;  hi = swap; }
    
    /* Find the histogram to be fit. */
    TH1F *h = (TH1F*)gROOT->FindObject(name);

    /* Retrieve the histogram name and use it as the base for file output. */
    string histoname = h->GetName();
    histoname = histoname + ".txt";
    
    /* Retrieve the width of a bin, use bin 2 just because. */
    /** NOTE: invalid for variable width bins. **/
    Double_t binwidth = h->GetBinWidth(2);

    /* Set up the linear background. */
    TF1 *backlin = new TF1("pol1","pol1",lo,hi);
    backlin->SetLineColor(2);
    
    /* Take range passed into fitting routine and convert into bin numbers. */
    Int_t blo = Int_t (lo / binwidth) + 1;
    Int_t bhi = Int_t (hi / binwidth) + 1;
    
    Float_t lowEdge=0, hiEdge=0;
    for (Int_t i=blo; i<blo+3; i++) {
      lowEdge += h->GetBinContent(i);
    } 
    for (Int_t i=bhi; i>bhi-3; i--) {
      hiEdge += h->GetBinContent(i);
    } 
    lowEdge /= 3; hiEdge /= 3;

    /* Figure out the linear background */
    Float_t slope = (lowEdge - hiEdge)/(h->GetBinCenter(blo+1) - h->GetBinCenter(bhi-1));
    Float_t offset = (hiEdge - slope*h->GetBinCenter(bhi-1));

    backlin->SetParameters(offset, slope);

    /* Calculate the integral - linear background contribution */
    Float_t integral = 0;
    for (Int_t binNum = blo; binNum <= bhi; binNum++) {
      integral += (h->GetBinContent(binNum) - backlin->Eval(h->GetBinCenter(binNum), 0, 0));
    }

    backlin->Draw("same");
    return 0;
}

// Write out results of the most recent fit to a file

void writeresults(){
  // Retrieve the 1D histogram that is contained in the canvas
  TH1F *temp;
  TObject *obj;
  TIter next (c9->GetListOfPrimitives());
  while( (obj=next()) ){
    if(obj->InheritsFrom("TH1")) { temp = (TH1F*)obj; }
  }
  
  // Retrieve a name for file output
  string histoname = temp->GetName();
  histoname = histoname + ".txt";
  histoname = "fits.out";
  
  ofstream efile(histoname.c_str(),ios_base::app);
  
  if(!efile){
    return 1;
  }
  else{
    efile << "# " << temp->GetTitle() << " peak at ch " << fitresults[3] << endl;
    efile << setw(12) << fitresults[0] << " " << setw(12) << fitresults_err[0] << " "
	  << setw(12) << fitresults[1] << " " << setw(12) << fitresults_err[1] << " "
	  << setw(12) << fitresults[2] << " " << setw(12) << fitresults_err[2] << " "
	  << setw(12) << fitresults[3] << " " << setw(12) << fitresults_err[3] << " "
	  << setw(12) << fitresults[4] << " " << setw(12) << fitresults_err[4] << " "
	  << setw(12) << fitresults[5] << " " << setw(12) << fitresults_err[5] << " "
	  << setw(12) << fitresults[6] << " " << setw(12) << fitresults_err[6] << " "
	  << setw(12) << fitresults[7] << " " << setw(12) << fitresults_err[7] << " "
	  << setw(12) << fitintegral   << " " << setw(12) << fitintegral_err << endl;
  }
  
  efile.close();  
}

void testobjects(Double_t lo=-1, Double_t hi=-1, char *write="n"){
  
  TH1F *temp;
  TIter next (gPad->GetListOfPrimitives());
  TObject *obj;
  while((obj=next())){
    if(obj->InheritsFrom("TH1")){
      temp = (TH1F*)obj;
    }
  }

  int binlo, binhi;
  Double_t minx, maxx;
  
  binlo = binhi = -1;
  minx = maxx = 0;
  
  minx = lo;
  maxx = hi;
  
  if((maxx - minx) > 1000){
    return 1;
  }
  
  dofit(minx, maxx,temp->GetName(), write, 1);
  gPad->Modified();
  gPad->Update();
}


void interactive()
{
   // Create a new canvas which will have interactive
   // properties attached to it
   c9 = new TCanvas("c9","Interactive Fitting",10,10,700,500);
   c9->SetFillColor(10);
   c9->SetFrameFillColor(10);

   c9->Divide(1,2);
   TPad *c9_1, *c9_2;
   c9_1 = (TPad*)c9->GetPrimitive("c9_1");
   c9_2 = (TPad*)c9->GetPrimitive("c9_2");
   c9_1->SetPad(0., 0.2, 1.0, 1.0);
   c9_2->SetPad(0., 0.0, 1.0, 0.2);

   c9->cd(1);

   variable = 0;
   definefitregion = 0;

   fitregionlo = 0;
   fitregionhi = 0;
   fitregionlo_use = 0;
   fitregionhi_use = 0;
   
   // Add a TExec object to the canvas
   c9_1->AddExec("dynamic1","DynamicExec1()");
}

void DynamicExec1()
{
    // Example of function called when a mouse event occurs in a pad.

    // This portion of code derives from DynamicSlice example in the
    // ROOT tutorials and will perform actions for different mouse
    // and keyboard actions as follows.
    // "i" - zoom in by a factor of 2 (defined by variable expandfactor)
    // "o" - zoom out by a factor of 2 (also defined by variable expandfactor)
    // "l" - move spectrum to left
    // "r" - move spectrum to right
    // "p" - mark peak location
    // "u" - remove marked peak
    // 1st left click - begin to define fit region, fit region will display as gray box
    // 2nd left click - end defining fit region, fit then executed
    // middle button  - save most recent fit results to file

    // Flag to check if mouse event is over an axis
    Int_t axiscom;

    // Get the object the mouse is over
    TObject *select = gPad->GetSelected();
    if(!select) return;
    //if (!select->InheritsFrom("TFrame")) {gPad->SetUniqueID(0); return;}
    if (select->InheritsFrom("TAxis")) { axiscom = 1; }
    else { axiscom = 0; }
    TH1 *h = (TH1*)select;
    gPad->GetCanvas()->FeedbackMode(kTRUE);
    
    // If a left click is detected that is not over an axis and a fit region
    // has been created fit the region and reset the flag "definefitregion"
    if(gPad->GetEvent() == kButton1Down && definefitregion == 1 && axiscom == 0){
      definefitregion = 0;
      cout << "I will fit now from " << fitregionlo << " " << fitregionhi << endl;
      testobjects(gPad->AbsPixeltoX(fitregionlo), gPad->AbsPixeltoX(fitregionhi));
      //dofit(gPad->AbsPixeltoX(fitregionlo),gPad->AbsPixeltoX(fitregionhi),h->GetName());
      
    }
    else if(gPad->GetEvent() == kButton1Down && definefitregion == 0 && axiscom == 0){
      // if a left click is detected that is not over an axis and the flag
      // "definefitregion" is zero start defining the fit region
      definefitregion = 1;
      gPad->SetUniqueID(0);
    }
    
    if(gPad->GetEvent() == 12/*kButton3Down*/){
      // if the middle mouse button is pressed wirte the output
      cout << "middle button, writing output " << endl;
      writeresults();
    }

    if(gPad->GetEvent() == 24){ // Keyboard pressed
      float epsilon = 0.000000001;
      float expandfactor = 2;
      
      // Retrieve the histogram
      TH1F *temp;
      TIter next (gPad->GetListOfPrimitives());
      TObject *obj;
      while((obj=next())){
	if(obj->InheritsFrom("TH1")){
	  temp = (TH1F*)obj;
	}
      }
      
      if(gPad->GetEventX() == 105){ // "i" on keyboard, zoom in by factor 2
	float uxmin = gPad->GetUxmin();
	float uxmax = gPad->GetUxmax();
	float totrange = uxmax - uxmin;
	float center = (totrange / 2) + uxmin;
	float range = (center - uxmin) / expandfactor;
	float rlo = center - range;
	float rhi = center + range;
	temp->GetXaxis()->SetRangeUser(rlo,rhi-epsilon);
	gPad->Modified();
	gPad->Update();
      }
      
      if(gPad->GetEventX() == 111){ // "o" on keyboard, zoom out by factor 2
	float uxmin = gPad->GetUxmin();
	float uxmax = gPad->GetUxmax();
	float totrange = uxmax - uxmin;
	float center = (totrange / 2) + uxmin;
	float range = (center - uxmin) * expandfactor;
	float rlo = center - range;
	float rhi = center + range;
	temp->GetXaxis()->SetRangeUser(rlo,rhi-epsilon);
	gPad->Modified();
	gPad->Update();   
      }
      
      if(gPad->GetEventX() == 108){ // "l" on keyboard, move to left
	float uxmin = gPad->GetUxmin();
	float uxmax = gPad->GetUxmax();
	float totrange = uxmax - uxmin;
	float rlo = uxmin - totrange;
	float rhi = uxmin;
	temp->GetXaxis()->SetRangeUser(rlo,rhi);
	gPad->Modified();
	gPad->Update();   
      }

      if(gPad->GetEventX() == 114){ // "r" on keyboard, move to right
	float uxmin = gPad->GetUxmin();
	float uxmax = gPad->GetUxmax();
	float totrange = uxmax - uxmin;
	float rlo = uxmax;
	float rhi = uxmax+totrange;
	temp->GetXaxis()->SetRangeUser(rlo,rhi);
	gPad->Modified();
	gPad->Update();   
      }
      
      Int_t numpeaks = 0; // FIX THIS - HLC 10/18/19

      if (gPad->GetEventX() == 112) { /* "p" on keyboard, mark peak */
	peakmark[numpeaks] = last_x;
	numpeaks++;
	/* Draw arrow at peak position, just above peak. */
      }
      
      if (gPad->GetEventX() == 117) { /* "u" on keyboard, unmark peak */
	Double_t diff[500] = {0};
	Double_t low_diff = 0;
	for (Int_t i=0; i<500; i++) {
	  diff[i] = abs((last_x) - peakmark[i]);
	  if (diff[i] < low_diff) {
	    low_diff = diff[i];
	  }
	}
	for (Int_t i=0; i<500; i++) {
	  if (diff[i] == low_diff) {
	    peakmark[i] = 0;
	  }
	}
	/* Delete arrow... */
      }
      
      if(gPad->GetEventX() == 115){
	writeresults();
      }

    } else { /* If not a keyboard press... */
      /* Update x,y coordinates for next keyboard press. */
      last_x = gPad->GetEventX();
      last_y = gPad->GetEventY();
    }
    
    // In between two left clicks draw a box which graphically
    // shows the fit region on the histogram
    if(definefitregion == 1 && gPad->GetEvent() == kMouseMotion){

      // Expanding box in x dimension
      // Erase old position and draw a box at current position
       int pxold = gPad->GetUniqueID();
       int px = gPad->GetEventX();
       int py = gPad->GetEventY();
       float uxmin = gPad->GetUxmin();
       float uxmax = gPad->GetUxmax();
       float uymin = gPad->GetUymin();
       float uymax = gPad->GetUymax();
       int pxmin = gPad->XtoAbsPixel(uxmin);
       int pxmax = gPad->XtoAbsPixel(uxmax);
       int pymin = gPad->YtoAbsPixel(uymin);
       int pymax = gPad->YtoAbsPixel(uymax);
       if(pxold == 0) {
           fitregionlo = px;
           fitregionlo_use = gPad->GetEventX();
       }
       if(pxold) gVirtualX->DrawBox(fitregionlo,pymin,pxold,pymax,TVirtualX::kFilled);
       gVirtualX->DrawBox(fitregionlo,pymin,px,pymax,TVirtualX::kFilled);
       fitregionhi = px;
       fitregionhi_use = gPad->GetEventX();
       gPad->SetUniqueID(px);
       Float_t upx = gPad->AbsPixeltoY(px);
       Float_t x = gPad->PadtoY(upx);
   }
}
