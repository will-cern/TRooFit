{
  //start by creating the data and signal histograms in each region (using the 2D histograms)
  TFile f1("$SourceArea/TRooFit/demo/ABCD/abcdInputHistograms.root");
  
  const char* regions[4] = {"A","B","C","D"};
  //Load data and signal histograms
  TH1* h_data[4];
  TH1* h_signal[4];

  //define the edges of the various regions
  double dLeft = 0;
  double dRight = 50;
  double aLeft = dRight;
  double aRight = 100;
  
  double dBottom = 0;
  double dTop = 60;
  double aBottom = dTop;
  double aTop = 100;

  for(int i=0;i<4;i++) {
    h_data[i] = new TH1D(Form("data_%s",regions[i]),Form("data_%s",regions[i]),1,(i<2)?aLeft:dLeft,(i<2)?aRight:dRight);
    h_signal[i] = new TH1D(Form("sig_%s",regions[i]),Form("sig_%s",regions[i]),1,(i<2)?aLeft:dLeft,(i<2)?aRight:dRight);
  }
  
  TH2* hdata_unblinded = (TH2*)f1.Get("data_unblinded");
  TH2* hsig = (TH2*)f1.Get("sig");
  
  //loop over bins of 2D histograms and fill 1D histograms with them
  for(int i=1;i<=hdata_unblinded->GetNbinsX();i++) {
    for(int j=1;j<=hdata_unblinded->GetNbinsY();j++) {
      double x = hdata_unblinded->GetXaxis()->GetBinLowEdge(i);
      double y = hdata_unblinded->GetYaxis()->GetBinLowEdge(j);
      int r;
      if(x>=aLeft && y>=aBottom && x<aRight && y<aTop) r = 0;      //region A
      else if(x>=aLeft && y>=dBottom && x<aRight && y<dTop) r = 1; //region B
      else if(x>=dLeft && y>=aBottom && x<dRight && y<aTop) r = 2; //region C
      else if(x>=dLeft && y>=dBottom && x<dRight && y<dTop) r = 3;  //region D
      else continue; //not in any of the regions
      
      if(hdata_unblinded->GetBinContent(i,j)) h_data[r]->Fill(x, hdata_unblinded->GetBinContent(i,j));
      if(hsig->GetBinContent(i,j)) h_signal[r]->Fill(x, hsig->GetBinContent(i,j));
      //if( (r==2||r==0) && hdata_unblinded->GetBinContent(i,j)) h_signal[r]->Fill(x,hdata_unblinded->GetBinContent(i,j));
    }
  }

  
  
  
  
  
  TRooABCD abcd("abcd","my abcd");
  

  //add data for regions B C and D
  //if we wanted to unblind, we would 'AddData' with region A (i=0) too
  for(int i=1;i<4;i++) {
    abcd.AddData(i,h_data[i]);
  }
  
  
  //if we want to add signal shape we do the following
  
  //add signal
  /*for(int i=0;i<4;i++) {
    abcd.AddSignal(i,h_signal[i]);
  }*/
  
  //if we want to use asimov data (corresponding to mu=0) in the signal region we do:
  //abcd.AddAsimovData(0,0 /*signal strength*/);
  
  //if we want this to be a validation, we can add 'signal region' data, will be displayed in blue
  //will not be used in the fit but will be used in computing the pvalues
  //abcd.AddValidationData(0,h_data[0]);
  
  //if we want to add a relative uncertainty  of e.g. 4% to bkg prediction in region A
  //we add a scale factor with value 1 and uncert = 0.04 to region 0
  //abcd.AddBkgScaleFactor(0, 1.0, 0.043); //adds a 4.3% uncert to bkg prediction in region A
  
  
  //You can exclude parameters from the fit by setting them constant
  //By default, m1 is constant, which is the transfer factor gradient parameter
  //abcd.GetParameter("m1")->setConstant(false); //switches on linear transfer factor model
  //abcd.GetParameter("m1")->setStringAttribute("fixedValue","0"); //will add a systematic variation corresponding to fixing m1=0
  //abcd.GetParameter("mu")->setVal(0);abcd.GetParameter("mu")->setConstant(true); //disable signal
  
  abcd.Draw(); //draws the state of the model before running the fit
  
  TRooFitResult* result = abcd.Fit(); //fit function will draw postfit results
  
  //TCanvas pullPlot; result->Draw(); //to draw the pull plot of constrained parameters (on its own canvas)
  
  //compare the result to simple ABCD arithmetic method prediction
  double NB = abcd.GetDataHistogram(1)->Integral();
  double NC = abcd.GetDataHistogram(2)->Integral();
  double ND = abcd.GetDataHistogram(3)->Integral();
  double simplePrediction = NB * NC/ND;
  double simplePredictionError = sqrt( pow(NC/ND,2)*NB + pow(NB/ND,2)*NC + pow(NB*NC/(ND*ND),2)*ND ); //standard error propagation
  cout << "simple Prediction = " << simplePrediction << " +/- " << simplePredictionError 
       << " vs TRooABCD = " << abcd.GetBkgIntegral(0) << " +/- " << abcd.GetBkgIntegralError(0) << std::endl;
  
  
}