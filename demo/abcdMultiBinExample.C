{
  //This example is for a multibin ABCD simultaneous fit
  //Each region will actually be set up as a single bin histogram ... but you could naturally extend it
  //This can be used to obtain a background (and even a signal) estimate in signal region A
  
  //This example is set up to give the same results as abcdExample.C ... but can be extended
  //to multibin distributions


  //Need our actual observable (this is not the same as the ABCD plane observables)
  RooRealVar x("x","x",0,10);
  
  //Define a signal and background histogram for each region
  const char* regionLabels[4] = {"A","B","C","D"};
  TRooH1D* s[4]; TRooH1D* b[4];
  
  for(int i=0;i<4;i++) {
    s[i] = new TRooH1D(Form("s_%s",regionLabels[i]),Form("region %s signal",regionLabels[i]),x,1);s[i]->SetFillColor(kRed);         //CHANGE NUMBER OF BINS HERE!
    b[i] = new TRooH1D(Form("b_%s",regionLabels[i]),Form("region %s background",regionLabels[i]),x,1);b[i]->SetFillColor(kCyan);    //CHANGE NUMBER OF BINS HERE!
  }
  
  //Set the signal 'yield' in each region
  s[0]->SetBinContent(1,94.2786);s[1]->SetBinContent(1,2.456);
  s[2]->SetBinContent(1,0.12);s[3]->SetBinContent(1,1.02);
  
  //Statistical uncertainties can be added like this:
  s[0]->SetBinError(1,3.3);s[1]->SetBinError(1,2.2);
  s[2]->SetBinError(1,0.05);s[3]->SetBinError(1,0.44);
  
  //A systematic uncertainty can be added by defining a 
  //nuisance parameter, and making s depend on it, like this:
  
  RooRealVar syst_np("syst_np","Systematic nuisance parameter",0,-5,5);
  syst_np.setStringAttribute("constraintType","normal"); //means syst_np = 0 +/- 1
  
  for(int i=0;i<4;i++) s[i]->addParameter(syst_np);
  
  syst_np = 1;  //corresponds to the up-fluctuation values
  s[0]->SetBinContent(1,96.);s[1]->SetBinContent(2,3.);
  s[2]->SetBinContent(1,0.10);s[3]->SetBinContent(1,0.5);
  syst_np = -1; //now set the down-fluctuations
  s[0]->SetBinContent(1,96.);s[1]->SetBinContent(1,3.);
  s[2]->SetBinContent(1,0.10);s[3]->SetBinContent(1,0.5);
  
  syst_np = 0; //just reset to nominal, for convenience
  
  //normFactor for signal .. will be determined by fit
  RooRealVar mu("mu","mu",1,0,100);
   for(int i=0;i<4;i++) s[i]->addNormFactor(mu);
  
  //example of global scale factor with uncertainty .. a 10% lumi uncertainty
  RooRealVar lumi("lumi","lumi",1,0,2);
  lumi.setStringAttribute("constraintType","gaussian(1,0.1)");
  for(int i=0;i<4;i++) s[i]->addNormFactor(lumi);
  
  
  //Background is to be determined by fitting to data
  //But we will need add 'shapeFactors' (bin-specific normFactors) 
  //which will ensure that A/B = C/D 
  
  //all bkg bins will be determined by shapeFactors
  for(int i=0;i<4;i++) b[i]->SetBinContent(1,1);

  //shapeFactors for bkg
  for(int i=1;i<=b[0]->GetXaxis()->GetNbins();i++) {
    RooRealVar* bkg_A = new RooRealVar(Form("bkg_A_bin%d",i),Form("Background in region A, bin %d",i),1,0,100); //bkg in signal region is a free parameter
    b[0]->addShapeFactor(i, *bkg_A);
    RooRealVar* tBA = new RooRealVar(Form("tBA_bin%d",i),Form("Transfer factor (A->B), bin %d",i),1,0,1000);
    b[1]->addShapeFactor(i, *tBA);
    RooRealVar* bkg_C = new RooRealVar(Form("bkg_C_bin%d",i),Form("Background in region C, bin %d",i),1,0,100); //bkg in region C is a free parameter
    b[2]->addShapeFactor(i, *bkg_C);
    b[3]->addShapeFactor(i, *tBA);
  }
  
  b[1]->addNormFactor( *b[0] ); //makes region B = tBA*bkg_A
  b[3]->addNormFactor( *b[2] ); //makes region D = tBA*bkg_C
  

  //Combine s and b into 4 stacks
  //and draw them too ...
  
  TCanvas c("Prefit","Prefit distributions");
  c.Divide(2,2);
  
  TRooHStack* h[4];
  for(int i=0;i<4;i++) {
    h[i] = new TRooHStack(Form("h_%s",regionLabels[i]),Form("signal+background in region %s",regionLabels[i]));h[i]->SetLineColor(kBlue);
    h[i]->Add(b[i]);h[i]->Add(s[i]);
    c.cd(i+1);
    h[i]->SetMinimum(0);h[i]->SetMaximum(120);
    h[i]->Draw("e3005"); //draws with shaded error bar
  }
  
  
  

  //data
  
  //create a RooCategory to say which region the event has landed in
  RooCategory cat("region","region");
  for(int i=0;i<4;i++) cat.defineType(regionLabels[i]); 
  
  RooRealVar ww("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(x,ww,cat),"weightVar");
  x=5; //for this example, we will put all events at x = 5 (doesn't matter because there was only 1 bin)
  cat.setLabel("A");data.add(RooArgSet(x,cat),41);//events in signal region ... 
  cat.setLabel("B");data.add(RooArgSet(x,cat),21);//events in region B
  cat.setLabel("C");data.add(RooArgSet(x,cat),6);//events in region C
  cat.setLabel("D");data.add(RooArgSet(x,cat),10);//events in region D
  
  
  //now we construct a RooSimultaneous to perform a simultaneous fit of our 4 TRooHStacks to the data
  //To do that, we need to build the model for each TRooHStack
  
  RooArgList models;
  for(int i=0;i<4;i++) models.add( *TRooFit::BuildModel(*h[i],data) ); //adds the nuisance parameter constraint terms to each model (stack)

  RooSimultaneous model("model","model",models,cat);
  
  //RooWorkspace w; w.import(model);
  
  RooFitResult* r = model.fitTo(data,RooFit::Save());

  //Draw the postfit distributions
  TCanvas c2("Postfit","Postfit distributions");
  c2.Divide(2,2);

  for(int i=0;i<4;i++) {
    c2.cd(i+1);
    h[i]->Draw("e3005"); //draw post-fit state
    TH1* data_hist = (TH1*)s[0]->GetHist(0)->Clone(Form("data_%s",regionLabels[i])); //clones one of the source histograms to get same binning
    data_hist->Reset();
    data.fillHistogram( data_hist, x, Form("region==%d",i) );
    data_hist->SetMarkerStyle(20); data_hist->Draw("same");
    
    c.cd(i+1); //draw on the prefit too
    data_hist->Draw("same");
    
  }
  
  
  TCanvas cc("Pull","Nuisance Parameter Fit Pulls");
  TRooFitResult fr(r); 
  fr.Draw(); //draw the nuisance parameter pull plot

  
  std::cout << "bkg in signal region:" << b[0]->GetBinContent(1) << " +/- " << b[0]->GetBinError(1) << std::endl;
  std::cout << "signal in region:" << s[0]->GetBinContent(1) << " +/- " << s[0]->GetBinError(1) << std::endl;

}