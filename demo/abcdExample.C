{
  //This example is for a simple cut-and-count ABCD simultaneous fit
  //This can be used to obtain a background (and even a signal) estimate in signal region A
  


  //Need a dummy observable, 
  RooRealVar x("x","x",4,0,4);
  
  //Define a signal and background histogram (4 bins)
  TRooH1D s("s","signal",x,4);s.SetFillColor(kRed);
  TRooH1D b("b","background",x,4);b.SetFillColor(kCyan);
  
  //Set the signal 'yield' in each region
  s.SetBinContent(1,94.2786);s.SetBinContent(2,2.456);
  s.SetBinContent(3,0.12);s.SetBinContent(4,1.02);
  
  //Statistical uncertainties can be added like this:
  s.SetBinError(1,3.3);s.SetBinError(2,2.2);
  s.SetBinError(3,0.05);s.SetBinError(4,0.44);
  
  //A systematic uncertainty can be added by defining a 
  //nuisance parameter, and making s depend on it, like this:
  
  RooRealVar syst_np("syst_np","Systematic nuisance parameter",0,-5,5);
  syst_np.setStringAttribute("constraintType","normal"); //means syst_np = 0 +/- 1
  
  s.addParameter(syst_np);
  
  syst_np = 1;  //corresponds to the up-fluctuation values
  s.SetBinContent(1,96.);s.SetBinContent(2,3.);
  s.SetBinContent(3,0.10);s.SetBinContent(4,0.5);
  syst_np = -1; //now set the down-fluctuations
  s.SetBinContent(1,96.);s.SetBinContent(2,3.);
  s.SetBinContent(3,0.10);s.SetBinContent(4,0.5);
  
  syst_np = 0; //just reset to nominal, for convenience
  
  //normFactor for signal .. will be determined by fit
  RooRealVar mu("mu","mu",1,0,100);
  s.addNormFactor(mu);
  
  //example of global scale factor with uncertainty .. a 10% lumi uncertainty
  RooRealVar lumi("lumi","lumi",1,0,2);
  lumi.setStringAttribute("constraintType","gaussian(1,0.1)");
  s.addNormFactor(lumi);
  
  
  //Background is to be determined by fitting to data
  //But we will need add 'shapeFactors' (bin-specific normFactors) 
  //which will ensure that A/B = C/D 
  
  //all bkg bins will be determined by shapeFactors
  b.SetBinContent(1,1);b.SetBinContent(2,1);
  b.SetBinContent(3,1);b.SetBinContent(4,1);

  //shapeFactors for bkg
  RooRealVar bkg_A("bkg_A","bkg_A",1,0,100);
  RooRealVar tBA("tBA","Bkg Ratio: B/A",1,0,1000);
  RooRealVar bkg_C("bkg_C","bkg_C",1,0,100);
  b.addShapeFactor(1,bkg_A); //bkg in signal region is a free parameter
  b.addShapeFactor(2,tBA);b.addShapeFactor(2,bkg_A); //make region B = tBA*bkg_A 
  b.addShapeFactor(3,bkg_C); //bkg in region C is a free parameter
  b.addShapeFactor(4,tBA);b.addShapeFactor(4,bkg_C); //make region D = tBA*bkg_C

  //Combine s and b
  TRooHStack h("h","signal+background");h.SetLineColor(kBlue);
  h.Add(&b);h.Add(&s);

  h.Draw("e3005"); //draws with shaded error bar

  //data
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(x,w),"weightVar");
  x=0.5;data.add(x,41);//events in signal region ... 
  x=1.5;data.add(x,21);//events in region B
  x=2.5;data.add(x,6);//events in region C
  x=3.5;data.add(x,10);//events in region D

  model = TRooFit::BuildModel(h,data); //this adds the nuisance parameter constraint terms to the model (h)
  RooFitResult* r = model->fitTo(data,RooFit::Save());

  h.Draw("e3005"); //draw post-fit state
  TH1* data_hist = data.createHistogram("x",4); //fills a histogram with content of RooFit RooDataSet
  data_hist->SetMarkerStyle(20);data_hist->Draw("same"); //draw data .. to confirm the fit worked!
  
  
  TCanvas cc;
  TRooFitResult fr(r); 
  fr.Draw(); //draw the nuisance parameter pull plot


  std::cout << "bkg in signal region:" << b.GetBinContent(1) << " +/- " << b.GetBinError(1) << std::endl;
  std::cout << "signal in region:" << s.GetBinContent(1) << " +/- " << s.GetBinError(1) << std::endl;

}