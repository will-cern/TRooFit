{
  //we will set up an s+b model with two bins
  //uncertainty on b will be a NormFactor constrained by a gaussian (e.g. a lumi measurement)
  
  //bkg = {950,0}
  //sig = {300k,200k};
  //data = {1050,0};
  
  //old
  //bkg = {10,0}
  //sig = {30,20};
  //data = {12,0};
  
  RooRealVar x("x","x",0,10);
  
  TRooH1D bkg("bkg","bkg",x,2);bkg.SetFillColor(kBlue); bkg.SetBinContent(1,10); //bkg.SetBinContent(2,1);
  TRooH1D sig("sig","sig",x,2); sig.SetFillColor(kRed); sig.SetBinContent(1,1); sig.SetBinContent(2,0.2);

  TRooHStack mean("mean","mean");mean.Add(&bkg);mean.Add(&sig);
  
  RooRealVar alpha("alpha","alpha",0,-100,100);
  bkg.addParameter(alpha);
  
  alpha = 1; bkg.SetBinContent(1,11);
  alpha = 0;
  
  RooRealVar mu("mu","mu",0,0,20);
  //RooRealVar mu("mu","mu",0,0,0.0001);
  sig.addNormFactor(mu);
  
  
 

  RooConstVar alpha_gobs("alpha_gobs","alpha_gobs",0);RooConstVar alpha_sigma("alpha_sigma","alpha_sigma",5);

  RooGaussian alpha_constraint("alpha_constraint","alpha_constraint",alpha_gobs,alpha,alpha_sigma); //5% uncert

  alpha.setConstant(1);

  RooProdPdf model("model","model",RooArgList(mean,alpha_constraint));

  //Data will be 20 events in bin 1, 0 event in bin 2
  TH1D* hdata = new TH1D("hdata","hdata",2,0,10); hdata->Sumw2(); 
  hdata->SetBinContent(1,12); //hdata->SetBinError(1,sqrt(20.)); 
  //hdata->SetBinContent(2,6);hdata->SetBinError(2,sqrt(6));
  RooDataHist data("data","data",x,hdata);
  
  
  mean.GetStack()->Draw(); //visualise model
  hdata->Draw("same");
  
  //model.fitTo(data);
  //mean.GetStack()->Draw(); //visualise model
  //hdata->Draw("same");

TCanvas dd;

//first do without systematic

  alpha.setConstant(1);
  nll = mean.createNLL(data);
  pll = nll->createProfile(mu);
  //alpha.setConstant(0);pll->Print();alpha.setVal(1);alpha.setConstant(1); .. must not do this because paramAbsMin would end up away from 1!!
  RooCurve c2(*pll,mu,mu.getMin()+1e-9,mu.getMax(),100);c2.SetLineColor(kRed);
  c2.Draw();
  
  
  //then repeat with it floating
  alpha.setConstant(0);
  RooProfileLL pll2("pll2","pll2",*nll,mu);
  RooCurve c(pll2,mu,mu.getMin()+1e-9,mu.getMax(),100);
  c.Draw("same");
 
  //cycle over mu values, conditional MLE
  TGraph gg;
  for(int i=0;i<c2.GetN();i++) {
    if(c2.GetX()[i] < mu.getMin() || c2.GetX()[i] > mu.getMax()) continue;
    mu.setVal(c2.GetX()[i]);pll2.getVal(); //will move alpha to CMLE
    gg.SetPoint(gg.GetN(),mu.getVal(),alpha.getVal());
  }
  gg.Draw("same LP");
 
  //get gradients from last 10 bins
  cout << endl;
  cout << " red (no syst) gradient = " << (c2.GetY()[c2.GetN()-6] - c2.GetY()[c2.GetN()-16])/(c2.GetX()[c2.GetN()-6] - c2.GetX()[c2.GetN()-16]) << std::endl;
  cout <<  " blue (floating) gradient = " << (c.GetY()[c.GetN()-6] - c.GetY()[c.GetN()-16])/(c.GetX()[c.GetN()-6] - c.GetX()[c.GetN()-16]) << std::endl;
 
  /*
  gradient is: s1+s2 - s1N1/(b1+s1*mu) - s2N2/(b2+s2*mu)  ... but N2 = 0 in our case
  
  so gradient = s2 + s1(1 - N1/(b1+s1*mu))
  if can fit bin 1, then gradient is s2
  
  If underfit bin1, then N1 > b1+s1*mu, leading to gradient < s2
  this can lead to impression of wider likelihood curve in this underfit (no syst) case
  
  -logL = b1+b2+mu*(s1+s2) - N1*log( (b1+s1*mu) )  (assuming N2 is 0 and neglect constants)
  
  this is a straight line with a log subtracted, i.e. an upside log curve with a rising slope to the right
  
  
  
  */
 
  //fr = mu.frame();
  //nll->plotOn(fr);
  //fr->Draw();

}