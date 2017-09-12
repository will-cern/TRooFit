{
  RooRealVar x("x","x",4,0,4);
  
  TRooH1D s("s","signal",x,4);s.SetFillColor(kRed);
  TRooH1D b("b","background",x,4);b.SetFillColor(kCyan);
  
  s.SetBinContent(1,94.2786);s.SetBinContent(2,2.456);
  s.SetBinContent(3,0.12);s.SetBinContent(4,1.02);
  
  //all bkg bins will be determined by shapeFactors
  b.SetBinContent(1,1);b.SetBinContent(2,1);
  b.SetBinContent(3,1);b.SetBinContent(4,1);

  //normFactors for signal
  RooRealVar mu("mu","mu",1,0,100);
  s.addNormFactor(mu);
  RooRealVar lumi("lumi","lumi",1,0,2);
  lumi.setStringAttribute("constraintType","gaussian(1,0.1)");
  s.addNormFactor(lumi);
  
  //shapeFactors for bkg
  RooRealVar bkg_A("bkg_A","bkg_A",1,0,100);
  RooRealVar tBA("tBA","Bkg Ratio: B/A",1,0,1000);
  RooRealVar bkg_C("bkg_C","bkg_C",1,0,100);
  b.addShapeFactor(1,bkg_A);
  b.addShapeFactor(2,tBA);b.addShapeFactor(2,bkg_A);
  b.addShapeFactor(3,bkg_C);
  b.addShapeFactor(4,tBA);b.addShapeFactor(4,bkg_C);

  //Combine s and b
  TRooHStack h("h","signal+background");h.SetLineColor(kBlue);
  h.Add(&b);h.Add(&s);

  h.Draw("e3005"); //draws with shaded error bar

  //data
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(x,w),"weightVar");
  x=0.5;data.add(x,41);
  x=1.5;data.add(x,21);//events in region B
  x=2.5;data.add(x,6);//events in region C
  x=3.5;data.add(x,10);//events in region D

  model = TRooFit::BuildModel(h,data);
  model->fitTo(data);

  h.Draw("e3005"); //draw post-fit state

}