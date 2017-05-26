{
  RooRealVar x("x","x",0,10);
  TRooH1D h("h","h",x,2); h.Fill(6,2); h.SetFillColor(kRed);
  TRooH1D h2("h2","h2",x,2); h2.Fill(2); h2.Fill(6); h2.SetFillColor(kBlue);
  
  TRooHStack hsum("hsum","hsum"); hsum.Add(&h); hsum.Add(&h2);

  hsum.GetStack()->Draw("hist");

  TCanvas cc;
  //try to stack a hist with a gaussian  
  TRooH1D h3("h3","h3",x,5);
  //h3.Fill(3,4);h3.Fill(7,2);h3.Fill(9,1);h3.SetFillColor(kBlue);
  h3.SetBinContent(2,4);h3.SetBinContent(4,2);h3.SetBinContent(5,1);h3.SetFillColor(kBlue);
  
  RooRealVar y("y","y",5);
  RooRealVar z("z","z",2);
  RooGaussian g("g","g",x,y,z);
  
  RooConstVar ch("1","1",1);
  RooRealVar cg("cg","cg",1,0,100);
  RooRealSumPdf dd("dd","dd",RooArgList(h,g),RooArgList(ch,cg));
  fr = x.frame();
  dd.plotOn(fr); 
  fr->Draw(); //has the right shape but cannot be 'extended' :-(
  
  //try instead a flat TRooH1D scaled by a gaussian!
  TRooH1D s("s","signal strength",x,1);
  /*s.Fill(3,5);*/s.SetBinContent(1,5);s.SetFillColor(kRed);
  s.AddNormFactor(g); //FIXME: when scaling by a continuous variable like this, integrations should no longer be performed binwise!!
  
  //Try a RooProdPdf ...
  TRooH1D s2("s2","signal strength",x,1);
  /*s.Fill(3,5);*/s2.SetBinContent(1,5);s2.SetFillColor(kRed);
  RooProdPdf s2g("s2g","s2g",RooArgList(s2,g));
  

  TRooHStack hs("hs","hs"); hs.Add(&s);hs.Add(&h3);
  
  RooCurve curve_hs(hs,x,0,10,100);
  curve_hs.Draw("ALP"); //this draws the DENSITY

  x.setBinning(RooUniformBinning(0,10,100),"myBinning");
  hs.SetRangeName("myBinning");
  hs.GetStack()->Draw(); //this draws the content (density*binVolume) following myBinning
  
  //generating toys efficiently means you must switch the binning first ...
  //otherwise will sample the pdf with the default (100 bin) binning
  x.setBinning(x.getBinning("myBinning"));
  toy = hs.generate( x , RooFit::Extended() ); //poisson fluctuates about "4"
  toy = hs.generate( x , RooFit::Extended() ); //poisson fluctuates about "4"
  
  
  TGraphErrors toy_data;
  for(int i=0;i<toy->numEntries();i++) {
    double x = toy->get(i)->getRealValue("x");
    toy_data.SetPoint(i,x,toy->weight());
    toy_data.SetPointError(i,0,sqrt(toy->weight()));
  }
  TCanvas toy_canvas;
  toy_data.Draw("AP");
  
  //add a signal norm factor and then fit the total model to the data, to extract fitted norm factor
  s.AddNormFactor(cg);
  
  model = hs.buildConstraints(*toy->get(0),"",true);
  
  
  
  
}