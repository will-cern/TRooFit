{
  //From the presentation given on 12th September 2017 in the stat forum

  RooRealVar x("x","x",0,8);
  TRooH1D h("h","h",x,4); //4 bin histogram
  h.Draw();
  h.Print(); //lists dependent args
  
  h.SetBinContent(3,3);
  h.SetLineColor(kBlue);
  h.SetLineWidth(2);
  
  h.GetBinContent(3); //returns 3
  x = 6; //move x to 6 (center of bin)
  h.getVal();  //RooFit method: returns 1.5
  h.getVal(x); //RooFit method: returns 0.5 
  
  h.Fill(3,2); //also creates a stat error
  
  h.GetBinContent(2); //= 2
  h.GetBinError(2);   //= 2
  h.GetBinContent(2,"h_stat_bin2=1.5"); //=3
  
  RooRealVar lumi("lumi", "lumi",1,0,2);
  lumi.setStringAttribute("constraintType",
                                          "gaussian(1,0.1)");
  
  h.addNormFactor(lumi);
  
  h.Draw("","lumi=0.5");
  
  RooRealVar alpha("alpha", "alpha",0,-5,5);
  alpha.setStringAttribute("constraintType",
                                          "normal");
  
  h.addParameter(alpha);
  alpha=1.0;
  h.SetBinContent(3,2);
  h.SetBinContent(2,3);
  
  h.Draw("","alpha=0.5"); //both bins=2.5


 RooRealVar mu("mu", "mu",1,0,100);

 h.addNormFactor(mu);
 
 //construct a RooFit dataset
 RooRealVar w("weightVar","weight",1);  
 RooDataSet data("data","data",
	RooArgSet(x,w),"weightVar");  
 x=3; data.add(x,52);
 x=5; data.add(x,60);


 model = TRooFit::BuildModel(h,data);

 res =  model->fitTo(data,RooFit::Save());

 h.Draw(); //post-fit histogram

 TRooFitResult r(res);
 r.Draw();

}
