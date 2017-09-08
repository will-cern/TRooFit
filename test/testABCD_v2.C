{
  //In this construction of the ABCD method
  //the degrees of freedom are muB and muD (content of region B and D)
  //and factor tCD (where prediction for region C = tCD*muD).
  //factor tABCD (where prediction for region A = muB*tCD*tABCD) is fixed at "1"
  
  //will do in a single histogram ..
  RooCategory cat("cat","cat");
  cat.defineType("A");
  cat.defineType("B");
  cat.defineType("C");
  cat.defineType("D");
  
  TRooH1D h("h","Regions",cat);
  
  //fill each cat with '1' for now, will use binwise scalefactors to define relationships
  h.SetBinContent("A",1);h.SetBinContent("B",1);h.SetBinContent("C",1);h.SetBinContent("D",1);
  
  RooRealVar muB("muB","muB",1,0,1000);
  RooRealVar muD("muD","muD",1,0,1000);
  
  h.addShapeFactor( 2, muB ); h.addShapeFactor( 4, muD ); //scales those bins by the given factors
  
  RooRealVar tCD("tCD","tCD",1,0,2);
  
  h.addShapeFactor( 3, tCD ); h.addShapeFactor( 3, muD ); //region C is now tCD*muD
  
  RooRealVar tABCD("tABCD","tABCD",1,0,2);
  tABCD.setStringAttribute("constraintType","gaussian(1,0.1)"); //let it float
  
  h.addShapeFactor( 1, tABCD); h.addShapeFactor(1, tCD); h.addShapeFactor(1, muB); //region A is now tABCD*tCD*muB
  
  
  //build the dataset
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(cat,w),"weightVar");
  cat.setLabel("A");data.add(cat,100);
  cat.setLabel("B");data.add(cat,22);//events in region B
  cat.setLabel("C");data.add(cat,3);//events in region C
  cat.setLabel("D");data.add(cat,100);//events in region D
  
  //now run a fit ... in only the range cat = [B,C,D]
  
  cat.setRange("CRs","B,C,D");
  h.setNormRange("CRs"); //necessary for discrete variables :-(
  
  model = TRooFit::BuildModel(h,data);
  
  res = model->fitTo(data,RooFit::Save(),RooFit::Range("CRs")); //should leave hA = 22 * 3/100 = 0.66
  
}