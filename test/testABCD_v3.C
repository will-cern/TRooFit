{
  //In this construction of the ABCD method
  //the degrees of freedom are muB and muD (content of region B and D)
  //and factor tCD (where prediction for region C = tCD*muD).
  //factor tABCD (where prediction for region A = muB*tCD*tABCD) is fixed at "1"
  
  //will do in a single histogram ..
  RooRealVar cat("cat","cat",0,4);

  
  TRooH1D h("h","Regions",cat,4);
  
  //fill each cat with '1' for now, will use binwise scalefactors to define relationships
  h.SetBinContent(1,1);h.SetBinContent(2,1);h.SetBinContent(3,1);h.SetBinContent(4,1);
  
  RooRealVar muB("muB","muB",22,0,1000);
  RooRealVar muD("muD","muD",100,0,1000);
  
  h.addShapeFactor( 2, muB ); h.addShapeFactor( 4, muD ); //scales those bins by the given factors
  
  RooRealVar tCD("tCD","tCD",0.03,0,2);
  
  h.addShapeFactor( 3, tCD ); h.addShapeFactor( 3, muD ); //region C is now tCD*muD
  
  RooRealVar tABCD("tABCD","tABCD",1); //constant=1
  
  h.addShapeFactor( 1, tABCD); h.addShapeFactor(1, tCD); h.addShapeFactor(1, muB); //region A is now tABCD*tCD*muB
  
  
  //build the dataset
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(cat,w),"weightVar");
  cat=0.5;data.add(cat,100);
  cat=1.5;data.add(cat,22);//events in region B
  cat=2.5;data.add(cat,3);//events in region C
  cat=3.5;data.add(cat,100);//events in region D
  
  //now run a fit ... in only the range cat = [B,C,D]
  
  cat.setRange("CRs",1,4);
  
  res = h.fitTo(data,RooFit::Save(),RooFit::Range("CRs")); //should leave hA = 22 * 3/100 = 0.66
  
}