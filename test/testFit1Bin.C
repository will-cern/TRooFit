{
  
  
  TRooH1 h_sr("h_sr","h_sr");
  h_sr.SetBinContent(1,5);
  RooRealVar mu("mu","mu",1,0,100);
  h_sr.AddNormFactor(mu);
  
  RooRealVar w("weightVar","weight",1);
  RooDataSet data("data","data",RooArgSet(w),"weightVar");
  data.add(RooArgSet(),10);
  
  h_sr.fitTo(data); //should give mu=2

}