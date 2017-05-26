{
  
  
  TRooH0D h_sr("h_sr","h_sr");
  h_sr.SetBinContent(5);
  RooRealVar mu("mu","mu",1,0,100);
  h_sr.addNormFactor(mu);
  
  TRooH0D h_data("h_data","h_data");
  h_data.setData();
  h_data.Fill(10);
  
  
  h_sr.fitTo(*h_data.data()); //should give mu=2

}