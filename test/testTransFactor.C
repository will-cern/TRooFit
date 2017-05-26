{
  
  
  TRooH0D h_sr("h_sr","h_sr");
  h_sr.Fill(5);
  
  
  //add an alpha nuisance parameter
  RooRealVar alpha("alpha","alpha",-5,5);
  h_sr.AddParameter(alpha);
  
  alpha.setVal(1);
  h_sr.Fill(7); //SR value is 7 if alpha = 1
  alpha.setVal(-1);
  h_sr.Fill(4); //SR value is 4 if alpha = -1
  alpha.setVal(0);alpha.setError(1);
  
  
  TRooH0D h_cr("h_cr","h_cr");
  h_cr.Fill(50);
  h_cr.AddParameter(alpha);
  alpha.setVal(1);h_cr.Fill(60);
  alpha.setVal(-1);h_cr.Fill(40);
  
  auto tFactor = h_sr.createTransFactor( &h_cr );
  
  RooRealVar mu("mu","mu",1,0,100);
  h_cr.AddNormFactor(mu);
  

  RooDataSet data("data","data",RooArgSet());
  data.add(RooArgSet(),2);
  

}