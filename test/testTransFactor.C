{
  
  
  TRooH0D h_sr("h_sr","h_sr");
  h_sr.Fill(5);
  
  
  //add an alpha nuisance parameter
  RooRealVar alpha("alpha","alpha",0,-5,5);
  h_sr.addParameter(alpha);
  
  alpha = 1;
  h_sr.Fill(7); //SR value is 7 if alpha = 1
  alpha = -1;
  h_sr.Fill(4); //SR value is 4 if alpha = -1
  alpha = 0;
  
  //normally constrain it:
  alpha.setStringAttribute("constraintType","normal");
  
  
  TRooH0D h_cr("h_cr","h_cr");
  h_cr.Fill(50);
  h_cr.addParameter(alpha);
  
  alpha=1;h_cr.Fill(60);
  alpha=-1;h_cr.Fill(40);
  
  auto tFactor = h_sr.createTransFactor( &h_cr ); //this makes sr = tFactor * cr_nom
  tFactor->Print(); //will have value of 0.1 for alpha<=0, and varying above that (alpha 7/60 at alpha=1)
  
  RooRealVar mu("mu","mu",1,0,100);
  h_cr.addNormFactor(mu);
  

  h_sr.Draw("mu=1,alpha=0",""); //should be 5
  h_sr.Draw("mu=1,alpha=-1",""); //still 5 (since value of cr_nom is used, alpha is ignored except in tFactor)
  h_sr.Draw("mu=2,alpha=-1",""); //is 10 ... cr_nom gets multiplied by 2, and then transferred
  

}