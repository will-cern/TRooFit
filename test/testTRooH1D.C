{
  RooCategory cat("cat","cat");cat.defineType("CR");cat.defineType("SR");
  
  TRooH1 h("h","h",cat,{2},{0});
  
  h.SetBinContent(2,5);//h.Fill("SR",5);//h.Fill("SR",4);
  
  //add an alpha nuisance parameter
  RooRealVar alpha("alpha","alpha",-5,5);
  h.AddParameter(alpha);
  
  alpha.setVal(1);
  h.Fill("SR",7); //SR value is 7 if alpha = 1
  alpha.setVal(-1);
  h.Fill("SR",4); //SR value is 4 if alpha = -1
  alpha.setVal(0);alpha.setError(1);
  h.GetHistogram()->Draw();

  fr = alpha.frame();
  cat.setIndex(1);
  h.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::Raw)); //seems to plot h(alpha)/(integral h(alpha))

  //test plotting of a TRooH1D
  RooRealVar x("x","x",0,10);
  TRooH1D h2("h2","h2",x,2);
  h2.Fill(1,2);
  fr2 = x.frame();
  h2.plotOn(fr2);

  //error propagation test
  TRooFitResult r(*h.getParameters(cat));
  r.Print(); //shows the values that will be used in the error propagation
  

}