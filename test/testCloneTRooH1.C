{
  RooCategory cat("cat","cat");cat.defineType("CR");cat.defineType("SR");
  TRooH1 h("h","h",cat,{2},{0});
  h.Fill("SR",5);
  RooRealVar alpha("alpha","alpha",-5,5);
  h.addParameter(alpha);
  
  alpha.setVal(1);
  h.Fill("SR",7); //SR value is 7 if alpha = 1
  alpha.setVal(-1);
  h.Fill("SR",4); //SR value is 4 if alpha = -1
  
  cat.setIndex(1);
  
  TRooH1* h2 = (TRooH1*)(h.clone("h2"));
  
  h2->GetHistogram()->Draw();

}