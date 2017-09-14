{
  //Demo of different interpolation modes
  //Also shows how TRooFit objects can be plotted like RooFit functions
  
  
  TRooH0D h("h","h");
  h.SetBinContent(5);
  
  RooRealVar alpha("alpha","alpha",0,-5,5);  alpha.setBins(10); //ensures alphaBinWidth = 1
  
  h.addParameter(alpha);
  
  alpha = -1;h.SetBinContent(4);
  alpha = 1;h.SetBinContent(7);
  
  //show flexibleInterpVar for comparison ... has a bug in it for code=4!!!
  RooStats::HistFactory::FlexibleInterpVar xx("xx","xx",alpha,5,{4},{7});
  xx.setInterpCode(alpha,4);


  //now we want to draw h as a function of alpha 
  auto fr = alpha.frame();
  
  xx.plotOn(fr); //FlexibleInterpVar's 6th order with log extrap
  h.setInterpCode("alpha",0); //piecewise linear (default)
  h.plotOn(fr,RooFit::Normalization(h.expectedEvents(alpha)),RooFit::LineColor(kBlack)); 
  h.setInterpCode("alpha",3); //6th order with linear extrap
  h.plotOn(fr,RooFit::Normalization(h.expectedEvents(alpha)),RooFit::LineColor(kGreen)); //plots: h.getVal(alpha)*h.expectedEvents(alpha) / alphaBinWidth
  h.setInterpCode("alpha",4); //6th order with log extrap
  h.plotOn(fr,RooFit::Normalization(h.expectedEvents(alpha)),RooFit::LineStyle(2),RooFit::LineColor(kRed));
  h.setInterpCode("alpha",2); //6th order with log extrap (alternative implementation, may be faster/slower?)
  h.plotOn(fr,RooFit::Normalization(h.expectedEvents(alpha)),RooFit::LineStyle(3),RooFit::LineColor(kViolet));
  
  
  fr->Draw();


}