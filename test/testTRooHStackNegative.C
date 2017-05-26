{
  //attempts to test behaviour when one of the subcomponents is 'negative'
  RooRealVar x("x","x",0,10);
  TRooH1D h("h","h",x,2); h.Fill(6,2); h.SetFillColor(kRed);
  TRooH1D h2("h2","h2",x,2); h2.Fill(2); h2.Fill(6,-1); h2.SetFillColor(kBlue);
  
  RooMsgService::instance().addStream( RooFit::DEBUG ); //turns on debug messages
  
  TRooHStack hsum("hsum","hsum"); hsum.Add(&h); hsum.Add(&h2);

}