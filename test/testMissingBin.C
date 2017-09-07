{
  RooRealVar bin1("bin1","bin1 truth content",0,100);
  RooRealVar bin2("bin2","bin2 truth content",0,100);

  RooRealVar x("x","x",0,10);
  TRooH1D h("h","h",x,2);h.addNormFactor(bin1); h.SetFillColor(kRed);
  TRooH1D h2("h2","h2",x,2);h2.addNormFactor(bin2); h2.SetFillColor(kBlue);
  h.SetBinContent(1,3); h.SetBinError(1,0.1);   h.SetBinContent(2,0.1); //bin migration
  h.SetMissingContent(5);
  h2.SetBinContent(2,4);  h2.SetBinContent(1,0.1); //bin migration
  h2.SetMissingContent(4);

  TRooHPdfStack hs("hs","hs");
  hs.Add(&h);hs.Add(&h2);
  
  RooDataHist data("data","data",RooArgSet(x),"h");
  x=3; data.add(x,11);
  x=7; data.add(x,17);
  
  TRooFitResult tres( hs.model().fitTo(data, RooFit::Save() ) );
  
  //truth content is given by complete integral of each:
  cout << " truth bin1 = " << h.expectedEvents(x) + h.missingEvents() << std::endl;
  cout << " truth bin2 = " << h2.expectedEvents(x) + h2.missingEvents() << std::endl;

  hs.Draw();

  tres.Draw();

}