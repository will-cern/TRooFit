{
  RooRealVar bin1("bin1","bin1 truth content",0,100);
  RooRealVar bin2("bin2","bin2 truth content",0,100);

  RooRealVar x("x","x",0,10);
  TRooH1D h("h","h",x,2);h.addNormFactor(bin1); h.SetFillColor(kRed);
  TRooH1D h2("h2","h2",x,2);h2.addNormFactor(bin2); h2.SetFillColor(kBlue);
  h.SetBinContent(1,3.0/8.1); h.SetBinError(1,0.1/8.2);   h.SetBinContent(2,0.1/8.1); //bin migration
  h2.SetBinContent(2,4/8.1);  h2.SetBinContent(1,0.1/8.1); //bin migration

  TRooHPdfStack hs("hs","hs");
  hs.Add(&h);hs.Add(&h2);
  
  RooDataHist data("data","data",RooArgSet(x),"h");
  x=3; data.add(x,11);
  x=7; data.add(x,17);
  
  TRooFitResult tres( hs.model().fitTo(data, RooFit::Save() ) );
  
  //in this setup, the bin content was filled with efficiencies
  //therefore the parameters of interest are themselves the bin contents
  std::cout << "bin1 truth content = " << bin1.getVal() << std::endl;
  std::cout << "bin2 truth content = " << bin2.getVal() << std::endl;

  hs.Draw();

  tres.Draw();

}