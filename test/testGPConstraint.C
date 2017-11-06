{
  RooRealVar x("x","x",0,10);
  
  TRooH1D h("h","h",x,5); //5 bins 
  h.SetBinContent(1,3);h.SetBinContent(2,2);
  
  TH1D* h_ref = new TH1D("x","x",5,0,10);
  h_ref->SetBinContent(5,4);
  for(int i=1;i<=5;i++) h_ref->SetBinError(i,1.);
  
  RooDataHist h_ref_d("h_ref","h_ref",x,h_ref);

  h.getVal(x);

  h_ref_d.attachBuffers(x);

  for(int i=0;i<h_ref_d.numEntries();i++) {
    const RooArgSet* a = h_ref_d.get(i);
    h.setValueDirty();
    std::cout << "x=" << x.getVal() << " ref = " << h_ref_d.weight() << " h = " << h.getVal(x) << std::endl;
  }


  TMatrixD k(5,5); for(int i=0;i<5;i++) k(i,i)=1.;


  TRooGPConstraint c("c","c",h,h_ref_d,k);
  std::cout << "CP constraint val = " << c.getLogVal() << std::endl;
  
  
  TRooChi2Constraint c2("c2","c2",h,h_ref_d,true);
  std::cout << "Chi2 constraint val = " << c2.getLogVal() << std::endl;

}