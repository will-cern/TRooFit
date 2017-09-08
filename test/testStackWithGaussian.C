{
  RooRealVar x("x","x",0,10);
  
  TRooH1D h3("h3","h3",x,5);
  h3.Fill(3,4);h3.Fill(7,2);h3.Fill(9,1);h3.SetFillColor(kBlue);
  
  RooRealVar y("y","y",5);
  RooRealVar z("z","z",2/*0.1*/);
  
  RooGaussian g("g","g",x,y,z);
  RooRealVar s("s","signal strength",5);
  RooExtendPdf gs("gs","gs",g,s);
  
  TRooHPdfStack hs("hs","hs");hs.Add(&h3);hs.Add(&gs);

  //compare stack to a RooAddPdf:
  RooAddPdf hs2("hs2","hs2",RooArgList(gs,h3));
/*
  fr = x.frame();
  //FIXME:: why does plotting only work correctly if plot gs before plotting hs2?
  hs2.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected));
  ///gs.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected),RooFit::LineStyle(kDashed));
  ///hs2.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected),RooFit::LineColor(kGreen));
  //hs2.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected),RooFit::Components("h3"),RooFit::LineStyle(kDashed));
  ///h3.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected),RooFit::LineColor(kRed));
  //hs2.plotOn(fr,RooFit::Normalization(1.0,RooAbsReal::RelativeExpected));
  fr->Draw();*/
  
  
  //to plot the density, you would do:
  //RooArgSet ss(x);
  //RooCurve c2(hs,x,0,10,100,hs.expectedEvents(x),&ss);
  
  hs.Draw();
  
  
  //create a fitresult with all parameters ... used to draw error bar
  RooArgSet* params = hs.getParams(x);
  TRooFitResult r(*params);
  delete params;
  
  hs.SetLineColor(kRed);hs.SetFillColor(kRed-9);hs.SetLineWidth(2);
  hs.Draw(r,"pdf A3L"); //draws the error bar of the pdf
  hs.Draw(r,"pdf LX same"); //draws the central value line
  
  
}