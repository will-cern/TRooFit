{
  //Region C and D are used to measure a transfer factor (C/D) 
  //which is then applied to the measurement in region B to get estimate in region A
  //Degrees of freedom are on region D, the transfer factor (C/D), and region B
  
  
  //create regions
  //normally would fill with the MC-based estimate 
  //but in absence of any estimate, we will just use uniform everywhere
  
  //RooFit currently struggles to cope when there are no observables
  RooCategory dummy("dummy","dummy");dummy.defineType("dummy");dummy.setConstant(1);
  
  
  TRooH1D hA("hA","Region A",dummy);hA.SetBinContent(1,1);
  TRooH1D hB("hB","Region B",dummy);hB.SetBinContent(1,1);
  TRooH1D hC("hC","Region C",dummy);hC.SetBinContent(1,1);
  TRooH1D hD("hD","Region D",dummy);hD.SetBinContent(1,1);hD.SetBinError(1,0.1);
  
  //add normFactors (which will float in the fit)
  RooRealVar muB("muB","muB",1,0,100);hB.addNormFactor(muB);
  RooRealVar muD("muD","muD",1,0,100);hD.addNormFactor(muD);
  
  //create transfer factors between B->A and D->C
  
  auto transAB = hA.createTransFactor( &hB ); //represents A/B ... makes hA return transAB*hB, and moves all uncerts out of hA and hB into transAB
  auto transCD = hC.createTransFactor( &hD ); //represents C/D
  
  //add a normFactor to transCD
  RooRealVar muCD("muCD","muCD",1,0,100); transCD->addNormFactor(muCD); // ... effectively like adding a norm factor on C
  
  //finally create a transfer factor between the transfactors: CD -> AB
  //this represents the 'rho' factor: (A/B)/(C/D)
  auto transABCD = transAB->createTransFactor( transCD );
  //since transABCD = 1 here, it basically constrains A/B = C/D
  
  //connect up all regions up with a RooSimultaneous
  RooCategory cat("cat","cat");
  //cat.defineType("A"); -- doing a blinded background estimate, so no A type
  cat.defineType("B");
  cat.defineType("C");
  cat.defineType("D");
  

  RooSimultaneous simPdf("simPdf","simPdf",RooArgList(hB.model(),hC.model(),hD.model()),cat);

//   //take some data and fit!
//   FIXME: need to find a good way to include the dummy observable into this!
//   TRooH1D data("data","data",cat);data.isData();
//   data.Fill("B",22);
//   data.Fill("C",3);
//   data.Fill("D",100);
  
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(dummy,cat,w),"weightVar");
  cat.setLabel("B");data.add(cat,22);//events in region B
  cat.setLabel("C");data.add(cat,3);//events in region C
  cat.setLabel("D");data.add(cat,100);//events in region D

  res = simPdf.fitTo(data,RooFit::Save()); //should leave hA = 22 * 3/100 = 0.66
  TRooFitResult tres(res);
  tres.Draw();
}