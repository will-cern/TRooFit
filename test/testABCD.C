{
  //Region C and D are used to measure a transfer factor (C/D) 
  //which is then applied to the measurement in region B to get estimate in region A
  //Degrees of freedom are on region D, the transfer factor (C/D), and region B
  
  
  //create regions
  //normally would fill with the MC-based estimate 
  //but in absence of any estimate, we will just use uniform everywhere
  TRooH0D hA("hA","Region A",1);
  TRooH0D hB("hB","Region B",1);
  TRooH0D hC("hC","Region C",1);
  TRooH0D hD("hD","Region D",1);
  
  //add normFactors (which will float in the fit)
  RooRealVar muB("muB","muB",1,0,100);hB.addNormFactor(muB);
  RooRealVar muD("muD","muD",1,0,100);hD.addNormFactor(muD);
  
  //create transfer factors between B->A and D->C
  
  auto transAB = hA.createTransFactor( &hB ); //represents A/B
  auto transCD = hC.createTransFactor( &hD ); //represents C/D
  
  //add a normFactor to transCD
  RooRealVar muCD("muCD","muCD",1,0,100); transCD->addNormFactor(muCD);
  
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
    
  
  RooSimultaneous simPdf("simPdf","simPdf",RooArgList(hB,hC,hD),cat);

  //take some data and fit!
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(cat,w),"weightVar");
  cat.setLabel("B");data.add(cat,22);//events in region B
  cat.setLabel("C");data.add(cat,3);//events in region C
  cat.setLabel("D");data.add(cat,100);//events in region D

  res = simPdf.fitTo(data,RooFit::Save()); //should leave hA = 22 * 3/100 = 0.66
  TRooFitResult tres(res);
  tres.Draw();
}