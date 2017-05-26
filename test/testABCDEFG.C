{
  //this example supposes there are two source of bkg
  //bkg1 and bkg2
  //bkg2 contributes to regions ABC EFG
  //bkg1 contributes only to regions ABCD
  
  //create regions
  //normally would fill with the MC-based estimate 
  //but in absence of any estimate, we will just use uniform everywhere
  TRooH0D hA1("hA1","Region A bkg1",1);
  TRooH0D hB1("hB1","Region B bkg1",1);
  TRooH0D hC1("hC1","Region C bkg1",1);
  TRooH0D hD1("hD1","Region D bkg1",1);
  
  TRooH0D hA2("hA2","Region A bkg2",1);
  TRooH0D hB2("hB2","Region B bkg2",1);
  TRooH0D hC2("hC2","Region C bkg2",1);
  
  TRooH0D hE2("hE2","Region E bkg2",1);
  TRooH0D hF2("hF2","Region F bkg2",1);
  TRooH0D hG2("hG2","Region G bkg2",1);
  
  //transfer factors: E->A, F->B, G->C ... will all be fixed @ '1' because of above values
  auto transAE = hA2.createTransFactor( &hE2 ); //represents A/E
  auto transBF = hB2.createTransFactor( &hF2 ); //represents B/F
  auto transCG = hC2.createTransFactor( &hG2 ); //represents C/G
  
  //add normFactors to B C D E F and G .. i.e to let float in fit
  //6 floats, 6 datapoints
  RooRealVar muB("muB","muB",1,0,100);hB1.AddNormFactor(muB);//not added to B2, since that comes from F2 and transfer factor
  //RooRealVar muC("muC","muC",1,0,100);hC2.AddNormFactor(muC); ... this dof is added onto the transCD factor instead
  RooRealVar muD("muD","muD",1,0,100);hD1.AddNormFactor(muD);
  RooRealVar muE("muE","muE",1,0,100);hE2.AddNormFactor(muE);
  RooRealVar muF("muF","muF",1,0,100);hF2.AddNormFactor(muF);
  RooRealVar muG("muG","muG",1,0,100);hG2.AddNormFactor(muG);
  
  
  
  //create transfer factors between B1->A1 and D1->C1
  
  auto transAB = hA1.createTransFactor( &hB1 ); //represents A/B
  auto transCD = hC1.createTransFactor( &hD1 ); //represents C/D ... because of below, C1 becomes = muCD*(1)*D1
  
  //add a normFactor to transCD
  RooRealVar muCD("muCD","muCD",1,0,100); transCD->AddNormFactor(muCD);
  
  //finally create a transfer factor between the transfactors: CD -> AB
  //this represents the 'rho' factor: (A/B)/(C/D)
  auto transABCD = transAB->createTransFactor( transCD );
  //since transABCD = 1 in this example (because of numbers above), this effectively makes transAB = transCD = muCD
  
  //create combined estimates from region A, B, C:
  TRooHStack hA("hA","hA");hA.Add(&hA1);hA.Add(&hA2);
  TRooHStack hB("hB","hB");hB.Add(&hB1);hB.Add(&hB2);
  TRooHStack hC("hC","hC");hC.Add(&hC1);hC.Add(&hC2);
  
  
  //connect up all regions up with a RooSimultaneous
  RooCategory cat("cat","cat");
  //cat.defineType("A"); -- doing a blinded background estimate, so no A type
  cat.defineType("B");
  cat.defineType("C");
  cat.defineType("D");
  cat.defineType("E");
  cat.defineType("F");
  cat.defineType("G");
  
  RooSimultaneous simPdf("simPdf","simPdf",RooArgList(hB,hC,hD1,hE2,hF2,hG2),cat);

  //take some data and fit!
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(cat,w),"weightVar");
  cat.setLabel("B");data.add(cat,22);//events in region B
  cat.setLabel("C");data.add(cat,3);//events in region C
  cat.setLabel("D");data.add(cat,100);//events in region D
  cat.setLabel("E");data.add(cat,50);//events in region E (transfers straight to A)
  cat.setLabel("F");data.add(cat,2);//events in region F (transfers straight to B)
  cat.setLabel("G");data.add(cat,1);//events in region G (transfers straight to C)

  //this gives C1 = C-G = 2
  //B1 = B-F = 20
  //C1/D = 2/100 = 0.02
  //A1 = B1*C1/D = 0.02*20 = 0.4
  //So estimate in A = 0.4 + 50 = 50.4

  auto fr = simPdf.fitTo(data,RooFit::Save()); 

}