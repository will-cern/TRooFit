{
  ///This is an example of a model for a simple
  ///s + b 'counting' model, with an uncertainty 
  ///on b, represented by a nuisance parameter: theta
  ///and a signal strength parameter of interest: mu
  
  
  //Start with a dummy observable, which RooFit needs for now
  //The observable we will call 'region', imagining it represented
  //a signal region
  RooCategory region("region","region");
  region.defineType("SR"); //signal region
  
  //define s and b components
  TRooH1D s("s","signal",region);
  TRooH1D b("b","background",region);
  
  //and parameter of interest 
  RooRealVar mu("mu","mu",1,-10,100); //note: we allow 'negative' values in the model
  //and nuisance parameter for uncertainty on b
  RooRealVar theta("theta","theta",0,-5,5);
  theta.setStringAttribute("constraintType","normal"); //normal gaussian constraint
  
  //set nominal values
  s.SetBinContent("SR",1); //nominal signal = 1 (with no stat uncert)
  b.SetBinContent("SR",4); //nominal bkg = 4
  
  //add variation of b, under theta
  b.addParameter(theta);
  theta = 1;
  b.SetBinContent("SR",5); //+1sigma bkg = 5
  theta = -1;
  b.SetBinContent("SR",3.5); //-1 sigma bkg=3.5 ... asymmetric uncertainty!
  theta = 0; //move back to nominal

  //add mu as a normFactor on s
  s.addNormFactor(mu);
  
  //Combine s and b, using a TRooHStack
  TRooHStack sb("sb","Signal+Background");
  sb.Add(&b);sb.Add(&s);
  sb.setFloor(true); //ensures overall pdf must be positive

  
  ///Done! ... 
  ///the complete model, including constraint term on theta
  ///Can be accessed with:
  sb.model(); //is a RooAbsPdf&
  
  
  //Model Visualization:
  //Here are some examples of how to visualise the model:
  b.SetFillColor(kBlue);
  s.SetFillColor(kRed);
  sb.SetLineColor(kGreen);
  
  sb.Draw(); //show the nominal stack
  sb.Draw("theta=2",""); //show what the model is when theta=2
  sb.Draw("mu=1",""); //show model with mu=1
  sb.SetMaximum(20);
  sb.Draw("mu=2","hist same"); //draw histogram on top corresponding to mu=2


  //Model fitting:
  //Here's how you would fit to some data
  RooRealVar w("weightVar","weightVar",1);
  RooDataSet data("data","data",RooArgSet(region,w),"weightVar");
  data.add(region,10); //suppose we have 10 events 
  
  auto res = sb.model().fitTo(data,RooFit::Save());
  TRooFitResult tres(res);
}