{
  //here we will construct the classic
  //mu + b single bin (i.e. cut and count) model
  //with a gaussian uncertainty on b
  //will set a limit on mu
  
  double nObs = 0;
  
  double b_nom = 1.; //nominal background
  double b_relUncert = 0.1; //10% uncert
  
  //need a roofit variable for the signal strength (mu)
  RooRealVar mu("mu","Signal Strength",1,0,100);
  
  //need a dummy observable (especially if going to use RooStats)
  RooRealVar x("x","x",0,1);
  
  //we want two single-bin TRooH1Ds, one for signal, one for bkg ...
  TRooH1D s("s","signal",x,1); s.SetFillColor(kRed);
  TRooH1D b("b","background",x,1); b.SetFillColor(kCyan);
  
  //set up signal ... we want it to just equal mu
  //lots of ways to achieve this. One way is to set the bin content
  //to 1, and then add mu as a normalization factor ...
  s.SetBinContent(1,1);
  s.addNormFactor(mu);
  
  //now add the nominal background
  b.SetBinContent(1,b_nom);
  
  //make b depend on a nuisance parameter that is 0 when b=nom, +1 for 1 sigma, -1 for -1sigma
  RooRealVar alpha_b("alpha_b","alpha_b",-5,5);
  alpha_b.setStringAttribute("constraintType","normal"); //means syst_np = 0 +/- 1 (i.e. normal distribution)

  b.addParameter(alpha_b,4); //4 is the interpolation mode: 6th order polynomial with linear extrapolation
  
  alpha_b = 1;b.SetBinContent(1,b_nom*(1.+b_relUncert));
  
  //could do an asymmetric error like this:
  //alpha_b = -1;b.SetBinContent(1,b_nom*(1.-0.5*b_relUncert));
  
  alpha_b = 0; //just return to nominal position
  
  
  //put in a stack
  TRooHStack sb("sb","s + b");
  sb.Add(&b);sb.Add(&s);
  
  //visualize the model
  sb.SetMinimum(0);
  sb.Draw("e3005"); //drawing with e3005 means an error band (e) is drawn, with shading pattern 3005

  //can use fillGraph method to inspect pdf value as function of variables
  //fillGraph will sample the TRooFit object at (by default, 100) points across the range of the variable
  TGraph g1;
  sb.fillGraph(&g1,mu);
  TCanvas cg1; g1.Draw("AL");
  
  TGraph g2;
  sb.fillGraph(&g2,alpha_b);
  TCanvas cg2; g2.Draw("AL");


  //Define a dataset ... the data will just be an event count
  //the the events will have had the same value of x (=0.5)
  RooRealVar weightVar("weightVar","weightVar",1); //needed to give entry a weight
  RooDataSet data("data","data",RooArgSet(x,weightVar),"weightVar");
  
  //can add data like this:
  data.add(x,nObs);


  //we use the dataset to call BuildModel ... which will identify the constrained and unconstrained PARAMETERS
  //constrained parameters will have constraint terms created for them
  //constrained parameters have 'constraintType' attribute set
  
  RooAbsPdf* model = TRooFit::BuildModel( sb, data );
  model->SetName("myModel");


  //can shove everything into a workspace now, if we want to save to a file ..
  RooWorkspace w("w","w");
  w.import(*model);
  w.import(data);

  //inspect what ended up in the workspace with:
  //w.Print()
  
  
  //to use the model in RooStats, you have to prepare a ModelConfig ...
  auto mc = TRooFit::CreateModelConfig(w, "myModel","data","mu");
  
  w.import(*mc);
  
  w.writeToFile("/tmp/myModelWorkspace.root");
  
  //can then run the standard ROOT limit setting macro:
  /*
  .L $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C+
  StandardHypoTestInvDemo("/tmp/myModelWorkspace.root","w","ModelConfig","","data",
      2,  //Asymptotic calculator
      3,  //One-sided, lower bound test statistic (i.e. q^tilde_mu
      true,  //do CLS
      100,0,10) //number of points to scan, and min and max mu

//Result is: (for b= 1 +/- 0.1)
//  The computed upper limit is: 2.33574 +/- 0
//Expected upper limits, using the B (alternate) model : 
// expected limit (median) 3.39876
// expected limit (-1 sig) 2.13919
// expected limit (+1 sig) 5.59963
// expected limit (-2 sig) 1.44866
// expected limit (+2 sig) 8.97334


  */

}