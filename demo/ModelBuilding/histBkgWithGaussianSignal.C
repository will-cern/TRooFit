{
  //this will be a background in a histogram
  //with a gaussian signal put on top
  //for the bkg we will just sample an exponential
  
  
  TH1D* hBkg = new TH1D("hBkg","background",10,0,10);
  TF1 efunc("efunc","exp([0]+[1]*x)",0.,10.);
  efunc.SetParameters(1,-0.3);
  hBkg->FillRandom("efunc",150); //fill with 150 random bkg events from exponential
  //remove any stat uncertainty, for this example ..
  hBkg->Sumw2();
  for(int i=1;i<=hBkg->GetNbinsX();i++) hBkg->SetBinError(i,0);
  
  //need a roofit variable for the signal strength (mu)
  RooRealVar mu("mu","Signal Strength",1,0,100);
  
  //our observable
  RooRealVar x("x","x",0,10);
  
  //construct the background
  TRooH1D b("b","background",x,10); b.SetFillColor(kCyan);
  b.Add(hBkg);
  
  
  //construct a gaussian signal ...
  //and put it inside a TRooH1D ...
  //will use a RooRealVar for the mean ...allows use to vary the position of the mean later
  RooRealVar mean("mean","gaus mean",4.7);
  RooFormulaVar gaus("gaus","gaus","mu*ROOT::Math::gaussian_pdf(x,0.5,mean)",RooArgList(mu,x,mean)); //mean = 4.7, sigma = 0.5
  TRooH1D s("s","s",x,10);s.SetFillColor(kRed);
  s.Add(gaus);
  //Technical note ... when you add a function to a TRooH1D, it will assume it is not a density and will divide it by bin width
  //Adding a pdf (inherit from RooAbsPdf) will not divide by bin width
  //we can force same behaviour for a function by setting the 'isDensity' attribute:
  gaus.setAttribute("isDensity",true);
  
  
  
  //put in a stack
  TRooHStack sb("sb","s + b");
  sb.Add(&b);sb.Add(&s);
  
  //visualize the model
  sb.SetMinimum(0);
  sb.Draw("e3005"); //drawing with e3005 means an error band (e) is drawn, with shading pattern 3005
  //can visualize with a different signal strength like this:
  sb.Draw("e3005","mu=50");

  //we can overlay the PDF value (i.e. the probability DENSITY)
  //this will be scanned more finely than the binning of the TRooH1D
  //"v" draws a TGraph
  sb.Draw("v L same","mu=50");
  sb.SetLineColor(kBlue);sb.SetLineWidth(2);
  
  
  //we could redraw with the mean in a different location like this:
  //sb.Draw("","mu=50,mean=3.1");

  //Define a dataset ... sample events from efunc again
  RooDataSet data("data","data",RooArgSet(x));
  for(int i=0;i<155;i++) {
    //can add data like this:
    x = efunc.GetRandom();
    data.add(x);
  }
  //can visualize datasets by filling a histogram with them ...
  TH1* hData = (TH1*)hBkg->Clone("hData");
  hData->Reset();
  data.fillHistogram(hData,x);
  hData->SetMarkerStyle(20);
  hData->Draw("same");


  //Build the full model (only would have effect if we had a nuisance parameter)
  RooAbsPdf* model = TRooFit::BuildModel(sb,data);
  
  model->Print();
  
  //We could run a fit now, if we like ...
  model->fitTo(data);

  sb.Draw("e3005"); //draws post-fit result
  hData->Draw("same");

  //could move mean and then refit ...
  mean = 7.5;
  model->fitTo(data);

  sb.Draw("e3005"); //draws post-fit result
  hData->Draw("same");
  
  //can save everything to a workspace for limit setting etc ...
   //can shove everything into a workspace now, if we want to save to a file ..
  RooWorkspace w("w","w");
  w.import(*model);
  w.import(data);
  auto mc = TRooFit::CreateModelConfig(w, model->GetName(),data.GetName(),mu.GetName());
  w.import(*mc);
  w.writeToFile("/tmp/myModelWorkspace.root");

}