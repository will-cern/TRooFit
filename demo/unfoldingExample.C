{
  gStyle->SetOptStat(0);

  //==============================================================
  // Maximum-Likelihood Unfolding example
  //  Unfolds some "data" using a folding matrix from "mc"
  //  An optional regularization chi2 penalty term can be included
  //==============================================================

  bool doRegularization(true); //should we add a regularization chi2 penalty term?

  TFile f1("histograms.root");
  //get three different truth and corresponding reco distributions
  //data, mc1, and mc2
  TH1D* hist_data_rec = (TH1D*)f1.Get("hist_data_rec_10bins");
  TH1D* hist_data_truth = (TH1D*)f1.Get("hist_data_truth");
  TH1D* hist_mc1_rec = (TH1D*)f1.Get("hist_mc1_rec_10bins");
  TH1D* hist_mc1_gen = (TH1D*)f1.Get("hist_mc1_gen");
  TH1D* hist_mc2_rec = (TH1D*)f1.Get("hist_mc2_rec_10bins");
  TH1D* hist_mc2_gen = (TH1D*)f1.Get("hist_mc2_gen");
  
  //Get the response matrix for mc1
  TH2D* hist_mc1_recgen = (TH2D*)f1.Get("hist_mc1_recgen_10bins");

  TCanvas myCanvas("unfoldingInputs","TRooFit Unfolding Inputs",900,600);
  myCanvas.Divide(2,2);
  myCanvas.cd(1);
  hist_data_truth->SetLineColor(kBlue);
  hist_data_truth->DrawCopy("hist");
  hist_mc1_gen->SetLineColor(kRed);
  hist_mc1_gen->DrawCopy("hist same");
  myCanvas.cd(2);
  hist_mc1_rec->SetLineColor(kRed);
  hist_mc1_rec->DrawCopy("hist");
  hist_data_rec->SetMarkerStyle(20);
  hist_data_rec->DrawCopy("hist same");
  hist_data_rec->DrawCopy("p same");
  
  myCanvas.cd(3);
  hist_mc1_recgen->Draw("COLZ");


  //switch these to point at whatever you want to use as data
  //and what you want to show at truth level to compare to unfolded
  TH1* theDataHist = hist_data_rec;
  TH1* theTruthHist = hist_data_truth;
//   TH1* theDataHist = hist_mc1_rec;
//   TH1* theTruthHist = hist_mc1_gen;
//   TH1* theDataHist = hist_mc2_rec;
//   TH1* theTruthHist = hist_mc2_gen;


  //TRooFit unfolding demo ...
  //hist_mc1_recgen contains events generated at truth level that made it to reco 
  //hist_mc1_gen contains the actual events generated ... 
  //hist_mc1_recgen events in the last x bin are actually bkg events
  
  //1. Assemble model from these two histograms 
  TRooHStack all("allbins","Model"); //stack of signal
  //need a var to represent the observable
  RooRealVar r("r","r_{rec}",0,1.);
  
  
  for(int ii=1;ii<=hist_mc1_gen->GetNbinsX();ii++) {
    int i = ii-1;if(i==0) i=hist_mc1_gen->GetNbinsX();//do the last bin first, just so 'bkg' is at bottom of stack
    TRooH1D* h = new TRooH1D(Form("bin%d",i),Form("bin%d",i),r,hist_mc1_recgen->GetNbinsY());
    RooRealVar* n = new RooRealVar(Form("norm%d",i),Form("norm%d",i),1,0,100);
    h->addNormFactor(*n);
    
    for(int j=1;j<=hist_mc1_recgen->GetNbinsY();j++) {
      h->SetBinContent(j,hist_mc1_recgen->GetBinContent(i,j));
      h->SetBinError(j,hist_mc1_recgen->GetBinError(i,j));
    }
    //add the missing amount, i.e. events not reconstructed: hist_mc1_gen - integral for this hist 
    h->SetMissingContent( hist_mc1_gen->GetBinContent(i) - h->Integral() );
    
    h->SetFillColor(i);
    
    all.Add(h);
  }
  
  TCanvas c("c","TRooFit Unfolding demo",900,600);c.Divide(2,2);
  
  c.cd(1);
  theDataHist->SetMarkerStyle(20);
  hist_mc1_rec->SetLineColor(kRed);
  hist_mc1_rec->SetLineWidth(2);
  hist_mc1_rec->DrawCopy("histe");
  theDataHist->Draw("histe same");
  
  
  
  c.cd(2);
  all.SetLineColor(kRed);all.SetLineWidth(2);
  all.Draw("e3005");
  

  //2. Assemble histogram to represent unfolded result 
  // Each bin is the integral of a single 'hist' in the stack
  // Must include the 'missing' content ... i.e. lost events 
  
  RooRealVar r_gen("r_gen","r_gen",0,1);
  TRooH1D truth("truth","Truth Distribution",r_gen,hist_mc1_gen->GetNbinsX());
  for(int i=1;i<=hist_mc1_gen->GetNbinsX();i++) {
    TRooH1D* h = (TRooH1D*)all.findServer(Form("bin%d",i)); //locates the ith hist
    truth.SetBinContent(i, *h->createIntegralWM(r)); //can set a bin content equal to a function!
  }
  
  
  

  //3. Assemble roofit dataset
  RooRealVar w("weightVar","weightVar",0,10000);
  RooDataSet data("data","data",RooArgSet(r,w),"weightVar");
  for(int i=1;i<=hist_data_rec->GetNbinsX();i++) {
    r = theDataHist->GetBinCenter(i);
    data.add(r,theDataHist->GetBinContent(i));
  }
  
  auto model = TRooFit::BuildModel(all,data); //this returns a RooProdPdf with additional constraint terms for NP uncerts
  
  
  //3.5 OPTIONAL Add a regularization 'penality term'
  // Penalty term is equivalent to a product of gaussian constraints
  // This is done to reduce oscillations that occur when significant bin migration
  // Penality is (truth - mc_truth)/sigma^2 ... where 1/sigma^2 is the regularization strength
  // when the strength is 0, then there is effectively 0 penalty for pulling the truth from the mc_truth
  //
  // modify the bin error in hist_mc1_gen so that every bin has an error = regularization strength
  //
  //Note that adding a regularization term like this will bias the fit!
  TRooChi2Constraint* penalty =0;
  TRooGPConstraint* gpPenalty = 0;
  if(doRegularization) {
    double tau=pow(10.,-2.25);
    for(int i=1;i<=hist_mc1_gen->GetNbinsX();i++) {
      hist_mc1_gen->SetBinError(i,1./tau);
    }
    
    RooDataHist* truth_ref = new RooDataHist("truth_ref","truth_ref",r_gen,hist_mc1_gen);
    penalty = new TRooChi2Constraint("penalty","penalty",truth,*truth_ref,true);
    
    TMatrixD k(hist_mc1_gen->GetNbinsX(),hist_mc1_gen->GetNbinsX());
    for(int i=0;i<hist_mc1_gen->GetNbinsX();i++) k(i,i)=2./(tau*tau);
    gpPenalty = new TRooGPConstraint("penalty","penalty",truth,*truth_ref,k);
    
    model = new RooProdPdf("allWithPenalty","Model with Regularization Penalty",RooArgList(*model,*penalty));
    
  }

  //4. Run the fit ...
  RooFitResult* res = model->fitTo(data,RooFit::Save());
  
  
  //5. Plot results
  
  //check that everything looks good, post-fit
  c.cd(3);
  all.SetLineColor(kBlack);all.SetLineWidth(1);
  all.Draw("e3005");
  all.SetLineColor(kRed);all.SetLineWidth(2);
  all.Draw("e3004same init hist",res);
  theDataHist->Draw("same");
  
  //draw the original truth hist (theTruthHist) and the final unfolded
  //distribution (truth)
  c.cd(4);
  theTruthHist->SetLineColor(kBlue);
  theTruthHist->Draw("hist"); //compare to truth
  
  //also draw the truth_ref 
  hist_mc1_gen->SetLineColor(kRed);
  hist_mc1_gen->DrawCopy("hist same");
  hist_mc1_gen->SetFillStyle(3005);
  hist_mc1_gen->SetFillColor(kRed);
  hist_mc1_gen->Draw("e2same");
  
  truth.SetMarkerStyle(20);
  truth.Draw("same");
  
  
  //draw correlation matrix on another canvas
  TCanvas corCanvas;
  RooArgSet allArgs; model->treeNodeServerList(&allArgs);
  RooArgList normFactors;
  for(int i=1;i<=10;i++) {
    normFactors.add( *allArgs.find(Form("norm%d",i)) );
  }
  TRooFitResult rr(res);
  rr.Draw("cor",normFactors);
  

}