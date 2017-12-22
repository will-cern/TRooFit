

#include "TRooFit/apps/TRooABCD.h"

#include "TRooFit/Utils.h"

#include "TLatex.h"
#include "TCanvas.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

TRooABCD::TRooABCD(const char* name, const char* title, double xlow, double xhigh) : TNamed(name,title) {
  m_xVar = new RooRealVar("x","x",xlow,xhigh);
  m_cat = new RooCategory("region","region");
  const char* regionLabels[4] = {"A","B","C","D"};
  for(int i=0;i<4;i++) m_cat->defineType(regionLabels[i]);
  m_weightVar = new RooRealVar("w","w",1);
  
  m_data = new RooDataSet("data","data",RooArgSet(*m_xVar,*m_cat,*m_weightVar),"w");
  
  m_signalStrength = new RooRealVar("mu","#mu",1,0,10);
  
}


Bool_t TRooABCD::AddData(int region, TH1* data) {
  m_cat->setIndex(region);
  
  for(int i=1;i<=data->GetNbinsX();i++) {
    *m_xVar = data->GetBinCenter(i);
    if(data->GetBinContent(i)) {
      double remain;
      if(std::modf(data->GetBinContent(i),&remain)==0) {
        for(int j=1;j<int(data->GetBinContent(i)+0.5);j++) m_data->add( RooArgSet(*m_xVar,*m_cat) ); //adds each event individually
      } else {
        m_data->add( RooArgSet(*m_xVar,*m_cat), data->GetBinContent(i)/*, data->GetBinError(i)*/ ); //errors don't store properly in RooDataSet :-(
     }
    }
  }
  
  if(m_bkg[region]==0) {
    m_bkg[region] = (data->GetXaxis()->IsVariableBinSize()) ? 
                                      new TRooH1D(Form("bkd%d",region),
                                       Form("Background in region %d",region),*m_xVar,
                                       data->GetNbinsX(),
                                       data->GetXaxis()->GetXbins()->GetArray())
                                  :
                                      new TRooH1D(Form("bkg%d",region),
                                       Form("Background in region %d",region),*m_xVar,
                                       data->GetNbinsX(),
                                       data->GetXaxis()->GetXmin(),
                                       data->GetXaxis()->GetXmax());
    
    m_bkg[region]->SetFillColor(kCyan);
    m_bkg[region]->SetLineColor(kBlue);
    m_bkg[region]->setFloor(true); //says that the value of this pdf histogram cannot be less than 0
    m_bkg[region]->SetMinimum(1e-9);//y-axis minimum ... for when we draw the histograms later
    
    //if region is B (1) or D (3) then add scale factors 
    const char* regionLabels[4] = {"A","B","C","D"};
    if(region==1||region==3) {
      for(int i=1;i<=data->GetNbinsX();i++) {
        RooRealVar* sf = new RooRealVar(Form("sf_%s_bin%d",regionLabels[region],i),
                                        Form("sf_{%d}^{%s}",i,regionLabels[region]),1,-1e-9,5); //scale factor range is between just under 0 and 5
        m_bkg[region]->addShapeFactor(i,*sf);
      }
      
      //we can also create the A and B region plots too, now ...
      if(m_bkg[region-1]==0) {
        m_bkg[region-1] = (data->GetXaxis()->IsVariableBinSize()) ? 
                                      new TRooH1D(Form("bkd%d",region-1),
                                       Form("Background in region %d",region-1),*m_xVar,
                                       data->GetNbinsX(),
                                       data->GetXaxis()->GetXbins()->GetArray())
                                  :
                                      new TRooH1D(Form("bkg%d",region-1),
                                       Form("Background in region %d",region-1),*m_xVar,
                                       data->GetNbinsX(),
                                       data->GetXaxis()->GetXmin(),
                                       data->GetXaxis()->GetXmax());
    
        m_bkg[region-1]->SetFillColor(kCyan);
        m_bkg[region-1]->SetLineColor(kBlue);
        m_bkg[region-1]->setFloor(true); //says that the value of this pdf histogram cannot be less than 0
        m_bkg[region-1]->SetMinimum(1e-9);//y-axis minimum ... for when we draw the histograms later
        m_bkg[region-1]->Fill( *m_bkg[region] );
      }
      
    } else {
      //need to Fill this component with the control region version
      if(m_bkg[region+1]==0) {
        //not existing yet ... so create it ...
        TH1* emptyHist = (TH1*)data->Clone("emptyHist");
        emptyHist->Reset();
        AddData(region+1,emptyHist);
        delete emptyHist;
      }
      
      m_bkg[region]->Fill( *m_bkg[region+1] );
      
    }
    
  }
  
  if(region==1||region==3) {
    //update the nominal value to match the data ...
    
    for(int i=1;i<=data->GetNbinsX();i++) {
      m_bkg[region]->SetBinContent(i,m_bkg[region]->GetBinContent(i)+data->GetBinContent(i));
      if(m_bkg[region]->GetBinContent(i)<2) m_bkg[region]->SetBinContent(i,2);
    }
    
  }
  
  //refresh the dataHist
  if(m_dataHist[region]==0) {
    m_dataHist[region] = (TH1*)data->Clone(Form("data%d",region));
    m_dataHist[region]->SetDirectory(0); m_dataHist[region]->Sumw2();
    m_dataHist[region]->SetMarkerStyle(20);
  }
  
  m_dataHist[region]->Reset();
  m_data->fillHistogram( m_dataHist[region], *m_xVar, Form("region==%d",region ) );
  
  

  return kTRUE;
}

 
const TH1* TRooABCD::GetDataHistogram(int region) {
  return m_dataHist[region];
}
 
Double_t TRooABCD::GetBkgBinContent(int region, int bin) {
  if(m_bkg[region]==0) return 0;
  return m_bkg[region]->GetBinContent(bin);
}

Double_t TRooABCD::GetBkgIntegral(int region, TRooFitResult* r) {
  if(m_bkg[region]==0) return 0;
  return m_bkg[region]->Integral("",(r==0)?m_lastFitResult:r);
}

Double_t TRooABCD::GetBkgIntegralError(int region, TRooFitResult* r) {
  if(m_bkg[region]==0) return 0;
  double out;
  m_bkg[region]->IntegralAndError(out,(r==0)?m_lastFitResult:r);
  return out;
}

Double_t TRooABCD::GetSignalIntegral(int region, TRooFitResult* r) {
  if(m_signal[region]==0) return 0;
  return m_signal[region]->Integral("",(r==0)?m_lastFitResult:r);
}

Double_t TRooABCD::GetSignalIntegralError(int region, TRooFitResult* r) {
  if(m_signal[region]==0) return 0;
  double out;
  m_signal[region]->IntegralAndError(out,(r==0)?m_lastFitResult:r);
  return out;
}

Double_t TRooABCD::GetDataBinContent(int region, int bin) {
  if(m_dataHist[region]==0) return 0;
  return m_dataHist[region]->GetBinContent(bin);
}

Bool_t TRooABCD::AddOther(int region, TH1* other) {
  if(m_other[region]==0) {
    m_other[region] = (other->GetXaxis()->IsVariableBinSize()) ? 
                                      new TRooH1D(Form("signal%d",region),
                                       Form("Signal in region %d",region),*m_xVar,
                                       other->GetNbinsX(),
                                       other->GetXaxis()->GetXbins()->GetArray())
                                  :
                                      new TRooH1D(Form("signal%d",region),
                                       Form("Signal in region %d",region),*m_xVar,
                                       other->GetNbinsX(),
                                       other->GetXaxis()->GetXmin(),
                                       other->GetXaxis()->GetXmax());
    m_other[region]->SetFillColor(kOrange);
    m_other[region]->SetLineColor(kOrange);
    m_other[region]->setFloor(true); //says that the value of this pdf histogram cannot be less than 0
    m_other[region]->SetMinimum(1e-9);//y-axis minimum ... for when we draw the histograms later
  }
  
  m_other[region]->Add(other);


  return kTRUE;
}

Bool_t TRooABCD::AddSignal(int region, TH1* signal) {
  if(m_signal[region]==0) {
    m_signal[region] = (signal->GetXaxis()->IsVariableBinSize()) ? 
                                      new TRooH1D(Form("signal%d",region),
                                       Form("Signal in region %d",region),*m_xVar,
                                       signal->GetNbinsX(),
                                       signal->GetXaxis()->GetXbins()->GetArray())
                                  :
                                      new TRooH1D(Form("signal%d",region),
                                       Form("Signal in region %d",region),*m_xVar,
                                       signal->GetNbinsX(),
                                       signal->GetXaxis()->GetXmin(),
                                       signal->GetXaxis()->GetXmax());
    m_signal[region]->SetFillColor(kRed);
    m_signal[region]->SetLineColor(kRed);
    m_signal[region]->setFloor(true); //says that the value of this pdf histogram cannot be less than 0
    m_signal[region]->SetMinimum(1e-9);//y-axis minimum ... for when we draw the histograms later
    m_signal[region]->addNormFactor(*m_signalStrength);
  }
  
  m_signal[region]->Add(signal);

  return kTRUE;

}



RooRealVar* TRooABCD::AddBkgScaleFactor(int region, double value, double uncert) {
  //adds a scale factor with uncertainty to the bkg term 
  if(m_bkg[region]==0) {
    Error("AddBkgScaleFactor","Please add data to region %d before adding a scale factor",region);
    return 0;
  }
  
  int i=1;
  while( m_bkg[region]->findServer(Form("bkg_sf%d_region%d",i,region)) ) {
    i++;
  }
  
  RooRealVar* sf = new RooRealVar(Form("bkg_sf%d_region%d",i,region),Form("#alpha_{%d}^{%s}",i,m_cat->lookupType(i)->GetName()),value,value-5*uncert,value+5*uncert);
  if(uncert) sf->setStringAttribute("constraintType",Form("GAUSSIAN(%f,%f)",value,uncert));
  else sf->setConstant(1);
  
  m_bkg[region]->addNormFactor( *sf );
  
  return sf;
  
}

RooRealVar* TRooABCD::AddSignalScaleFactor(int region, double value, double uncert) {
  //adds a scale factor with uncertainty to the bkg term 
  if(m_signal[region]==0) {
    Error("AddSignalScaleFactor","Please add signal to region %d before adding a scale factor",region);
    return 0;
  }
  
  int i=1;
  while( m_signal[region]->findServer(Form("sig_sf%d_region%d",i,region)) ) {
    i++;
  }
  
  RooRealVar* sf = new RooRealVar(Form("sig_sf%d_region%d",i,region),Form("#alpha_{%d}^{%s}",i,m_cat->lookupType(i)->GetName()),value,value-5*uncert,value+5*uncert);
  if(uncert) sf->setStringAttribute("constraintType",Form("GAUSSIAN(%f,%f)",value,uncert));
  else sf->setConstant(1);
  
  m_signal[region]->addNormFactor( *sf );
  
  return sf;
  
}


TRooFitResult* TRooABCD::Fit(int modelType, bool floatSignal) {
  if(m_dataHist[1]==0) {
    Error("Fit","No data provided for region B. You must at least construct at empty histogram and pass to AddData");
    return 0;
  }
  if(m_dataHist[2]==0) {
    Error("Fit","No data provided for region C. You must at least construct at empty histogram and pass to AddData");
    return 0;
  }
  if(m_dataHist[3]==0) {
    Error("Fit","No data provided for region D. You must at least construct at empty histogram and pass to AddData");
    return 0;
  }
  
  bool hasSignal = false;
  for(int i=0;i<4;i++) if(m_signal[i]) hasSignal=true;
  
  if(hasSignal && m_dataHist[0]==0 && floatSignal && m_dataHist[1]->GetNbinsX()==1 && m_dataHist[3]->GetNbinsX()==1) {
    Error("Fit","Cannot float signal when signal region is blinded and only have 1 bin in the control regions");
    return 0;
  }
  
  m_signalStrength->setConstant(!floatSignal);
  
  


  //FIXME: need to remove these from bkg terms before deleting
  for(auto t : m_transferFactors) {
    for(int i=0;i<4;i++) m_bkg[i]->removeNormFactor( *t );
  
    delete t; 
  }
  m_transferFactors.clear(); 


  if(modelType==0) {
    //doing flat model     
    //A = mB
    //C = mD
    
    //estimate parameter from C over D
    double parEstimate = m_dataHist[2]->Integral()/m_dataHist[3]->Integral();
    
    if(m_modelPars.size()<1) {
      m_modelPars.push_back(new RooRealVar("m1","#tilde{m}_{1}",parEstimate,0,parEstimate*100));
    } else {
      m_modelPars[0]->setRange(0,parEstimate*100);
      m_modelPars[0]->setVal(parEstimate);
    }
    
    
  
    m_transferFactors.push_back( new TRooHF1D("transfer","Transfer Factor",*m_xVar) );
    m_transferFactors[0]->SetLineColor(kBlue);
    m_transferFactors[0]->Fill( *m_modelPars[0] );
    
    m_bkg[2]->addNormFactor( *m_transferFactors[0] );
    m_bkg[0]->addNormFactor( *m_transferFactors[0] );
  } else if(modelType==1) {
    if( m_dataHist[1]->GetNbinsX()==1 && m_dataHist[3]->GetNbinsX()==1 ) {
      Error("Fit","Cannot use linear model when you only have 1 bin per region");
      return 0;
    }
  
    //doing linear model
    //A = (m1 + m2x)B
    //C = (m1 + m2x)D
    
    
    //estimate parameter from C over D
    double parEstimate = m_dataHist[2]->Integral()/m_dataHist[3]->Integral();
    
    if(m_modelPars.size()<1) {
      m_modelPars.push_back(new RooRealVar("m1","#tilde{m}_{1}",parEstimate,-5,parEstimate*100));
    } else {
      m_modelPars[0]->setRange(-5,parEstimate*100);
      m_modelPars[0]->setVal(parEstimate);
    }
    
    if(m_modelPars.size()<2) {
      m_modelPars.push_back(new RooRealVar("m2","#tilde{m}_{2}",0,-5,5));
    }
    
    RooFormulaVar* m_linearFormula = new RooFormulaVar("transferFunc","Transfer Factor as function of x","(m2*x + m1)",RooArgList(*m_modelPars[0],*m_modelPars[1],*m_xVar));
    
    m_transferFactors.push_back( new TRooHF1D("transfer1","Transfer Factor D->C",*m_xVar) );
    m_transferFactors[0]->SetLineColor(kBlue);
    m_transferFactors[0]->Fill( *m_linearFormula );
    
    m_transferFactors[0]->setFloor(true);
    
    m_bkg[2]->addNormFactor( *m_transferFactors[0] );
    m_bkg[0]->addNormFactor( *m_transferFactors[0] );
    
  } 
  
  
  //assemble stacks
  for(auto t : m_stacks) delete t.second; 
  m_stacks.clear(); 
  
  for(int i=0;i<4;i++) {
    m_stacks[i] = new TRooHStack(Form("stack%d",i),Form("Model in region %d",i));
    m_stacks[i]->SetMinimum(1e-9);
    m_stacks[i]->SetStats(0);
    if(m_other[i]) m_stacks[i]->Add(m_other[i],false); //do not acquire stat factors
    m_stacks[i]->Add(m_bkg[i]);
    if(m_signal[i]) m_stacks[i]->Add(m_signal[i],false);
  }
  
  TCanvas* cc = new TCanvas(Form("%s_prefit",GetName()),Form("%s Pre-fit",GetTitle()),800,600);;
  cc->Divide(2,2);
  cc->cd(1);if(m_lastFitResult) m_stacks[2]->Draw("e3005",m_lastFitResult); else  m_stacks[2]->Draw("e3005"); m_dataHist[2]->Draw("same");
  cc->cd(2);if(m_lastFitResult) m_stacks[0]->Draw("e3005",m_lastFitResult); else  m_stacks[0]->Draw("e3005");if(m_dataHist[0]) m_dataHist[0]->Draw("same");
  cc->cd(3);if(m_lastFitResult) m_stacks[3]->Draw("e3005",m_lastFitResult); else  m_stacks[3]->Draw("e3005");m_dataHist[3]->Draw("same");
  cc->cd(4);if(m_lastFitResult) m_stacks[1]->Draw("e3005",m_lastFitResult); else  m_stacks[1]->Draw("e3005");m_dataHist[1]->Draw("same");
  
  cc->Modified(1);
  cc->Update();
  
  //put stacks into a model
  if(m_model) delete m_model;
  
  m_model = new RooSimultaneous("model","model",*m_cat);
  for(int i=0;i<4;i++) {
    auto modelWithConstraints = TRooFit::BuildModel(*m_stacks[i],*m_data);
    m_model->addPdf( *modelWithConstraints , m_cat->lookupType(i)->GetName() );
  }
  
  ///now run the fit
  m_lastFitResult = 0;
  
  RooMsgService::instance().setGlobalKillBelow(m_printLevel);
  if(m_printLevel>=RooFit::ERROR) {
    Info("Fit","Running Fit ... (use SetPrintLevel(RooFit::INFO) to show roofit output)...");
  }
  
  RooAbsData* theData = m_data;
  if(m_dataHist[0]==0) {
    m_xVar->setRange("myRange",m_dataHist[1]->GetXaxis()->GetXmin(),m_dataHist[1]->GetXaxis()->GetXmax());
    m_stacks[0]->setBlindRange("myRange"); //the blind range will make the value of this component appear to be 0 in this range
    theData = m_data->reduce("region!=0"); //removes the region A data
  }
  
  double prediction_mu0 = 0;
  
  if(floatSignal && hasSignal) {
    //first run the fit with a signal strength of 0
    m_signalStrength->setVal(0);
    m_signalStrength->setConstant(1);
    
    
    Info("Fit","Running mu=0 Fit (to compute signal strength systematic) ... ");
    m_lastFitResult = new TRooFitResult(m_model->fitTo(*theData,RooFit::Save()));
    
    //obtain the estimate.. must first remove blinding
    if(m_dataHist[0]==0) m_stacks[0]->setBlindRange(""); 
    prediction_mu0 = GetBkgIntegral(0);
    if(m_dataHist[0]==0) m_stacks[0]->setBlindRange("myRange"); //put it back
    
    Info("Fit","Background Signal region prediction = %g ... ",prediction_mu0);
    
    delete m_lastFitResult;
    
    m_signalStrength->setVal(1.); //for getting the bin content
    
    double totSignal = 0;
    double totData = 0;
    //scale signal up so that yield in control regions equals data in the control regions
    for(int i=0;i<4;i++) {
      if(!m_signal[i]) continue;
      if(!m_dataHist[i]) continue;
      
      for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
        double t = m_signal[i]->GetBinContent(j);
        if(!t) continue;
        totSignal += t;
        totData += m_dataHist[i]->GetBinContent(j);
      }
      
    }
    
    Info("Fit","Setting signal strength to %g to test for maximal signal", (totData/totSignal));
    
    m_signalStrength->setConstant(0);
    m_signalStrength->setRange(0,(totData/totSignal)+1);
    m_signalStrength->setVal(totData/totSignal);
    
    Info("Fit","Re-running fit with floating signal...");
    
  }
  
  
  m_lastFitResult = new TRooFitResult(m_model->fitTo(*theData,RooFit::Save()));
  
  if(m_dataHist[0]==0) {
    m_stacks[0]->setBlindRange(""); //remove the blinding range for the post-fit plotting
    delete theData;
  }
  
  
  
  ///draw the post-fit results
  
  TCanvas* cc2 = new TCanvas(Form("%s_postfit",GetName()),Form("%s Post-fit",GetTitle()),800,750);
  cc2->Divide(2,3);
  
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    double maxData = m_dataHist[i]->GetBinContent( m_dataHist[i]->GetMaximumBin() );
    if(maxData) m_stacks[i]->SetMaximum( maxData*1.1 ); //ensures data is in axis range
  }
  
  cc2->cd(1);m_stacks[2]->Draw("e3005",m_lastFitResult);m_dataHist[2]->Draw("same");
  cc2->cd(2);m_stacks[0]->Draw("e3005",m_lastFitResult);if(m_dataHist[0]) m_dataHist[0]->Draw("same");
  cc2->cd(3);m_stacks[3]->Draw("e3005",m_lastFitResult);m_dataHist[3]->Draw("same");
  cc2->cd(4);m_stacks[1]->Draw("e3005",m_lastFitResult);m_dataHist[1]->Draw("same");
  
  TH1* data_ratio = (TH1*)m_dataHist[2]->Clone("data_ratio"); data_ratio->SetTitle("Transfer factor");
  data_ratio->Divide(m_dataHist[3]);
  data_ratio->SetStats(0);
  cc2->cd(5);data_ratio->Draw();m_transferFactors[0]->Draw("val L same e3005",m_lastFitResult);

  cc2->cd(6);
  TLatex t;
  double prediction = GetBkgIntegral(0);
  double predictionError = GetBkgIntegralError(0);
  if(floatSignal && hasSignal) {
    t.DrawLatex(0.05,0.8,Form("SR Bkg Predicted = %g #pm %g (syst.) #pm %g (mu syst.)",prediction,predictionError,fabs(prediction-prediction_mu0)));
  } else {
    t.DrawLatex(0.05,0.8,Form("SR Bkg Predicted = %g #pm %g",prediction,predictionError));
  }
  if(m_dataHist[0]) {
    double signif = RooStats::NumberCountingUtils::BinomialObsZ(m_dataHist[0]->Integral(),prediction,predictionError/prediction);
    
    if(signif < 1) t.SetTextColor(kGreen);
    else if(signif < 2.5) t.SetTextColor(kOrange);
    else t.SetTextColor(kRed);
    
    t.DrawLatex(0.05,0.75,Form("SR Observed = %g (%.1f#sigma)",m_dataHist[0]->Integral(),signif));
  }
  else t.DrawLatex(0.05,0.75,Form("SR Observed = BLINDED"));
  
  //compute chi^2 for the available data .
  
  double chi2 = 0;
  int ndof=0;
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    TH1* model = m_stacks[i]->GetHistogram(m_lastFitResult,true);
    for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
      if(m_dataHist[i]->GetBinContent(j) != model->GetBinContent(j)) {
        double val = pow( (m_dataHist[i]->GetBinContent(j) - model->GetBinContent(j)) , 2 )/( pow(model->GetBinError(j),2)+model->GetBinContent(j) ); //using error = model poisson error + model error (in quadrature)
        if(m_printLevel<=RooFit::INFO) Info("Fit","chi2 for bin %d in region %s = %g",j,m_cat->lookupType(i)->GetName(),val);
        chi2 += val;
      }
      ndof++;
    }
    //pVal *= m_dataHist[i]->Chi2Test( model, "UW" );
    delete model;
  }
  
  //correct ndof to account for free parameters
  if(floatSignal && hasSignal) ndof--;
  ndof -= m_dataHist[1]->GetNbinsX();
  ndof -= m_dataHist[3]->GetNbinsX();
  if(modelType==0) ndof--;
  else if(modelType==1) ndof -= 2;
  
  if(ndof>0) {
    double pVal = TMath::Prob(chi2,ndof);
    if(RooStats::PValueToSignificance(pVal) < 1) t.SetTextColor(kGreen);
    else if(RooStats::PValueToSignificance(pVal) < 2.5) t.SetTextColor(kOrange);
    else t.SetTextColor(kRed);
    t.DrawLatex(0.05,0.65,Form("Fit #chi^{2} p-value = %g (%.1f#sigma)",pVal,RooStats::PValueToSignificance(pVal)));
  }
  
  
  //compute signal contributions in unblinded regions ... any bigger than ??? just recommend
  //that user rerun Fit with fixed signal strength = 0 (SetSignalStrength(0);Fit(x,false) and compare bkg predictions
  
  

  return m_lastFitResult;
}

