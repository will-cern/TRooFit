

#include "TRooFit/apps/TRooABCD.h"

#include "TRooFit/Utils.h"

#include "TLatex.h"
#include "TCanvas.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

TRooABCD::TRooABCD(const char* name, const char* title) : TNamed(name,title) {
  m_xVar = new RooRealVar("x","x",0,1);m_xVar->setAttribute("isInit",true);
  m_cat = new RooCategory("region","region");
  const char* regionLabels[4] = {"A","B","C","D"};
  for(int i=0;i<4;i++) m_cat->defineType(regionLabels[i]);
  m_weightVar = new RooRealVar("w","w",1);
  
  m_data = new RooDataSet("data","data",RooArgSet(*m_xVar,*m_cat,*m_weightVar),"w");
  
  m_signalStrength = new RooRealVar("mu","#mu",1,0,10);
  m_signalStrength->setStringAttribute("fixedValue","0"); //whenever this parameter is floating, we will first run the fit with it fixed at this value, and quote diff
  m_signalStrength->setConstant(true); //becomes floating as soon as any signal is added (see AddSignal)
  
  //create two model parameters
  m_modelPars.push_back(new RooRealVar("m0","#tilde{m}_{0}",0,-5,5));
  m_modelPars.push_back(new RooRealVar("m1","#tilde{m}_{1}",0,-5,5)); //fixme ... should we change this range?
  
  //by default, second is constant, so that we are doing the standard ABCD method
  m_modelPars[1]->setConstant(true);
  

  m_allParameters.add(*m_signalStrength);
  m_allParameters.add(*m_modelPars[0]);
  m_allParameters.add(*m_modelPars[1]);
  
}

void TRooABCD::checkRangeChange(TH1* hist) {
  bool changedXRange(false);
  if(m_xVar->getAttribute("isInit")) {
    //first hist load 
    changedXRange=true;
    m_xVar->setRange(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    m_xVar->setAttribute("isInit",false);
  }
  if(hist->GetXaxis()->GetXmin() < m_xVar->getMin()) {
    m_xVar->setRange(hist->GetXaxis()->GetXmin(), m_xVar->getMax());
    changedXRange=true;
  }
  if(hist->GetXaxis()->GetXmax() > m_xVar->getMax()) {
    m_xVar->setRange(m_xVar->getMin(),hist->GetXaxis()->GetXmax());
    changedXRange=true;
  }
  
  if(changedXRange) {
    //have to recreate the dataset, since otherwise RooFit will complain ranges dont match
    RooDataSet* newData = new RooDataSet("data","data",m_data,RooArgSet(*m_xVar,*m_cat,*m_weightVar),0,"w");
    delete m_data;
    m_data = newData;
  }
}

Bool_t TRooABCD::AddAsimovData(int region, double signalStrength) {
  if(region!=0) {
    Error("AddAsimovData","Can only add to region 0 at the moment");
    return false;
  }
  m_useAsimov=true;
  m_asimovSignalStrength=signalStrength;
  return true;
}

Bool_t TRooABCD::addAsimovData(int region, TH1* data) {
  m_cat->setIndex(region);
  
  checkRangeChange(data);
  
  if(!m_asimovData) {
    m_asimovData = new RooDataSet("asimovdata","asimovdata",RooArgSet(*m_xVar,*m_cat,*m_weightVar),"w");
  }
  
  for(int i=1;i<=data->GetNbinsX();i++) {
    *m_xVar = data->GetBinCenter(i);
    if(data->GetBinContent(i)) {
      double remain;
      if(std::modf(data->GetBinContent(i),&remain)==0) {
        for(int j=0;j<int(data->GetBinContent(i)+0.5);j++) m_asimovData->add( RooArgSet(*m_xVar,*m_cat) ); //adds each event individually
      } else {
        m_asimovData->add( RooArgSet(*m_xVar,*m_cat), data->GetBinContent(i)/*, data->GetBinError(i)*/ ); //errors don't store properly in RooDataSet :-(
     }
    }
  }
  
  
  
  //refresh the dataHist
  if(m_asimovDataHist[region]==0) {
    m_asimovDataHist[region] = (TH1*)data->Clone(Form("asimovdata%d",region));
    m_asimovDataHist[region]->SetDirectory(0); m_asimovDataHist[region]->Sumw2();
    m_asimovDataHist[region]->SetMarkerStyle(24);m_asimovDataHist[region]->SetLineStyle(2);
  }
  
  m_asimovDataHist[region]->Reset();
  m_asimovData->fillHistogram( m_asimovDataHist[region], *m_xVar, Form("region==%d",region ) );
  
  //histograms filled from RooDataSet dont retain errors properly
  for(int i=1;i<=m_asimovDataHist[region]->GetNbinsX();i++) m_asimovDataHist[region]->SetBinError(i, sqrt( m_asimovDataHist[region]->GetBinContent(i) ) );
  

  return kTRUE;
}

Bool_t TRooABCD::AddValidationData(int region, TH1* data) {
  if(region!=0) {
    Error("AddValidationData","Validation data can only be added to the signal region at this time");
    return false;
  }
  if(m_dataHist[0] && !m_isVR) {
    Error("AddValidationData","Regular data was already added to region %d so cannot add validation data",region);
    return false;
  }
  bool out = AddData(region,data);
  SetValidationMode(true);
  return out;
}

Bool_t TRooABCD::AddData(int region, TH1* data) {
  m_cat->setIndex(region);
  
  checkRangeChange(data);
  
  for(int i=1;i<=data->GetNbinsX();i++) {
    *m_xVar = data->GetBinCenter(i);
    if(data->GetBinContent(i)) {
      /*double remain;
      if(std::modf(data->GetBinContent(i),&remain)==0) {
        for(int j=0;j<int(data->GetBinContent(i)+0.5);j++) m_data->add( RooArgSet(*m_xVar,*m_cat) ); //adds each event individually
      } else {*/
        m_data->add( RooArgSet(*m_xVar,*m_cat), data->GetBinContent(i)/*, data->GetBinError(i)*/ ); //errors don't store properly in RooDataSet :-(
     //}
    }
  }
  
  if(m_bkg[region]==0) {
    m_bkg[region] = (data->GetXaxis()->IsVariableBinSize()) ? 
                                      new TRooH1D(Form("bkg%d",region),
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
        m_allParameters.add(*sf);
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
        m_bkg[region-1]->Add( *m_bkg[region] );
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
      
      m_bkg[region]->Add( *m_bkg[region+1] );
      
    }
    
  }
  
  if(region==1||region==3) {
    //update the nominal value to match the data ...
    
    for(int i=1;i<=data->GetNbinsX();i++) {
      m_bkg[region]->SetBinContent(i,/*m_bkg[region]->GetBinContent(i)+*/data->GetBinContent(i));
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
  //histograms filled from RooDataSet dont retain errors properly
  for(int i=1;i<=m_dataHist[region]->GetNbinsX();i++) m_dataHist[region]->SetBinError(i, sqrt( m_dataHist[region]->GetBinContent(i) ) );
  
  //change the color of the 'signal region' data histogram if we are
  //in validation mode
  if(region==0) {
    m_dataHist[0]->SetLineColor((m_isVR)?kBlue:kBlack);
    m_dataHist[0]->SetMarkerColor((m_isVR)?kBlue:kBlack);
  }

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
    checkRangeChange(other);
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
    checkRangeChange(signal);
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
  
  if(signal->Integral() && m_signalStrength->isConstant()) m_signalStrength->setConstant(false);
  
  m_signal[region]->Add(signal);

  return kTRUE;

}

RooRealVar* TRooABCD::AddTransferFactorUncertainty(double uncert) {
  //adds a uncertainty to the transfer factor
  if(m_transferFactor2==0) {
    BuildModel();
  }
  
  int i=1;
  while( m_transferFactor2->findServer(Form("tf_sf%d",i)) ) {
    i++;
  }
  
  RooRealVar* sf = new RooRealVar(Form("tf_sf%d",i),Form("#alpha_{%d}^{%s}",i,"TF"),1.,0.,1.+5*uncert);
  if(uncert) sf->setStringAttribute("constraintType",Form("GAUSSIAN(%f,%f)",1.,uncert));
  else sf->setConstant(1);
  
  m_transferFactor2->addNormFactor( *sf );
  m_allParameters.add(*sf);
  
  return sf;
  
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
  m_allParameters.add(*sf);
  
  return sf;
  
}

void  TRooABCD::AddBkgScaleFactor(int region, RooRealVar* sf) {
  if(m_bkg[region]==0) {
    Error("AddBkgScaleFactor","Please add data to region %d before adding a scale factor",region);
    return;
  }
  m_bkg[region]->addNormFactor( *sf );
  m_allParameters.add(*sf);
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
  m_allParameters.add(*sf);
  return sf;
  
}

void  TRooABCD::AddSignalScaleFactor(int region, RooRealVar* sf) {
  if(m_signal[region]==0) {
    Error("AddSignalScaleFactor","Please add signal to region %d before adding a scale factor",region);
    return;
  }
  m_signal[region]->addNormFactor( *sf );
  m_allParameters.add(*sf);
}

RooRealVar* TRooABCD::AddOtherScaleFactor(int region, double value, double uncert) {
  //adds a scale factor with uncertainty to the bkg term 
  if(m_other[region]==0) {
    Error("AddOtherScaleFactor","Please add other to region %d before adding a scale factor",region);
    return 0;
  }
  
  int i=1;
  while( m_other[region]->findServer(Form("other_sf%d_region%d",i,region)) ) {
    i++;
  }
  
  RooRealVar* sf = new RooRealVar(Form("other_sf%d_region%d",i,region),Form("#alpha_{%d}^{%s}",i,m_cat->lookupType(i)->GetName()),value,value-5*uncert,value+5*uncert);
  if(uncert) sf->setStringAttribute("constraintType",Form("GAUSSIAN(%f,%f)",value,uncert));
  else sf->setConstant(1);
  
  m_other[region]->addNormFactor( *sf );
  m_allParameters.add(*sf);
  return sf;
  
}

void  TRooABCD::AddOtherScaleFactor(int region, RooRealVar* sf) {
  if(m_other[region]==0) {
    Error("AddOtherScaleFactor","Please add other to region %d before adding a scale factor",region);
    return;
  }
  m_other[region]->addNormFactor( *sf );
  m_allParameters.add(*sf);
}

bool TRooABCD::BuildModel() {
  if(m_model) {
    Error("BuildModel","Model already built");
    return false;
  }

  RooFormulaVar* m_linearFormula = new RooFormulaVar("transferFunc","Transfer Factor as function of x","(m1*x + m0)",RooArgList(*m_modelPars[0],*m_modelPars[1],*m_xVar));
  m_transferFactor2 = new TRooHF1D("transfer2","Transfer factor",*m_xVar); m_transferFactor2->Add( *m_linearFormula );
  //m_transferFactor->setFloor(true); //can't go negative, would be unphysical
  m_bkg[2]->addNormFactor( *m_transferFactor2 );
  m_bkg[0]->addNormFactor( *m_transferFactor2 );
  
  RooFormulaVar* transFactor = new RooFormulaVar("ratio","Transfer Factor as function of x","(@0+@1)/(@2+@3)",RooArgList(*m_bkg[0],*m_bkg[2],*m_bkg[1],*m_bkg[3]));
  m_transferFactor = new TRooHF1D("transfer1","Transfer Factor D->C",*m_xVar);
  m_transferFactor->SetLineColor(kBlue);
  m_transferFactor->Add( *transFactor );
  
  //assemble stacks
  for(int i=0;i<4;i++) {
    m_stacks[i] = new TRooHStack(Form("stack%d",i),Form("Model in region %d",i));
    m_stacks[i]->SetMinimum(1e-9);
    m_stacks[i]->SetStats(0);
    if(m_other[i]) m_stacks[i]->Add(m_other[i],false); //do not acquire stat factors
    m_stacks[i]->Add(m_bkg[i]);
    if(m_signal[i]) m_stacks[i]->Add(m_signal[i],false);
  }
  //put stacks into a model
  m_model = new RooSimultaneous("model","model",*m_cat);
  for(int i=0;i<4;i++) {
    auto modelWithConstraints = TRooFit::BuildModel(*m_stacks[i],*m_data);
    m_model->addPdf( *modelWithConstraints , m_cat->lookupType(i)->GetName() );
  }
  
  //initialize transfer factor as flat with 0 gradient
  double parEstimate = m_dataHist[2]->Integral()/m_dataHist[3]->Integral();
  if(parEstimate==0) parEstimate=1;
  m_modelPars[0]->setRange(0,parEstimate*100);
  m_modelPars[0]->setVal(parEstimate);
  
  return true;
  
}



TRooFitResult* TRooABCD::Fit(bool drawPostFit) {
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
  
  //if(hasSignal && (m_dataHist[0]==0||m_isVR) && !m_signalStrength->isConstant()) {
  if( (m_useAsimov) && !m_fitForAsimov ) { //if validating or blinded ...
    //when we are blinded and have floating signal, we will run a fit with signal strength = 0
    //i.e. we will obtain the expected data under the bkg-only hypothesis (i.e. asimov data)
    //then we will run the full fit using the asimov data
    
    double val = m_signalStrength->getVal();
    bool cstate = m_signalStrength->isConstant();
    m_signalStrength->setConstant(true);m_signalStrength->setVal(m_asimovSignalStrength);
    m_fitForAsimov=true;
    TRooFitResult* res = Fit(false);
    TH1* pseudo_data = m_stacks[0]->GetHistogram(res);
    addAsimovData(0,pseudo_data);
    auto args = m_model->getObservables(res->floatParsInit());
    *args = res->floatParsInit(); //resets arg values to prefit values
    delete args;
    delete res;
    m_signalStrength->setConstant(cstate);m_signalStrength->setVal(val);
    /*
    if(m_dataHist[1]->GetNbinsX()==1 && m_dataHist[3]->GetNbinsX()==1) {
      Error("Fit","Cannot float signal when signal region is blinded and only have 1 bin in the control regions");
      Error("Fit","Please fix the signal strength (with GetParameter(\"mu\")->setConstant(true))");
      return 0;
    }
    */
  }
  
  
 
  
  if(!m_model) BuildModel();

  double parEstimate = m_dataHist[2]->Integral()/m_dataHist[3]->Integral();
  if(parEstimate==0) parEstimate=1;
  if(m_modelPars[1]->isConstant()) {
    //doing flat model     
    //A = mB
    //C = mD
    
    m_modelPars[0]->setRange(0,parEstimate*100);
    m_modelPars[0]->setVal(parEstimate);
    m_modelPars[1]->setVal(0); //no gradient
    
    
  } else if(!m_modelPars[1]->isConstant()) {
    if( m_dataHist[1]->GetNbinsX()==1 && m_dataHist[3]->GetNbinsX()==1 ) {
      Error("Fit","Cannot use linear model when you only have 1 bin per region");
      return 0;
    }
  
    //doing linear model
    //A = (m1 + m2x)B
    //C = (m1 + m2x)D

    m_modelPars[0]->setRange(-5,parEstimate*100);
    m_modelPars[0]->setVal(parEstimate);
    m_modelPars[1]->setVal(0); //start at no gradient
  }
   
   
  //ensure the range of mu parameter is big enough to cover all data with signal 
  
  if(hasSignal && !m_signalStrength->isConstant()) {
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
    if(totSignal!=0) {
      m_signalStrength->setRange(0,(totData/totSignal)+1);
      Info("Fit","Setting signal strength to %g to test for maximal signal", (totData/totSignal));
      m_signalStrength->setVal((totData/totSignal));
    }
  }
  
  //blind the signal region if necessary
  RooDataSet* theData = m_data;
  if(m_dataHist[0]==0 || m_useAsimov || m_isVR) {
    theData = static_cast<RooDataSet*>(m_data->reduce("region!=0")); //removes the region A data
    //if we have asimovdata then don't blind, just use that ..
    if(!m_asimovData) {
      m_xVar->setRange("myRange",m_dataHist[1]->GetXaxis()->GetXmin(),m_dataHist[1]->GetXaxis()->GetXmax());
      m_stacks[0]->setBlindRange("myRange"); //the blind range will make the value of this component appear to be 0 in this range
    } else {
      theData->append( *m_asimovData ); 
    }
    
  }
  
  
  ///now run the fit
  m_lastFitResult = 0;
  
  RooMsgService::instance().setGlobalKillBelow(m_printLevel);

  
  
  //check which floating parameters have 'fixedValue' string attribute set. 
  //these parameters will have a fit run with 
  RooArgSet* params = m_model->getParameters( *theData );
  RooFIter itr(params->fwdIterator());
  RooAbsArg* p = 0;
  
  std::map<TString,TRooFitResult*> afr;
  while( (p = itr.next()) ) {
    if(!p->getStringAttribute("fixedValue")) continue;
    RooRealVar* v = dynamic_cast<RooRealVar*>(p);
    if(!v) continue;
    if(v->isConstant()) continue;
    Info("Fit","Running %s=%s Fit...",v->GetName(),v->getStringAttribute("fixedValue"));
    
    double oldVal = v->getVal();
    v->setVal( TString(v->getStringAttribute("fixedValue")).Atof() );
    v->setConstant(true);
    
    
    auto x = m_model->fitTo(*theData,RooFit::Save(),RooFit::Minimizer("Minuit2"),RooFit::Minos());
    TRooFitResult* res = new TRooFitResult(x);
    
    //special retry of fit with restricted mu range, since we played with the mu ranges above
    if(res->status()!=0 && !m_signalStrength->isConstant() && hasSignal) { 
      double oldMax = m_signalStrength->getMax();
      m_signalStrength->setRange(0,m_signalStrength->getVal()*2.);
      Warning("Fit","Floating signal fit failed, retrying with restricted range [%g,%g]...",0.,m_signalStrength->getMax());
      auto args = m_model->getObservables(m_lastFitResult->floatParsInit());
      *args = m_lastFitResult->floatParsInit(); //resets arg values to mu=0 postfit values
      m_signalStrength->setVal(m_signalStrength->getMax()/2.); //but put the signal back at its post fit value
      delete args;
      delete res;delete x;
      x = m_model->fitTo(*theData,RooFit::Save(),RooFit::Minimizer("Minuit2"));
      res = new TRooFitResult(x);
      m_signalStrength->setMax(oldMax);
    }
    
    //restore all parameters to pre-fit values ..
    auto args = m_model->getObservables(res->floatParsInit());
    *args = res->floatParsInit(); //resets arg values to mu=0 postfit values
    delete args;
    
    v->setVal(oldVal);
    v->setConstant(false);
    
    afr[Form("%s=%s",v->GetName(),v->getStringAttribute("fixedValue"))] = res;
    
  }
  
  
  
  if(m_printLevel>=RooFit::ERROR) {
    Info("Fit","Running Fit ... (use SetPrintLevel(RooFit::INFO) to show roofit output)...");
  }
  
  auto x = m_model->fitTo(*theData,RooFit::Save(),RooFit::Minimizer("Minuit2"),RooFit::Minos(),RooFit::PrintLevel(3));
  m_lastFitResult = new TRooFitResult(x);
  
  //special retry of fit with restricted mu range, since we played with the mu ranges above
  if(m_lastFitResult->status()%10!=0 && !m_signalStrength->isConstant() && hasSignal) { 
    m_signalStrength->setRange(0,m_signalStrength->getVal()*2.);
    Warning("Fit","Floating signal fit failed, retrying with restricted range [%g,%g]...",0.,m_signalStrength->getMax());
    auto args = m_model->getObservables(m_lastFitResult->floatParsInit());
    *args = m_lastFitResult->floatParsInit(); //resets arg values to mu=0 postfit values
    m_signalStrength->setVal(m_signalStrength->getMax()/2.); //but put the signal back at its post fit value
    delete args;
    delete x; delete m_lastFitResult;
    m_lastFitResult = new TRooFitResult(m_model->fitTo(*theData,RooFit::Save(),RooFit::Minimizer("Minuit2")));
  }
  
  
  if(m_useAsimov || m_dataHist[0]==0 || m_isVR) {
    m_stacks[0]->setBlindRange(""); //remove the blinding range for the post-fit plotting
    delete theData;
  }
  
  //add zero fit result as an associated fit result
  //m_lastFitResult->AddAssociatedFitResult("mu=0",zeroFitResult);
  
  for(auto& a : afr) {
    m_lastFitResult->AddAssociatedFitResult(a.first,a.second);
  }
  
  //check the postFit errors and ensure the ranges on parameters are adjusted to cover them ... otherwise we can underestimate uncerts 
  for(uint i=0;i<m_modelPars.size();i++) {
  // if( m_modelPars[i]->getVal()+m_modelPars[i]->getErrorLo() < m_modelPars[i]->getMin() ) m_modelPars[i]->setMin( m_modelPars[i]->getVal()+m_modelPars[i]->getErrorLo() - 1. );
  // if( m_modelPars[i]->getVal()+m_modelPars[i]->getErrorHi() > m_modelPars[i]->getMax() ) m_modelPars[i]->setMax( m_modelPars[i]->getVal()+m_modelPars[i]->getErrorHi() + 1. );
  }
  
  if(drawPostFit) {
    ///draw the post-fit results
    TCanvas* cc2 = new TCanvas(Form("%s_postfit",GetName()),Form("%s Post-fit",GetTitle()),800,750);
    Draw(m_lastFitResult);
    cc2->Modified(1);
    cc2->Update();
  }
  
  return m_lastFitResult;
}


void TRooABCD::Draw(TRooFitResult* fr, const char* canvasName) {
  if(!m_model) BuildModel(); //need all model components available to draw
  
  TCanvas* cc = 0;
  if(strlen(canvasName) || !gPad) {
     if(strlen(canvasName)==0) cc = new TCanvas(GetName(),GetTitle(),800,750);
     else cc = new TCanvas(canvasName,canvasName,800,750);
  }
  
  if(fr==0) fr = m_lastFitResult;
  
  gPad->Clear();
  
  auto currPad = gPad;
  
  gPad->Divide(2,3);
  
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    double maxData = m_dataHist[i]->GetBinContent( m_dataHist[i]->GetMaximumBin() );
    double maxPred(0);
    for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
      double v = m_stacks[i]->GetBinContent(j,fr);
      if(v>maxPred) maxPred=v;
    }
    if(maxData > maxPred) m_stacks[i]->SetMaximum( maxData*1.1 ); //ensures data is in axis range
    else m_stacks[i]->SetMaximum(-1111); //removes max
  }
  
  TText tt;
  
  currPad->cd(1);
  if(fr) {m_stacks[2]->Draw("e3005",fr);}
  else {m_stacks[2]->Draw("e3005");}
  m_dataHist[2]->Draw("same");tt.DrawTextNDC(gPad->GetLeftMargin(),1.-gPad->GetTopMargin()+.02,Form("Region %s",m_cat->lookupType(2)->GetName()));
  currPad->cd(2);
  if(fr) {m_stacks[0]->Draw("e3005",fr);}
  else {m_stacks[0]->Draw("e3005");}
  if(m_dataHist[0]) m_dataHist[0]->Draw("same");
  if(m_asimovDataHist[0]) m_asimovDataHist[0]->Draw("same");
  tt.DrawTextNDC(gPad->GetLeftMargin(),1.-gPad->GetTopMargin()+0.02,Form("Region %s (%s Region)",m_cat->lookupType(0)->GetName(),(m_isVR)?"Validation":"Signal"));
  currPad->cd(3);
  if(fr) {m_stacks[3]->Draw("e3005",fr); }
  else {m_stacks[3]->Draw("e3005");}
  m_dataHist[3]->Draw("same");tt.DrawTextNDC(gPad->GetLeftMargin(),1.-gPad->GetTopMargin()+0.02,Form("Region %s",m_cat->lookupType(3)->GetName()));
  currPad->cd(4);
  if(fr) {m_stacks[1]->Draw("e3005",fr); }
  else {m_stacks[1]->Draw("e3005");}
  m_dataHist[1]->Draw("same");tt.DrawTextNDC(gPad->GetLeftMargin(),1.-gPad->GetTopMargin()+0.02,Form("Region %s",m_cat->lookupType(1)->GetName()));
  
  //we want to draw the ratio of the data, after signal and other subtraction
  
  TH1* data_ratio = (TH1*)m_dataHist[2]->Clone("data_ratio"); data_ratio->SetTitle("Transfer factor"); data_ratio->SetDirectory(0);
  if(m_signal[2]) data_ratio->Add( m_signal[2]->GetHistogram(), -1. );
  if(m_other[2]) data_ratio->Add( m_other[2]->GetHistogram(), -1. );
  TH1* tmpHist = (TH1*)m_dataHist[3]->Clone("data_denom");
  if(m_signal[3]) tmpHist->Add( m_signal[3]->GetHistogram(), -1. );
  if(m_other[3]) tmpHist->Add( m_other[3]->GetHistogram(), -1. );
  data_ratio->Divide(tmpHist);
  //only show it for bins where there was data in both parts
  for(int i=1;i<=data_ratio->GetNbinsX();i++) {
    if(m_dataHist[2]->GetBinContent(i)==0 || m_dataHist[3]->GetBinContent(i)==0) { data_ratio->SetBinContent(i,0);data_ratio->SetBinError(i,0); }
  }
  
  delete tmpHist;
  
  if(m_dataHist[0]) {
    //unblinded, so will draw transfer factor with full range ..
    std::vector<double> binBoundaries;
    for(int i=1;i<=m_dataHist[3]->GetNbinsX();i++) binBoundaries.push_back( m_dataHist[3]->GetBinLowEdge(i) );
    for(int i=1;i<=m_dataHist[1]->GetNbinsX()+1;i++) binBoundaries.push_back( m_dataHist[1]->GetBinLowEdge(i) );
    
    TH1* data_ratio2 = (TH1*)m_dataHist[0]->Clone("data_ratio2");
    if(m_signal[0]) data_ratio2->Add( m_signal[0]->GetHistogram(), -1. );
    if(m_other[0]) data_ratio2->Add( m_other[0]->GetHistogram(), -1. );
    TH1* tmpHist = (TH1*)m_dataHist[1]->Clone("data_denom");
    if(m_signal[1]) tmpHist->Add( m_signal[1]->GetHistogram(), -1. );
    if(m_other[1]) tmpHist->Add( m_other[1]->GetHistogram(), -1. );
    data_ratio2->Divide(tmpHist);
    
    //only show it for bins where there was data in both parts
    for(int i=1;i<=data_ratio->GetNbinsX();i++) {
      if(m_dataHist[0]->GetBinContent(i)==0 || m_dataHist[1]->GetBinContent(i)==0) { data_ratio2->SetBinContent(i,0);data_ratio2->SetBinError(i,0); }
    }
    
    delete tmpHist;
    
    TH1* data_ratio1 = data_ratio;
    
    data_ratio = new TH1D("data_ratio","Transfer factor",binBoundaries.size()-1,&binBoundaries[0]); data_ratio->SetDirectory(0);
    *static_cast<TAttLine*>(data_ratio) = *static_cast<TAttLine*>(data_ratio1);
    *static_cast<TAttMarker*>(data_ratio) = *static_cast<TAttMarker*>(data_ratio1);
    data_ratio->Sumw2();
    
    for(int i=1;i<=data_ratio1->GetNbinsX();i++) {
      data_ratio->SetBinContent( data_ratio->FindFixBin( data_ratio1->GetBinCenter(i) ), data_ratio1->GetBinContent(i) );
      data_ratio->SetBinError( data_ratio->FindFixBin( data_ratio1->GetBinCenter(i) ), data_ratio1->GetBinError(i) );
    }
    for(int i=1;i<=data_ratio2->GetNbinsX();i++) {
      data_ratio->SetBinContent( data_ratio->FindFixBin( data_ratio2->GetBinCenter(i) ), data_ratio2->GetBinContent(i) );
      data_ratio->SetBinError( data_ratio->FindFixBin( data_ratio2->GetBinCenter(i) ), data_ratio2->GetBinError(i) );
    }
    delete data_ratio1; delete data_ratio2;
  }
  
  data_ratio->SetXTitle(m_xVar->GetTitle());
  data_ratio->SetStats(0);
  currPad->cd(5);data_ratio->Draw();
  if(fr) m_transferFactor->Draw("val L same e3005",fr);
  else  m_transferFactor->Draw("val L same e3005");
  
  TLatex t;
  
  if(fr && fr->floatParsFinal().find("m0")) {
    RooRealVar* m = static_cast<RooRealVar*>(fr->floatParsFinal().find("m0"));
    t.DrawLatexNDC(0.3,1.-gPad->GetTopMargin()+0.01,Form("%s = %g #pm %g",m->GetTitle(),m->getVal(),m->getError()));
  }
  
  if(fr && fr->floatParsFinal().find("m1")) {
    RooRealVar* m = static_cast<RooRealVar*>(fr->floatParsFinal().find("m1"));
    t.DrawLatexNDC(0.3,1.-gPad->GetTopMargin()+0.05,Form("%s = %g #pm %g",m->GetTitle(),m->getVal(),m->getError()));
  }

  currPad->cd(6);
  
  
  
  
  
  double prediction = GetBkgIntegral(0,fr);
  double predictionError = GetBkgIntegralError(0,fr);
  
  
  //check for fixedValue fit results ...
  //these are associated fit results with a single floating parameter missing ...
  std::map<RooRealVar*,TRooFitResult*> fixedFitResults;
  
  if(fr) {
    for(auto& afr : fr->GetAssociatedFitResults()) {
      if(!afr.second) continue; //just in case!
      RooArgList l( fr->floatParsFinal() );
      RooArgList l2( afr.second->floatParsFinal() );
      if(l.getSize()!=l2.getSize()+1) continue; //looking for fits with a single parameter lost
      l.remove(l2,true,true/*remove by name*/);
      if(l.getSize()!=1) continue;
      fixedFitResults[static_cast<RooRealVar*>(l.at(0))] = afr.second;
    }
  }
  
 
  double latexY = 0.8;
  
  double extraErr(0);
  
  t.DrawLatex(0.05,latexY,Form("%sR Bkg Predicted = %g #pm %g (syst.)",(m_isVR)?"V":"S",prediction,predictionError));latexY-=0.05;
  for(auto& ffr : fixedFitResults) {
    double prediction_ffr = GetBkgIntegral(0,ffr.second);
    t.DrawLatex(0.05,latexY,Form("                   #pm %g (%s syst.)",fabs(prediction_ffr-prediction),ffr.first->GetTitle()));latexY-=0.05;
    extraErr += pow( prediction-prediction_ffr, 2);
  }
  
  
  if(m_signal[0]!=0) {
    //fit had signal, so print signal prediction in SR
    t.DrawLatex(0.05,latexY,Form("%sR Signal Predicted = %g #pm %g",(m_isVR)?"V":"S",GetSignalIntegral(0),GetSignalIntegralError(0)));latexY-=0.05;
  }
  
  
  if(m_dataHist[0]) {
    predictionError = sqrt( pow(predictionError,2) + extraErr );
    double signif = RooStats::NumberCountingUtils::BinomialObsZ(m_dataHist[0]->Integral(),prediction,predictionError/prediction);
    
    if(signif < 1) t.SetTextColor(kGreen);
    else if(signif < 2.5) t.SetTextColor(kOrange);
    else t.SetTextColor(kRed);
    
    t.DrawLatex(0.05,latexY,Form("%sR Observed = %g (%.1f#sigma)",(m_isVR)?"V":"S",m_dataHist[0]->Integral(),signif));latexY-=0.05;
  }
  else {
    t.DrawLatex(0.05,latexY,Form("%sR Observed = BLINDED",(m_isVR)?"V":"S"));latexY-=0.05;
  }
  
  
  //exit here if we don't have a fit result
  if(!fr) {
    if(cc) {
      cc->Modified(1);
      cc->Update();
    }
    return;
  }

  
  const char* statusMeanings[6] = {"OK", 
                                  "Covariance was made pos definite", 
                                  "Hesse is invalid", 
                                  "Edm is above max", 
                                  "Reached call limit", 
                                  "Any other failure"};
  int statusCols[6] = {kGreen,kOrange,kRed,kRed,kRed,kRed};
  
  latexY-=0.05;
  
  
  
  t.SetTextColor(statusCols[fr->status()%10]);
  t.DrawLatex(0.05,latexY,Form("Fit Status = %d (%s)",fr->status(),statusMeanings[fr->status()%10]));
  latexY-=0.05;
  //compute chi^2 for the available data .
  double pVal = Get2LLRPValue(fr); //GetChi2PValue(fr);
  if(pVal<1) {
    if(RooStats::PValueToSignificance(pVal) < 1) t.SetTextColor(kGreen);
    else if(RooStats::PValueToSignificance(pVal) < 2.5) t.SetTextColor(kOrange);
    else t.SetTextColor(kRed);
    //t.DrawLatex(0.05,latexY,Form("Fit #chi^{2} p-value = %g (%.1f#sigma)",pVal,RooStats::PValueToSignificance(pVal)));
    t.DrawLatex(0.05,latexY,Form("Fit 2LLR p-value = %g (%.1f#sigma)",pVal,RooStats::PValueToSignificance(pVal)));
    latexY-=0.05;
  }
  
  for(auto& ffr : fixedFitResults) {
    int stat = ffr.second->status()%10;
    if(!(stat>=0 && stat<=5)) stat=5;
    t.SetTextColor(statusCols[stat]);
    t.DrawLatex(0.05,latexY,Form("%s=%s Fit Status = %d (%s)",ffr.first->GetTitle(),ffr.first->getStringAttribute("fixedValue"),ffr.second->status(),statusMeanings[stat]));
    latexY-=0.05;
    
    pVal = Get2LLRPValue(ffr.second); //GetChi2PValue(ffr.second);
    if(pVal<1) {
      if(RooStats::PValueToSignificance(pVal) < 1) t.SetTextColor(kGreen);
      else if(RooStats::PValueToSignificance(pVal) < 2.5) t.SetTextColor(kOrange);
      else t.SetTextColor(kRed);
      //t.DrawLatex(0.05,latexY,Form("%s=%s Fit #chi^{2} p-value = %g (%.1f#sigma)",ffr.first->GetTitle(),ffr.first->getStringAttribute("fixedValue"),pVal,RooStats::PValueToSignificance(pVal)));
      t.DrawLatex(0.05,latexY,Form("%s=%s Fit 2LLR p-value = %g (%.1f#sigma)",ffr.first->GetTitle(),ffr.first->getStringAttribute("fixedValue"),pVal,RooStats::PValueToSignificance(pVal)));
      latexY-=0.05;
    }
    
  }
  
   
    
    
  
  
  if(cc) {
    cc->Modified(1);
    cc->Update();
  }



}

double TRooABCD::GetChi2PValue(TRooFitResult* fr) {
  if(!fr) fr = m_lastFitResult;
  double chi2 = 0;
  int ndof=0;
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    TH1* model = m_stacks[i]->GetHistogram(fr,true);
    for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
      if(m_dataHist[i]->GetBinContent(j) != model->GetBinContent(j)) {
        double val = pow( (m_dataHist[i]->GetBinContent(j) - model->GetBinContent(j)) , 2 )/( pow(model->GetBinError(j),2)+model->GetBinContent(j) ); //using error = model poisson error + model error (in quadrature)
        if(m_printLevel<=RooFit::INFO) Info("Fit","chi2 for bin %d in region %s = %g",j,m_cat->lookupType(i)->GetName(),val);
        chi2 += val;
      }
      ndof++;
    }
    //pVal *= m_dataHist[i]->Chi2Test( model, "UW" );
  }
  
  //correct ndof to account for free parameters
  if(fr && fr->floatParsFinal().find("mu")) ndof--;
  ndof -= m_dataHist[1]->GetNbinsX();
  ndof -= m_dataHist[3]->GetNbinsX();
  
  if(fr && fr->floatParsFinal().find("m0")) ndof--;
  if(fr && fr->floatParsFinal().find("m1")) ndof--;
  
  
  if(ndof<=0) return 1.;
  
  return TMath::Prob(chi2,ndof);
}

double TRooABCD::Get2LLRPValue(TRooFitResult* fr, int nToys) {
  double obsLLR = Get2LLR(fr,GetDataSet());
  
  if(nToys==0) {
    //use chi2 approx with ndof = number of bins
    int nBins=0;
    for(int i=0;i<4;i++) {
      if(m_dataHist[i]) {
        nBins += m_dataHist[i]->GetNbinsX();
      }
    }
    return TMath::Prob(obsLLR,nBins);
  }
  
  double low=0;
  for(int i=0;i<nToys;i++) {
    auto toy = createToyDataSet();
    if( Get2LLR(fr,toy)>obsLLR) low++;
    delete toy;
  }
  
  if(low==0) low=1; //to avoid returning 0 pValue
  
  return low/nToys;
}

double TRooABCD::Get2LLR(TRooFitResult* fr, RooAbsData* data) {
  //log likelihood ratio 'test statistic' as studied by Baker-Cousins
  //The denominator of the ratio is the likelihood for the 'perfect' binned model
  
  //NOTE: 2*LLR is approximately chi2 distribution, ndof ~= nBins (comes out a little higher) ??
  
  if(!fr) fr = m_lastFitResult;
  
  //move parameter values onto fr values
  //and constant parameters too
  RooArgList floatAndConst; floatAndConst.add(fr->floatParsFinal()); floatAndConst.add(fr->constPars());
  auto params = m_model->getObservables(floatAndConst);
  auto snap = static_cast<RooArgSet*>(params->snapshot());
  *params = floatAndConst;

  double nll = 0;
  double nll2 = 0; //denominator
  
  //first activate blinding for regions without data
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) {
      m_xVar->setRange("myRange",m_xVar->getMin(),m_xVar->getMax());
      m_stacks[i]->setBlindRange("myRange");
    }
  }
  
  //now evaluate the NLL of the data
  RooAbsPdf* model = GetModel();
  if(data==0) data = GetDataSet();
  
  RooArgSet* obs = model->getObservables(data);
  for(int i=0;i<data->numEntries();i++) {
    *obs = *data->get(i); //move observables values of ith entry
    nll -= model->getLogVal( obs )*data->weight(); //get probability of ith entry (scale by weight)
  }
  //if model is extended, add extended term ..
  nll += model->extendedTerm(data->sumEntries(),obs);
  
  //ll2 is the likelihood from the 'saturated' model where each bin prediction equals the data
  //use data histograms to estimate this ...
  //formula is:
  // nll2 = - sum_i [ n_i * log(n_i/(w_i*N)) ] + N - NlogN = -sum_i [ n_i*log(n_i/w_i) ] + N
  //last equality comes from sum_i [n_i * log(N)] = Nlog N 
  
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
      double nObs = data->sumEntries(Form("region==%d&&x>=%f&&x<%f",i,m_dataHist[i]->GetBinLowEdge(j),m_dataHist[i]->GetBinLowEdge(j+1)));
      if(nObs<=0) continue;
      nll2 -= nObs*log(nObs/m_dataHist[i]->GetBinWidth(j));
      nll2 += nObs; //adds the 'N' term
    }
  }
  
  
  //remove blind ranges ..
  for(int i=0;i<4;i++) {
    m_stacks[i]->setBlindRange("");
  }
  
  
  //restore parameter values in model
  *params = *snap;
  delete params;
  delete snap;
  
  return 2.*(nll - nll2);
  
}

#include "RooRandom.h"

RooDataSet* TRooABCD::createToyDataSet(bool expected) {
  //create a toy dataset from the current state of the model, using the data hist binning ...
  //Question: should we be flucuating the underlying parameters, hybrid style?
  //frequentist approach would say not to, since we just use the best fit values as they are
  
  RooAbsPdf* model = GetModel();
  RooAbsData* data = GetDataSet();
  
  RooArgSet* obs = model->getObservables(data);
  
  double nExpected = model->expectedEvents(obs);
  RooRealVar weight("weight","weight",0,1e9) ;
  RooArgSet tmp(*obs) ;
  tmp.add(weight) ;
  
  RooDataSet* binnedData = new RooDataSet("toy","toy",tmp,RooFit::WeightVar("weight"));
  
  for(int i=0;i<4;i++) {
    if(!m_dataHist[i]) continue;
    m_cat->setIndex(i);
    for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
      m_xVar->setVal( m_dataHist[i]->GetBinCenter(j) );
      double binExpectation = model->getVal( data->get() )*nExpected*m_dataHist[i]->GetBinWidth(j); //exploits flatness of pdf across bin width
      double w = (expected) ? binExpectation : RooRandom::randomGenerator()->Poisson( binExpectation );
      binnedData->add( *obs, w);
    }
  }
  
  delete obs;
  
  return binnedData;

}
