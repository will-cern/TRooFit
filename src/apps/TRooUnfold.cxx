

#include "TRooFit/apps/TRooUnfold.h"

#include "TRooFit/Utils.h"

#include "TLatex.h"
#include "TCanvas.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"
#include "RooBinning.h"

#include "RooDataHist.h"

#include "RooMinimizer.h"

TRooUnfold::TRooUnfold(const char* name, const char* title,int nBinsReco,double recoLow,double recoHigh,int nBinsTruth,double truthLow,double truthHigh) : TNamed(name,title) {
  m_recoVar = new RooRealVar("reco","reco",recoLow,recoHigh);
  m_truthVar = new RooRealVar("truth","truth",truthLow,truthHigh);
  m_recoVar->setBins(nBinsReco);
  m_truthVar->setBins(nBinsTruth);
  
}

TRooUnfold::TRooUnfold(const char* name, const char* title,TH2* nominalMigrationMatrix) : TNamed(name,title) {
  m_recoVar = new RooRealVar("reco","reco",nominalMigrationMatrix->GetXaxis()->GetXmin(),nominalMigrationMatrix->GetXaxis()->GetXmax());
  m_truthVar = new RooRealVar("truth","truth",nominalMigrationMatrix->GetYaxis()->GetXmin(),nominalMigrationMatrix->GetYaxis()->GetXmax());
  m_recoVar->setBins(nominalMigrationMatrix->GetNbinsX());
  m_truthVar->setBins(nominalMigrationMatrix->GetNbinsY());
  
  if(nominalMigrationMatrix->GetXaxis()->GetXbins()->GetArray()) {
    m_recoVar->setBinning(RooBinning(nominalMigrationMatrix->GetNbinsX(),nominalMigrationMatrix->GetXaxis()->GetXbins()->GetArray()));
  }
  if(nominalMigrationMatrix->GetYaxis()->GetXbins()->GetArray()) {
    m_truthVar->setBinning(RooBinning(nominalMigrationMatrix->GetNbinsY(),nominalMigrationMatrix->GetYaxis()->GetXbins()->GetArray()));
  }
  
  AddMigrationMatrix(nominalMigrationMatrix);
  
}

Bool_t TRooUnfold::AddMigrationMatrix(TH2* migrationMatrix, const char* variation) {
  TString sVar(variation);
  if(m_migrationMatrix[sVar]==0) {
    m_migrationMatrix[sVar] = (TH2*)migrationMatrix->Clone(Form("migration_%s",sVar.Data()));
    m_migrationMatrix[sVar]->SetDirectory(0);
  } else {
    m_migrationMatrix[sVar]->Add(migrationMatrix);
  }
  return kTRUE;
}

Bool_t TRooUnfold::AddTruthSignal(TH1* truth, const char* variation) {
  if(m_truthSignal[variation]==0) {
    m_truthSignal[variation] = (TH1*)truth->Clone(Form("truthSignal_%s",variation));
    m_truthSignal[variation]->SetDirectory(0);
  } else {
    m_truthSignal[variation]->Add(truth);
  }
  return kTRUE;
}

Bool_t TRooUnfold::AddRecoSignal(TH1* reco, const char* variation) {
  if(m_recoSignal[variation]==0) {
    m_recoSignal[variation] = (TH1*)reco->Clone(Form("recoSignal_%s",variation));
    m_recoSignal[variation]->SetDirectory(0);
  } else {
    m_recoSignal[variation]->Add(reco);
  }
  return kTRUE;
}

Bool_t TRooUnfold::AddBackground(TH1* reco, const char* bkgName, const char* variation) {
  if(m_bkg.find(bkgName)==m_bkg.end()) {
    m_bkgOrder.push_back(bkgName);
  } 
  if(m_bkg[bkgName][variation]==0) {
    m_bkg[bkgName][variation] = (TH1*)reco->Clone(Form("%s_%s",bkgName,variation));
    m_bkg[bkgName][variation]->SetDirectory(0);
  } else {
    m_bkg[bkgName][variation]->Add(reco);
  }
  return kTRUE;
}

Bool_t TRooUnfold::AddScaleFactor(TH1* reco, const char* sfName, const char* variation) {
  if(m_recoSF[sfName][variation]==0) {
    m_recoSF[sfName][variation] = (TH1*)reco->Clone(Form("%s_%s",sfName,variation));
    m_recoSF[sfName][variation]->SetDirectory(0);
  } else {
    m_recoSF[sfName][variation]->Add(reco);
  }
  return kTRUE;
}

Bool_t TRooUnfold::AddNormFactor(const char* sfName, const char* sfTitle) {
  if(m_normFactor[sfName]) {
    Error("AddNormFactor","%s already exists",sfName);
    return kFALSE;
  }
  m_normFactor[sfName] = new RooRealVar(sfName,(sfTitle)?sfTitle:sfName,1,0,10);
  return kTRUE;
}

Bool_t TRooUnfold::ApplyFactorToComponent(const char* sfName, const char* compName) {
  m_sfApplied[sfName].insert(compName);
  return kTRUE;
}

Bool_t TRooUnfold::AddData(TH1* data) {
  if(!m_data) {
    m_data = (TH1*)data->Clone("data");
    m_data->SetDirectory(0);
    m_data->SetMarkerStyle(20); 
    m_data->SetMarkerSize(0.7); 
    m_data->SetLineColor(kBlack);
  } else {
    m_data->Add(data);
  }
  return kTRUE;
}

RooRealVar* TRooUnfold::GetSigNormVar(int bin) {
  int nDigits = int(std::log10(m_truthVar->numBins()))+1;
  TString formString(Form("mu_%%0%dd",nDigits));
  RooRealVar* out = m_normFactor[Form(formString,bin)];
  if(!out) {
    AddNormFactor(Form(formString,bin),Form("#mu_{%d}",bin));
    m_poi.add( *m_normFactor[Form(formString,bin)] );
    out = m_normFactor[Form(formString,bin)];
  }
  return out;
  
}

Bool_t TRooUnfold::AddRegularization(TH1* hist) {
  if(!m_regularizationHist) {
    m_regularizationHist = (TH1*)hist->Clone("regularization");
    m_regularizationHist->SetDirectory(0);
    m_regularizationHist->SetLineColor(kRed);
    m_regularizationHist->SetFillColor(kRed);
    m_regularizationHist->SetFillStyle(3004);
  } else {
    m_regularizationHist->Add(hist);
  }
  return kTRUE;
}

Bool_t TRooUnfold::BuildModel() {

  if(m_builtModel) return kFALSE;

  //define unconstrained scale factors for each truth bin, we will later need to adjust the ranges of these parameters to ensure we can cover the data ...
  for(int i=1;i<m_truthVar->numBins()+1;i++) {
    /*if(!m_normFactor[Form("signalNorm_bin%d",i)]) {
      AddNormFactor(Form("signalNorm_bin%d",i),Form("SigNorm_%d",i));
      m_poi.add( *m_normFactor[Form("signalNorm_bin%d",i)] );
    }*/
    RooRealVar* f = GetSigNormVar(i);
    ApplyFactorToComponent(f->GetName(),Form("signal_bin%d",i));
  }

  //start by deriving fiducial correction factors, if required (i.e. convert recoSignal histograms into purity histograms and add as scale factors for each of the signals)
  
  if(m_recoSignal.size()) {
    for(auto myPair : m_recoSignal) {
      TH1* reco = myPair.second;
      if(!reco) continue;
      TH2* migrationMatrix = m_migrationMatrix[myPair.first];
      if(!migrationMatrix) {
        Error("BuildModel","No migration matrix available for variation %s despite reco distribution given",myPair.first.Data());
        return kFALSE;
      }
      //project migration matrix onto x-axis .. this will be the numerator of the purity 
      TH1* recodTruth = migrationMatrix->ProjectionX();
      recodTruth->SetDirectory(0);
      recodTruth->Divide( reco );
      recodTruth->SetTitle(Form("Fiducial Purity (%s)",myPair.first.Data()));
      m_fidPurity[myPair.first] = recodTruth;
      
      //check if any of the bins are > 1 ... which would indicate a problem ..
      if(recodTruth->GetBinContent( recodTruth->GetMaximumBin() ) > 1.000001) {
        Error("BuildModel","Fiducial Purity > 1 (%f) found in bin %d", recodTruth->GetBinContent( recodTruth->GetMaximumBin() ), recodTruth->GetMaximumBin());
        return kFALSE;
      }
      
      //need to register inverse purity as a scale factor ...
      TH1* fidCorrection = (TH1*)recodTruth->Clone("fidCorrection");
      fidCorrection->Reset();
      for(int i=1;i<=fidCorrection->GetNbinsX();i++) { fidCorrection->SetBinContent(i,1.); }
      fidCorrection->Divide(recodTruth); //is now the inverse purity, add as a scale factor on all the signal components 
      //only do so if the invPurity in any bin is sufficiently different from 1 
      bool nonUnity(false);
      for(int i=1;i<=fidCorrection->GetNbinsX();i++) { if(fabs(fidCorrection->GetBinContent(i)-1.)>1e-5) {nonUnity=true; break; } }
      if(nonUnity) {
        AddScaleFactor(fidCorrection,"fidCorrection",myPair.first);
      } 
      delete fidCorrection;
      
    }
    Info("BuildModel","Built purity histograms");
    //was invPurity ever defined?
    if(m_recoSF.find("fidCorrection")!=m_recoSF.end()) {
      for(int i=1;i<=m_truthVar->numBins()+1;i++) ApplyFactorToComponent("fidCorrection",Form("signal_bin%d",i));
    } else {
      Info("BuildModel","Skipping fiducial correction factor");
    }
  }


  //Start by converting the migration matrices into response matrices 
  //response matrices are first normalized, then scaled by the individual bin global efficiencies (if defined)
  for(auto myPair : m_migrationMatrix) {
    //check for corresponding truthSignal histogram 
    TH1* truth = m_truthSignal[myPair.first];
    if(truth==0) truth = m_truthSignal["Nominal"];
    TH2* migrationMatrix = myPair.second;
    if(!migrationMatrix) continue;
    for(int i=1;i<=migrationMatrix->GetNbinsY();i++) {
      double rowTotal(0); 
      for(int j=1;j<=migrationMatrix->GetNbinsX();j++) {rowTotal += migrationMatrix->GetBinContent(j,i); }
      if(truth) {
        double truthBinTotal = truth->GetBinContent(i);
        if(truthBinTotal < rowTotal) {
          Error("BuildModel","TruthSignal in truth bin %d is less than yield in the migration matrix (variation=%s)",i,myPair.first.Data());
          return kFALSE;
        }
        rowTotal = truthBinTotal;
      }
      for(int j=1;j<=migrationMatrix->GetNbinsX();j++) { 
        migrationMatrix->SetBinContent(j,i,migrationMatrix->GetBinContent(j,i)/rowTotal);
        migrationMatrix->SetBinError(j,i,migrationMatrix->GetBinError(j,i)/rowTotal);
      }
      
      //since we converted migration into response matrix, we want to upscale the unconstrained scale factors to match ... 
      //only do this if doing Nominal ..
      if(myPair.first=="Nominal") {
        GetSigNormVar(i)->setRange(0,rowTotal*5.);
        GetSigNormVar(i)->setVal(rowTotal);
      }
      
    }
    
    migrationMatrix->SetTitle(Form("Response Matrix (%s)",myPair.first.Data()));
    
    //save efficiency ... 
    TH1* recodTruth = migrationMatrix->ProjectionY(); //migration matrix has already been normalized
    recodTruth->SetTitle(Form("Efficiency (%s)",myPair.first.Data()));
    recodTruth->SetDirectory(0);
    m_efficiency[myPair.first] = recodTruth;
    
    
  }
  
  
  
  //now for each migration matrix we will create the signal components of each bin
  for(auto myPair : m_migrationMatrix) {
    TH2* migrationMatrix = myPair.second;
    if(!migrationMatrix) continue;
    
    for(int i=1;i<=migrationMatrix->GetNbinsY();i++) {
      TH1* tmpHist = migrationMatrix->ProjectionX("myTmpHist",i,i); //create a temporary histogram for the response of this truth bin .. 
      if(!AddBackground(tmpHist,Form("signal_bin%d",i),myPair.first)) {
        Error("BuildModel","Failed to build signal term for truth bin %d in variation %s",i,myPair.first.Data());
        return kFALSE;
      }
      delete tmpHist;
    }
  }
  
  
  //Time to build scale factor functions ...
  for(auto myPair : m_recoSF) {
    //load the "Nominal" for this SF ...
    TH1* nomHist = myPair.second["Nominal"];
    if(!nomHist) {
      Error("BuildModel","%s does not have a nominal histogram",myPair.first.Data());
      return kFALSE;
    }
    m_sfFunctions[myPair.first] = new TRooHF1D(myPair.first,myPair.first,*m_recoVar,nomHist);
    //go through variations, adding to function ..
    for(auto v : myPair.second) {
      if(v.first=="Nominal") continue;
      if(!v.second) continue;
      TString vName = (v.first.Index("__")==-1) ? v.first : v.first(0,v.first.Index("__")); //format is "VARIATION__1" or "VARIATION__-1" for example
      if(m_systNP[vName]==0) {
        m_systNP[vName] = new RooRealVar(vName,vName,0,-5,5);
        m_systNP[vName]->setStringAttribute("constraintType","normal");
        if(m_systGroup.find(vName)==m_systGroup.end()) m_systGroup[vName]=vName;
        m_systNP[vName]->setAttribute(m_systGroup[vName],kTRUE);
        m_systNP[vName]->setError(1);
      }
      double vVal = (vName==v.first) ? 1. : TString(v.first(v.first.Index("__")+2,v.first.Length())).Atof();
      m_sfFunctions[myPair.first]->AddVariation(*m_systNP[vName],vVal,v.second);
    }
  }
  
  
  //now build each component and scale by required scale or normalization factors ...
  for(auto myPair : m_bkg) {
    //load the "Nominal" for this bkg ...
    TH1* nomHist = myPair.second["Nominal"];
    if(!nomHist) {
      Error("BuildModel","%s does not have a nominal histogram",myPair.first.Data());
      return kFALSE;
    }
    m_bkgPdfs[myPair.first] = new TRooH1D(myPair.first,myPair.first,*m_recoVar,nomHist);
    //go through variations, adding to pdf ..
    for(auto v : myPair.second) {
      if(v.first=="Nominal") continue;
      if(!v.second) continue;
      TString vName = (v.first.Index("__")==-1) ? v.first : v.first(0,v.first.Index("__"));
      if(m_systNP[vName]==0) {
        m_systNP[vName] = new RooRealVar(vName,vName,0,-5,5);
        m_systNP[vName]->setStringAttribute("constraintType","normal");
        if(m_systGroup.find(vName)==m_systGroup.end()) m_systGroup[vName]=vName;
        m_systNP[vName]->setAttribute(m_systGroup[vName],kTRUE);
        m_systNP[vName]->setError(1);
      }
      double vVal = (vName==v.first) ? 1. : TString(v.first(v.first.Index("__")+2,v.first.Length())).Atof();
      m_bkgPdfs[myPair.first]->AddVariation(*m_systNP[vName],vVal,v.second);
    }
    //go through scale factors, adding as norm factors 
    for(auto sf : m_sfApplied) {
      if(sf.second.find(myPair.first)!=sf.second.end()) {
        if(m_normFactor.find(sf.first)!=m_normFactor.end()) {
          m_bkgPdfs[myPair.first]->addNormFactor(*m_normFactor[sf.first]);
        } else {
          m_bkgPdfs[myPair.first]->addNormFactor(*m_sfFunctions[sf.first]);
        }
      }
    }
  }
  
  //go through bkg in order they were added and add to our master stack ..
  m_stack = new TRooHStack("recoStack","reco model");
  int cols[15] = {kRed,kCyan,kGreen+1,kOrange,kAzure-1,kGray,kPink+1,kYellow-7,kBlue+2,kRed-2,kViolet,kAzure-9,kViolet-9,kSpring-9,kCyan+2};
  int i=0;
  for(auto bkgName : m_bkgOrder) {
    m_bkgPdfs[bkgName]->SetFillColor(cols[i%15]);
    i++;
    m_stack->Add( m_bkgPdfs[bkgName] );
  }
  if(m_floor>=0) m_stack->setFloor(true,1e-9);
  
  
 
  
  //build a 'result' histogram, which is just a histogram function containing the values of the signal norm terms
  m_result = new TRooH1D("sigNorm","Signal Normalizations",*m_truthVar);
  for(int i=1;i<=m_truthVar->numBins();i++) {
    m_result->Fill(m_result->GetXaxis()->GetBinCenter(i),*GetSigNormVar(i));
  }
  
  //save as prefit distribution
  m_prefitTruth = (TH1*)m_result->GetHistogram()->Clone("prefit"); m_prefitTruth->SetDirectory(0); m_prefitTruth->SetLineColor(kRed);
  
  //build the data from the data hist 
  m_dataSet = new RooDataHist("obsData","obsData",*m_recoVar,m_data);
  m_fullModel = TRooFit::BuildModel( *m_stack, *m_dataSet );
  
  
  //get all floating parameters and any that have constraintType attribute, decorate 'OTHER' attribute ..
  RooArgSet* allPars = m_fullModel->getParameters(RooArgSet());
  RooFIter itrr(allPars->fwdIterator());
  while( RooAbsArg* arg = itrr.next() ) {
    if(arg->getStringAttribute("constraintType")) {
      arg->setAttribute("OTHER",kTRUE);
      if(m_systGroup.find(arg->GetName())==m_systGroup.end()) {
        TString aName(arg->GetName());
        if(aName.BeginsWith("recoStack_stat_bin")) {
          m_systGroup[arg->GetName()]="MCSTAT";
          arg->setAttribute("MCSTAT",kTRUE);
        } else {
          m_systGroup[arg->GetName()]="OTHER";
        }
      }
    }
  }
  delete allPars;
  
  //if regularization given, construct a regularization constraint term ...
  if(m_regularizationHist) {
    // Penalty term is equivalent to a product of gaussian constraints
    // This is done to reduce oscillations that occur when significant bin migration
    // Penality is (truth - mc_truth)/sigma^2 ... where 1/sigma^2 is the regularization strength
    // when the strength is 0, then there is effectively 0 penalty for pulling the truth from the mc_truth
    //
    // regularization strength is taken from the bin errors
    //
    //Note that adding a regularization term like this will bias the fit!
    RooDataHist* truth_ref = new RooDataHist("truth_ref","truth_ref",*m_truthVar,m_regularizationHist);
    m_regStrength = new RooRealVar("regStrength","Regularization strength",1);
    TMatrixD k(m_regularizationHist->GetNbinsX(),m_regularizationHist->GetNbinsX());
    for(int i=0;i<m_regularizationHist->GetNbinsX();i++) k(i,i)=pow(m_regularizationHist->GetBinError(i+1),2);
    m_regularizationConstraint = new TRooGPConstraint("penalty","penalty",*m_result,*truth_ref,k,false,m_regStrength);
    m_fullModel = new RooProdPdf(Form("%s_with_reg",m_fullModel->GetName()),Form("%s with regularization",m_fullModel->GetTitle()),RooArgList(*m_fullModel,*m_regularizationConstraint));
  }
  
  
  
  
  m_npAndUnconstrained.removeAll();
  for(auto np : m_systNP) {
    m_npAndUnconstrained.add(*np.second);
  }
  for(auto normFactor : m_normFactor) {
    m_npAndUnconstrained.add(*normFactor.second);
  }
  
  m_nll = TRooFit::createNLL(m_fullModel,m_dataSet,0);
  

  m_builtModel = true;

  return kTRUE;

}

void TRooUnfold::WriteModelInputs(const char* file,const char* opt) {
  TFile f(file,opt);
  
  for(auto myPair : m_fidPurity) {
    if(!f.GetDirectory(myPair.first)) f.mkdir(myPair.first);
    f.GetDirectory(myPair.first)->cd();
    myPair.second->Write("fidPurity");
  }
  for(auto myPair : m_efficiency) {
    if(!f.GetDirectory(myPair.first)) f.mkdir(myPair.first);
    f.GetDirectory(myPair.first)->cd();
    myPair.second->Write("efficiency");
  }
  for(auto myPair : m_recoSF) {
    for(auto myVariation : myPair.second) {
      if(!f.GetDirectory(myVariation.first)) f.mkdir(myVariation.first);
      f.GetDirectory(myVariation.first)->cd();
      myVariation.second->Write(myPair.first);
    }
  }
  for(auto myPair : m_bkg) {
    //want to record list of normFactors and scaleFactors that multiply each term 
    TString sfApplied("");
    for(auto mySF : m_sfApplied) {
      if( mySF.second.find(myPair.first)!=mySF.second.end() ) {
        if(sfApplied!="") sfApplied+=",";
        sfApplied += mySF.first;
      }
    }
    for(auto myVariation : myPair.second) {
      myVariation.second->SetTitle(Form("%s [%s]",myPair.first.Data(),sfApplied.Data()));
      if(!f.GetDirectory(myVariation.first)) f.mkdir(myVariation.first);
      f.GetDirectory(myVariation.first)->cd();
      myVariation.second->Write(myPair.first);
    }
  }
  
  f.Close();
}

void TRooUnfold::WriteHistograms(const char* file,const char* opt) {
  BuildModel(); //ensures model is built
  
  TFile f(file,opt);
  
  m_data->Write("data");
  
  
  TH1* h = m_stack->GetHistogram(m_fitResult,true);
  h->SetFillStyle(3005);h->SetFillColor(kBlue);h->Write("reco");
  m_stack->GetStack(m_fitResult)->Write("recoStack");
  
  if(m_fitResult) {
    auto resultHists = GetPostfitTruthHistograms();
    for(auto h : resultHists) {
      h->Write();
      delete h;
    }
  }
  
  //m_result->GetHistogram(m_fitResult,true)->Write("unfolded");
  
  f.Close();
  
}

Bool_t TRooUnfold::ClearSignalMCUncert() {

  //clear migration matrices ..
  for(auto myPair : m_migrationMatrix) {
    TH2* migrationMatrix = myPair.second;
    if(!migrationMatrix) continue;
    
    for(int i=1;i<=migrationMatrix->GetNbinsY();i++) {
      for(int j=1;j<=migrationMatrix->GetNbinsX();j++) {
        migrationMatrix->SetBinError(j,i,0);
      }
    }
  }
  
  //clear the recoSignal hists
  for(auto myPair : m_recoSignal) {
      TH1* reco = myPair.second;
      if(!reco) continue;
      for(int j=1;j<=reco->GetNbinsX();j++) {
        reco->SetBinError(j,0);
      }
  }
  
  //clear the truth signal hists
  
  for(auto myPair : m_truthSignal) {
      TH1* truth = myPair.second;
      if(!truth) continue;
      for(int j=1;j<=truth->GetNbinsX();j++) {
        truth->SetBinError(j,0);
      }
  }
  
  return kTRUE;

}

//NLL without regularization
Double_t TRooUnfold::GetNLLMin() {
  double tmp = m_regStrength->getVal();
  m_regStrength->setVal(0);
  auto result = TRooFit::minimize(m_nll);
  
  

  m_nllMin = m_nll->getVal();
  
  m_regStrength->setVal(1e-9); //trick for getting hold of ConstraintMax ...
  m_constraintMax = -m_regularizationConstraint->getLogVal() / 1e-9;
  
  m_regStrength->setVal(tmp);
  
  //restore initial float par values ...
  RooArgSet* pars = m_nll->getObservables( result->floatParsInit() );
  *pars = result->floatParsInit();
  delete pars;
  delete result;
  
  return m_nllMin;
  
}

//NLL when fixing POI to regularization (prior) values
//we do this by running a fit where all the POI are fixed at prior values ...
Double_t TRooUnfold::GetNLLMax() {
  std::vector<double> tmp;
  
  for(int i=1;i<=m_truthVar->numBins();i++) {
    RooRealVar* poi = GetSigNormVar(i);
    tmp.push_back(poi->getVal());
    poi->setVal( m_regularizationHist->GetBinContent(i) );
    poi->setConstant(true);
  }
  
  
  RooFitResult* result = 0;
  RooArgSet* pars = m_nll->getParameters(RooArgSet());
  auto fpars = pars->selectByAttrib("Constant",0);
  if(fpars->getSize()!=0) {
    result = TRooFit::minimize(m_nll);
  }
  delete pars;delete fpars;
  
  m_nllMax = m_nll->getVal();
  
  //restore initial float par values ...
  if(result) {
    RooArgSet* pars = m_nll->getObservables( result->floatParsInit() );
    *pars = result->floatParsInit();
    delete pars;
    delete result;
  }
  
  //and refloat the poi
  for(int i=1;i<=m_truthVar->numBins();i++) {
    RooRealVar* poi = GetSigNormVar(i);
    poi->setVal( tmp[i-1] );
    poi->setConstant(false);
  }

  return m_nllMax;
}

//constraint min by definition is 0
Double_t TRooUnfold::GetConstraintNLLMax() {
  return m_constraintMax;
}

Double_t TRooUnfold::GetNLL() {
  return m_nll->getVal() + m_regularizationConstraint->getLogVal();
}

Double_t TRooUnfold::GetConstraintNLL() {
  return -m_regularizationConstraint->getLogVal()/m_regStrength->getVal();
}



Double_t TRooUnfold::ScanRegularizationStrength(TGraph** g) {
  
  GetNLLMax();
  GetNLLMin();
  
  if(m_nllMax <= m_nllMin) {
    Warning("ScanRegularizationStrength","Equivalent or better fit found from full regularization, returning 0");  //happens if fit was 'perfect' e.g. closure test
    return 0;
  }
  
  if(m_constraintMax <= 0) {
    Warning("ScanRegularizationStrength","Constraint term is 0 in unregularized fit, returning 0"); //happens if fit was 'perfect' e.g. closure test
    return 0;
  }
  
  //we scan values of regularizationStrength, at each point evaluating the NLL and the Constraint term 
  
  if(g) *g = new TGraph;
  
  //double lastNllNorm = 0;
  //double lastConstraintNorm = 1;
  double lastDiff = 1;
  double lastTau = 0;
  
  for(float i=-10;i<=10;i+=0.5) {
    m_regStrength->setVal( std::pow(10,i) );
    auto uFit = TRooFit::minimize(m_nll); 
    delete uFit;
    double nllNorm = (GetNLL()-m_nllMin)/(m_nllMax-m_nllMin);
    double constraintNorm = GetConstraintNLL()/m_constraintMax;
    double diff =  nllNorm-constraintNorm ;
    if(g) (*g)->SetPoint((*g)->GetN(),i,diff);
    if(diff < 0) {
     lastTau = std::pow(10,i);
      //lastNllNorm = nllNorm;
      //lastConstraintNorm = constraintNorm;
      lastDiff = nllNorm - constraintNorm;
     continue;
    }
    
    return lastTau - lastDiff*(std::pow(10,i)-lastTau)/(diff - lastDiff);
  }
  
  return 99999;
  
}


RooFitResult* TRooUnfold::Fit(TH1* data, bool doMinos) {
  
  if(!m_data) AddData(data);
  if(!m_builtModel) BuildModel();
  
 
  if(m_fitResult) delete m_fitResult;
  
  if(data) {
    //rebuild the dataset
    if(m_dataSet) delete m_dataSet;
    m_dataSet = new RooDataHist("obsData","obsData",*m_recoVar,m_data);
    if(m_nll) delete m_nll;
    m_nll = TRooFit::createNLL(m_fullModel,m_dataSet,0);
  }


  TRooFit::setRecommendedDefaultOptions();
  
  Info("Fit","Running unconditional global fit ... ");
  auto fitResult = TRooFit::minimize(m_nll); //run unconditional fit 
  
  //now run minos errors for poi
  if(doMinos && fitResult->status()%1000 == 0) {
    Info("Fit","Computing Minos errors for parameters of interest ... ");
    TRooFit::minos(m_nll,m_poi,fitResult);
  }
  
  
  
  //m_fitResult = new TRooFitResult( fitResult ); //extends functionality of RooFitResult
  //delete fitResult; //no longer needed because TRooFitResult should copy all it needs
  m_fitResult = fitResult;

  return m_fitResult;

}


TH1* TRooUnfold::GetExpectedRecoHistogram(RooFitResult* fr) {
  return (TH1*)m_stack->GetHistogram((fr==0)?m_fitResult:fr,true,0)->Clone("expected");
}

TH1* TRooUnfold::GetSignificanceHistogram() {
  //build the reco-level histogram containing the significance in each bin, evaluated
  //using the errors from standard (symmetric) error propagation of the covariance matrix
  //Note that the covariance matrix should already have been 'inflated' to cover minos errors
  TH1* out = new TH1D("significance","significance",m_recoVar->numBins(),m_recoVar->getBinning().array());
  out->SetDirectory(0);
  TH1* expectation = GetExpectedRecoHistogram();
  
  for(int i=1;i<=out->GetNbinsX();i++) {
    double expected = expectation->GetBinContent(i);
    double err = expectation->GetBinError(i);
    double data = m_data->GetBinContent(i);
    if(expected) out->SetBinContent(i,RooStats::NumberCountingUtils::BinomialObsZ(data,expected,err/expected));
  }
  
  delete expectation;
  return out;
}

TH1* TRooUnfold::GetResidualHistogram(RooFitResult* fr) {
  //the residual is the difference between the data and the prediction at reco level, with the error coming from the reco prediction (the data is taken to be 'exact') ...
  TH1* out = new TH1D("residual","residual",m_recoVar->numBins(),m_recoVar->getBinning().array());out->Sumw2();
  out->SetDirectory(0);
  TH1* expectation = GetExpectedRecoHistogram(fr);
  
  for(int i=1;i<=out->GetNbinsX();i++) {
    double expected = expectation->GetBinContent(i);
    double err = expectation->GetBinError(i);
    double data = m_data->GetBinContent(i);
    out->SetBinContent(i,data-expected);
    out->SetBinError(i,err);
  }
  
  delete expectation;
  return out;
  
  
}

TH1* TRooUnfold::GetPrefitTruthHistogram() {
  return m_prefitTruth;
}

std::vector<TH1*> TRooUnfold::GetPostfitTruthHistograms(const std::vector<TString>&& systGroups, bool doSTATCORR) {

  std::vector<TString> myGroups(systGroups);

  if(myGroups.size()==0) {
    /*RooArgSet* allPars = m_nll->getObservables(m_fitResult->floatParsFinal());
    RooFIter itr = allPars->fwdIterator();
    RooAbsArg* arg;
    bool hasMCSTAT=false;
    bool hasSYST=false;
    while( (arg = itr.next()) ) {
      TString aName(arg->GetName());
      if(aName.BeginsWith("recoStack_stat_bin")) {arg->setAttribute("MCSTAT"); hasMCSTAT=true;}
      else if(!aName.BeginsWith("mu_")) {arg->setAttribute("SYST"); hasSYST=true;}
    }
    if(hasMCSTAT) myGroups.push_back("MCSTAT");
    if(hasSYST) myGroups.push_back("SYST");
    */
    
    for(auto& s : m_systGroup) {
      if(std::find(myGroups.begin(),myGroups.end(),s.second)==myGroups.end()) myGroups.push_back(s.second);
    }
    
  }
  
  std::map<TString,RooArgList*> uncert_breakdown = TRooFit::breakdown(m_nll,m_poi,myGroups,m_fitResult,false,doSTATCORR);

  std::vector<TH1*> out;
  
  
  std::vector<TString> allGroups;
  allGroups.push_back("nominal");allGroups.push_back("TOTAL");allGroups.push_back("STAT");
  if(doSTATCORR) allGroups.push_back("STATCORR");
  
  for(auto syst : myGroups) allGroups.push_back(syst);
  
  
  for(auto group : allGroups) {
    out.push_back(new TH1D(Form("unfold_%s",group.Data()),Form("%s",group.Data()),m_truthVar->numBins(),m_truthVar->getBinning().array()));
    out.back()->Sumw2();out.back()->SetDirectory(0);
    if(group!="nominal"&&group!="STAT") {
      out.back()->SetFillColor(out.back()->GetLineColor());out.back()->SetFillStyle(3004);
    }
    
    
    for(int i=1;i<=out.back()->GetNbinsX();i++) {
      RooRealVar* v = static_cast<RooRealVar*>(uncert_breakdown[(group=="nominal")?"TOTAL":group]->find(GetSigNormVar(i)->GetName()));
      if(group=="nominal") {
        out.back()->SetBinContent(i,v->getVal());
      } else {
        out.back()->SetBinContent(i,v->getVal()+(v->getErrorHi()+v->getErrorLo())/2.);
        out.back()->SetBinError(i,(v->getErrorHi()-v->getErrorLo())/2.);
      }
    }
  }
 
  
  return out;
}


void TRooUnfold::Print(Option_t* opt) const {
  TNamed::Print(opt);
  std::cout << "  # Truth bins = " << m_truthVar->numBins() << std::endl;
  std::cout << "  # Reco bins  = " << m_recoVar->numBins() << std::endl;
  std::cout << std::endl;
  std::cout << "  Model " << ((m_builtModel) ? "IS" : "is NOT") << " Built" << std::endl;

  if(m_builtModel) {
    //print the model, with all the normalization and scale factor terms ...
    std::cout << "   Backgrounds: ";
    bool doneFirst(false);
    for(auto b : m_bkgOrder) {
      if(b.BeginsWith("signal_bin")) continue;
      if(doneFirst) std::cout << ", ";
      doneFirst=true;
      std::cout << b;
    }
    std::cout << std::endl;
    
    std::cout << "   Scale/Normalization Factors: ... (multiplies) ..." << std::endl;
    
    for(auto sf : m_sfApplied) {
      std::cout << "    " << sf.first << " ";
      if(m_recoSF.find(sf.first)!=m_recoSF.end()) std::cout << "(S)";
      else std::cout << "(N)";
      std::cout << ": ";
      bool doneFirst(false);
      for(auto s : sf.second) {
        if(doneFirst) std::cout << ", ";
        doneFirst=true;
        std::cout << s;
      }
      std::cout << std::endl;
    }
    
    std::cout << "   Systematics: ... (affects) ..." << std::endl;
    
    for(auto np : m_systNP) {
      std::cout << "    " << np.first << " (" << np.second->getStringAttribute("constraintType") << "):";
      bool doneFirst(false);
      for(auto sf : m_sfFunctions) {
        if(doneFirst) std::cout << ", ";
        doneFirst=true;
        if(sf.second->getParameter(np.first)) std::cout << sf.first;
      }
      
      for(auto sf : m_bkgPdfs) {
        if(doneFirst) std::cout << ", ";
        doneFirst=true;
        if(sf.second->getParameter(np.first)) std::cout << sf.first;
      }
      
      std::cout << std::endl;
    }
    
  }

}