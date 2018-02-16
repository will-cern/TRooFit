

#include "TRooFit/apps/TRooMatrixMethod.h"

#include "TRooFit/Utils.h"

#include "TLatex.h"
#include "TCanvas.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

TRooMatrixMethod::TRooMatrixMethod(const char* name, const char* title, uint nObjects, TH1* realEfficiency, TH1* fakeEfficiency ) : 
  TNamed(name,title), m_nObjects(nObjects) {
  
  //FIXME: check real and fake rate histograms are compatible binning
  
  m_efficiencyDummyHist = (TH1*)realEfficiency->Clone("dummy");
  m_efficiencyDummyHist->SetDirectory(0);
  m_efficiencyDummyHist->Reset();
  
  
  ///1. create a RooCategory for the 2^nObjects regions and a RooRealVar for the nBins^nObj bins per region
  
  int nEffBins = realEfficiency->GetNbinsX()*realEfficiency->GetNbinsY();
  int nBinsTot = std::pow(nEffBins,nObjects);
  
  
  //when doing the averaging method, will need one TRooH1D for each term of the matrix
  //i.e. will need nRegions^2 terms
  bool doAveraging=false;
  
  if(nBinsTot>10) {
    doAveraging=true;
    m_doAveraging=true;
    nBinsTot = 1;
  }
  
  int nRegions = std::pow(2,nObjects);
  
  m_xVar = new RooRealVar("x","BinNumber",1.-0.5,nBinsTot+1.-0.5);
  m_cat = new RooCategory("region","region");
  for(int i=0;i<nRegions;i++) {
    std::string label;
    for(uint j=0;j<nObjects;j++) label += ( (i&(1<<j)) ? "T" : "t" ); //makes labels = ttt, Ttt, tTt, TTt, ttT, TtT, tTT, TTT
    m_cat->defineType(label.c_str());
  }
  
  
  ///2. Create the 2^nObject (=nRegions) truth-level components, with shapeFactors in each bin so that whole thing can float
  
  for(int i=0;i<nRegions;i++) {
    std::string label;
    for(uint j=0;j<nObjects;j++) label += ( (i&(1<<j)) ? "F" : "R" ) ; //makes first label = RRR, FRR, RFR, ...
  
    m_truthComps[i] = new TRooH1D(label.c_str(),label.c_str(),*m_xVar,nBinsTot);
    
    for(int j=1;j<=nBinsTot;j++) {
       RooRealVar* sf = new RooRealVar(Form("sf_%s_bin%d",label.c_str(),j),
                                        Form("sf_%s_bin%d",label.c_str(),j),1.,-1e-9,1); //scale factor range is between just under 0 and 1
        m_truthComps[i]->SetBinContent(j,2);
        m_truthComps[i]->addShapeFactor(j,*sf);
     }
  }
  
  if(!doAveraging) {
  
    ///3. Create individual TRooHF1D for each bin of the real and fake efficiencies
    /// These will be used to build up efficiency histograms for each object
    /// r1 repeats every N bins, r2 events N*N bins, and so on 
    
    
    int k=0;
    for(int i=1;i<=realEfficiency->GetNbinsX();i++) {
      for(int j=1;j<=realEfficiency->GetNbinsY();j++) {
        k++;
        TRooHF1* binFunc = new TRooHF1(Form("r_bin%d",k),Form("r_{%d}",k));
        binFunc->SetBinContent(1,realEfficiency->GetBinContent(i,j));
        if(realEfficiency->GetBinError(i,j)) binFunc->SetBinError(1,realEfficiency->GetBinError(i,j));
        m_realBins.push_back(binFunc);
        
        binFunc = new TRooHF1(Form("f_bin%d",k),Form("f_{%d}",k));
        binFunc->SetBinContent(1,fakeEfficiency->GetBinContent(i,j));
        if(fakeEfficiency->GetBinError(i,j)) binFunc->SetBinError(1,fakeEfficiency->GetBinError(i,j));
        m_fakeBins.push_back(binFunc);
      }
    }
    
    
    for(uint k=0;k<nObjects;k++) {
      m_realEff[k] = new TRooHF1D(Form("r%d",k+1),"Real Efficiency",*m_xVar,nBinsTot);
      m_fakeEff[k] = new TRooHF1D(Form("f%d",k+1),"Fake Efficiency",*m_xVar,nBinsTot);
      
      int currBin = 0;
      
      int repeatLength = std::pow(nEffBins,k);
      for(int i=0;i<nBinsTot;) {
        
        for(int j=0;j<repeatLength;j++) {
          m_realEff[k]->Fill( m_realEff[k]->GetXaxis()->GetBinCenter(i+1), *m_realBins[currBin] );
          m_fakeEff[k]->Fill( m_fakeEff[k]->GetXaxis()->GetBinCenter(i+1), *m_fakeBins[currBin] );
          i++;
        }
        currBin++;
        if(currBin==nEffBins) currBin=0;
      }
      
    }
  } else {
    //averaging method, create a single TRooHF1D for each term of the matrix
    for(uint i=0;i<nRegions;i++) { //the reco-level region
      
      for(uint j=0;j<nRegions;j++) { //the truth-level region
        std::string label="eff_";
        std::string label2="<";
        for(uint k=0;k<nObjects;k++) {
          if(j&(1<<k)) { //is F
            if(i&(1<<k)) { //is anti-tight (T)
              label += "F";
              label2 += Form("(1-f%d)",k);
            } else {
              label += "f";
              label2 += Form("f%d",k);
            }
          } else { //is R
            if(i&(1<<k)) { //is anti-tight (T)
              label += "R";
              label2 += Form("(1-r%d)",k);
            } else {
              label += "r";
              label2 += Form("r%d",k);
            }
          }
        }
      
        m_realEff[i*nRegions+j] = new TRooHF1D(label.c_str(),(label2+">").c_str(),*m_xVar,nBinsTot);
        
        
      }
    }
  }
  
  
  
  ///4. Construct the stacks for each region ... each stack will have a scaled version of each truth-level component
  /// e.g. the first region (ttt) will have RRR scaled by r1r2r3 and RRF scaled by r1r2f3 and so on
  
  for(uint i=0;i<nRegions;i++) {
    m_stacks[i] = new TRooHStack( Form("stack_%s",m_cat->lookupType(i)->GetName()),m_cat->lookupType(i)->GetName() );
    m_stacks[i]->SetMinimum(1e-9);
    for(int j=0;j<nRegions;j++) { //loop over components (RRR, etc)
      //for each object determine if it is real or fake and tight or anti-tight, use to choose which coefficients to use
      
      TString factorFormula;
      RooArgList terms;
      
      std::string label;
      
      for(uint k=0;k<nObjects;k++) {
        bool isReal = !(j&(1<<k));
        bool isTight = !(i&(1<<k));
        
        std::string f;
        if(isReal && isTight) f = Form("r%d",k+1);
        else if(isReal && !isTight) f = Form("(1-r%d)",k+1);
        else if(!isReal && isTight) f = Form("f%d",k+1);
        else f = Form("(1-f%d)",k+1);
        
        if(!doAveraging) {
          if(k>0) factorFormula += "*";
          factorFormula += f.c_str();
          if(isReal) terms.add( *m_realEff[k] );
          else terms.add( *m_fakeEff[k] );
        }
        
        label += (isReal) ? "R" : "F";
        
      }
      
      TRooH1D* hist = new TRooH1D(Form("h_%s_%s",m_cat->lookupType(i)->GetName(),label.c_str()),Form("h_%s_%s",m_cat->lookupType(i)->GetName(),label.c_str()),*m_xVar,nBinsTot);
      hist->Add(*m_truthComps[j]);
      if(!doAveraging) {
        RooFormulaVar* c = new RooFormulaVar(Form("c_%s",hist->GetName()),factorFormula.Data(),terms);
        hist->Scale(*c);
      } else {
        hist->Scale(*m_realEff[i*nRegions+j]);
      }
      
      hist->SetFillColor(j+50);
      m_stacks[i]->Add(hist);
      
      m_comps[i][j] = hist;
      
      hist->setFloor(true);
      
    }
    
  }
  
  ///5. Create histogram for each region, which will hold the data (used just before fitting to create a RooDataHist)
  
  for(int i=0;i<nRegions;i++) {
    m_dataHist[i] = new TH1D(Form("data_%s",m_cat->lookupType(i)->GetName()),
                             Form("Data in region %s",m_cat->lookupType(i)->GetName()),
                             nBinsTot,0.5,0.5+nBinsTot);
    m_dataHist[i]->SetDirectory(0);m_dataHist[i]->Sumw2();
    m_dataHist[i]->SetMarkerStyle(20);
  }
}

Int_t TRooMatrixMethod::Fill(int region, int bin) {

  if(m_doAveraging) {
    //need to lookup the right efficiencies for the region, and then fill the product of efficiencies into the right m_realEff
  }

  return m_dataHist[region]->Fill(m_dataHist[region]->GetBinCenter(bin));
}

TRooFitResult* TRooMatrixMethod::Fit() {

  int nRegions = std::pow(2,m_nObjects);

  ///Use the data hists to adjust the values of the TRooH1D components
  ///Start each truth comp off at 2*the data value of the 'corresponding' reco-level region
  ///e.g. RR starts off equal to 2*Data in the tt region
  for(int i=0;i<nRegions;i++) {
      for(int j=1;j<=m_dataHist[i]->GetNbinsX();j++) {
        if(m_truthComps[i]->GetBinContent(j) < 2.*m_dataHist[i]->GetBinContent(j)) m_truthComps[i]->SetBinContent(j,2.*m_dataHist[i]->GetBinContent(j));
      }
  }
  ///build the RooDataHist from the histograms
  
  std::map<std::string,TH1*> input;
  for(int i=0;i<nRegions;i++) {
    input[m_cat->lookupType(i)->GetName()] = m_dataHist[i];
  }
  
  m_data = new RooDataHist("data","data",RooArgList(*m_xVar),*m_cat,input);
  
  ///Build the models for each stack and put into a simultaneous
  
  if(m_model) delete m_model;
  m_model = new RooSimultaneous("model","model",*m_cat);
  for(int i=0;i<nRegions;i++) {
    auto modelWithConstraints = TRooFit::BuildModel(*m_stacks[i],*m_data);
    m_model->addPdf( *modelWithConstraints , m_cat->lookupType(i)->GetName() );
  }

  ///Draw prefit distributions ...
  TCanvas* cc = new TCanvas(Form("%s_prefit",GetName()),Form("%s Pre-fit",GetTitle()),800,600);
  int ii = sqrt(nRegions);
  cc->Divide(nRegions/ii,ii);
  
  for(int i=0;i<nRegions;i++) {
    cc->cd(i+1);
    if(m_lastFitResult) m_stacks[i]->Draw("e3005",m_lastFitResult); 
    else  m_stacks[i]->Draw("e3005"); 
    m_dataHist[i]->Draw("same");
  }
  
  cc->Modified(1);
  cc->Update();

  return 0;

  m_lastFitResult = new TRooFitResult(m_model->fitTo(*m_data,RooFit::Save()));
  
  ///Draw postfit distributions ...
  TCanvas* cc2 = new TCanvas(Form("%s_postfit",GetName()),Form("%s Post-fit",GetTitle()),800,600);
  cc2->Divide(nRegions/ii,ii);
  
  for(int i=0;i<nRegions;i++) {
    cc2->cd(i+1);
    m_stacks[i]->Draw("e3005",m_lastFitResult); 
    m_dataHist[i]->Draw("same");
  }
  cc2->Modified(1);
  cc2->Update();

  //print the background in the signal region (region 0)
  TRooHStack allFake("allFake","allFake");
  for(int i=1;i<nRegions;i++) {
    allFake.Add(m_comps[0][i]);
  }
  
  double bkgErr(0);
  double bkg = allFake.IntegralAndError(bkgErr);
  std::cout << "Uncorrelated error = " << bkgErr << std::endl;
  bkg = allFake.IntegralAndError(bkgErr,m_lastFitResult); //accounts for correlations between the free parameters
  
  std::cout << "Background prediction = " << bkg << " +/- " << bkgErr << std::endl;

  return m_lastFitResult;

}
