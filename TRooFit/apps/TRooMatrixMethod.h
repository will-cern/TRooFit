
#ifndef TROOMATRIXMETHOD
#define TROOMATRIXMETHOD

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1D.h"
#include "TRooFit/TRooHF1D.h"
#include "TNamed.h"
#include "RooDataHist.h"
#include "RooCategory.h"

#include "RooSimultaneous.h"


class TRooMatrixMethod : public TNamed  {

  public:
    TRooMatrixMethod() : TNamed() { }
    TRooMatrixMethod(const char* name, const char* title, unsigned int nObjects, TH1* realEfficiency, TH1* fakeEfficiency);
    
    Int_t GetDimension() const {
      return m_efficiencyDummyHist->GetDimension();
    }
    
    Int_t FindFixRegion(const std::vector<bool>&& isTight) const {
      int out(0);
      for(uint i=0;i<isTight.size();i++) {
        if(!isTight[i]) out += (1<<i);
      }
      return out;
    }
    
    Int_t FindFixBin(const std::vector<double>&& vals) const {
      if(GetDimension()!=1) {
        Error("FindFixBin","Dimensionality of efficiency paramterization is %d, not 1",GetDimension());
        return -1;
      }
      if(vals.size()!=m_nObjects) {
        Error("FindFixBin","Number of values (%d) not equal to number of objects (%d)", int(vals.size()), m_nObjects);
        return -1;
      }
      int out = 0; int binFactor = 1; 
      for(uint i=0;i<m_nObjects;i++) {
        int bin = m_efficiencyDummyHist->FindFixBin(vals[i]);
        if(i!=0) bin--;
        bin *= binFactor;
        binFactor *= m_efficiencyDummyHist->GetNbinsX();
        out += bin;
      }
      return out;
    }
    
    Int_t Fill(const std::vector<bool>&& isTight, const std::vector<double>&& obs) {
      return Fill(std::move(isTight), FindFixBin(std::move(obs)));
    }
    
    Int_t FillBin(const std::vector<bool>&& isTight, const std::vector<int>&& bins) {
      int out = 0; int binFactor = 1;
      for(uint i=0;i<m_nObjects;i++) { 
        int bin = bins[i];
        if(i!=0) bin--;
        bin *= binFactor;
        binFactor *= m_efficiencyDummyHist->GetNbinsX() * m_efficiencyDummyHist->GetNbinsY();
        out += bin;
      }
      return Fill(std::move(isTight), out);
    }
    
    Int_t Fill(const std::vector<bool>&& isTight, int bin) {
      return Fill(FindFixRegion(std::move(isTight)),bin);
    }
    
    Int_t Fill(int region, int bin);  //add a data event to given region and bin
    
    TRooFitResult* Fit();
    
  
    ClassDef(TRooMatrixMethod,1);


    bool m_doAveraging = false; //in this mode, there is only one bin per region, and the coefficient is determined by averaging 
    
    unsigned int m_nObjects;
    TH1* m_efficiencyDummyHist = 0; //holds an emptied clone of one of the efficiency histograms
  
    RooRealVar* m_xVar = 0;
    RooRealVar* m_weightVar = 0;
    RooCategory* m_cat = 0;
    RooDataHist* m_data = 0;
  
    std::map<int, TRooH1D*> m_truthComps; //truth level distributions, i.e. N_RR, N_RF, etc
  
    std::map<int, std::map<int, TRooH1D* >> m_comps; //first index is region, second is comp number (nComps = nRegions)
    std::map<int, TRooHStack*> m_stacks; //one stack per region
    
    std::map<int, TRooHF1D*> m_realEff, m_fakeEff;
    std::vector<TRooHF1*> m_realBins, m_fakeBins; //bins of the real and fake efficiencies
    
    RooSimultaneous* m_model = 0;

    std::map<int, TH1*> m_dataHist; //data histograms, constructed from m_data

    TRooFitResult* m_lastFitResult = 0;

    RooFit::MsgLevel m_printLevel = RooFit::ERROR;

};

#endif