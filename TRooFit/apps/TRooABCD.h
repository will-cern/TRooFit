
#ifndef TROOABCD
#define TROOABCD

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1D.h"
#include "TRooFit/TRooHF1D.h"
#include "TNamed.h"
#include "RooDataSet.h"
#include "RooCategory.h"

#include "RooSimultaneous.h"


class TRooABCD : public TNamed  {

  public:
    TRooABCD() : TNamed() { }
    TRooABCD(const char* name, const char* title);
    
    
    Bool_t AddSignal(int region, TH1* signal);
    Bool_t AddData(int region, TH1* data);
    Bool_t AddOther(int region, TH1* other);
    
    RooRealVar* AddBkgScaleFactor(int region, double value, double uncert);
    RooRealVar* AddSignalScaleFactor(int region, double value, double uncert);
    RooRealVar* AddOtherScaleFactor(int region, double value, double uncert);
    
    void AddBkgScaleFactor(int region, RooRealVar* sf);
    void AddSignalScaleFactor(int region, RooRealVar* sf);
    void AddOtherScaleFactor(int region, RooRealVar* sf);
    
    TRooFitResult* Fit(int modelType=0, bool floatSignal=true);
    
  
    ClassDef(TRooABCD,1);

    Double_t GetBkgBinContent(int region, int bin);
    Double_t GetDataBinContent(int region, int bin);
    const TH1* GetDataHistogram(int region);

    Double_t GetBkgIntegral(int region, TRooFitResult* r=0);
    Double_t GetBkgIntegralError(int region, TRooFitResult* r=0);
    Double_t GetSignalIntegral(int region, TRooFitResult* r=0);
    Double_t GetSignalIntegralError(int region, TRooFitResult* r=0);

    RooDataSet* GetDataSet() { return m_data; }
    TRooH1* GetBkgHistogram(int region) { return m_bkg[region]; }

    RooSimultaneous* GetModel() { return m_model; }

    void SetPrintLevel(RooFit::MsgLevel level) { m_printLevel = level; } //use RooFit::INFO etc etc

    void SetSignalStrength(double in) { m_signalStrength->setVal(in); }


  private:
    void checkRangeChange(TH1* hist);
  
    RooRealVar* m_xVar = 0;
    RooRealVar* m_weightVar = 0;
    RooCategory* m_cat = 0;
    RooDataSet* m_data = 0;
    RooRealVar* m_signalStrength = 0;
  
    std::map<int, TRooH1D*> m_signal;
    std::map<int, TRooH1D*> m_bkg;
    std::map<int, TRooH1D*> m_other;
    
    std::map<int, TRooHStack*> m_stacks;
    
    RooSimultaneous* m_model = 0;

    std::map<int, TH1*> m_dataHist; //data histograms, constructed from m_data

    std::vector<RooRealVar*> m_modelPars;

    std::vector<TRooHF1D*> m_transferFactors; //these will be normFactors on m_bkg in region A and B

    TRooFitResult* m_lastFitResult = 0;

    RooFit::MsgLevel m_printLevel = RooFit::ERROR;

};

#endif