
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
    
    void SetXTitle(const char* title) { m_xVar->SetTitle(title); }
    
    Bool_t AddSignal(int region, TH1* signal);
    Bool_t AddData(int region, TH1* data);
    Bool_t AddOther(int region, TH1* other);
    
    RooRealVar* AddBkgScaleFactor(int region, double value, double uncert);
    RooRealVar* AddSignalScaleFactor(int region, double value, double uncert);
    RooRealVar* AddOtherScaleFactor(int region, double value, double uncert);
    
    void AddBkgScaleFactor(int region, RooRealVar* sf);
    void AddSignalScaleFactor(int region, RooRealVar* sf);
    void AddOtherScaleFactor(int region, RooRealVar* sf);
    
    TRooFitResult* Fit(bool drawPostFit=true);
    
  
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

    

    void SetPrintLevel(RooFit::MsgLevel level) { m_printLevel = level; } //use RooFit::INFO etc etc

    //use this method to blind the SR from the fit, even if there is data!
    void SetValidationRegion(bool in) { m_isVR = in; }

    bool BuildModel(); //Finishes constructing the TRooFit model
    RooSimultaneous* GetModel() { return m_model; }

    RooRealVar* GetParameter(const char* name) {
      if(!m_model) BuildModel();
      RooArgSet s; m_model->treeNodeServerList(&s);
      return dynamic_cast<RooRealVar*>(s.find(name));
    }
    void PrintParameters() const { m_allParameters.Print("v"); }
    
    virtual void        Draw(Option_t *option="") { Draw(0,option); }
    void Draw(TRooFitResult* fr, const char* canvasName="");
    
    double GetChi2PValue(TRooFitResult* fr=0);
    double GetBinomialPValue(TRooFitResult* fr=0);

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

    TRooHF1D* m_transferFactor = 0; //this will be a normFactor on m_bkg in region A and B

    TRooFitResult* m_lastFitResult = 0;

    RooFit::MsgLevel m_printLevel = RooFit::ERROR;

    RooArgList m_allParameters; //!
    
    bool m_isVR = false;

};

#endif