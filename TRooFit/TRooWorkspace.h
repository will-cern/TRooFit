/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * Modified version of a 
 * RooWorkspace ... Helps with model building                          * 
 *****************************************************************************/

#ifndef TROOWORKSPACE
#define TROOWORKSPACE

//next lines to hack RooFitResult so that statusHistory can be presered when fitTo
#define protected public
#include "RooFitResult.h"
#undef protected

#include "RooWorkspace.h"

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1.h"
#include "TRooFit/TRooHF1.h"

#include "TTree.h"
#include "TLegend.h"

#include "RooSimultaneous.h"
 
class TRooWorkspace : public RooWorkspace {
public:
  //using RooWorkspace::RooWorkspace; 
  
  TRooWorkspace() : RooWorkspace() { }
  TRooWorkspace(const char* name, const char* title = 0) : RooWorkspace(name,title) { }
  TRooWorkspace(const RooWorkspace& other);
  ~TRooWorkspace() { fNll.removeAll(); }
  

  bool definePoi(const char* poi) { return defineSet("poi",poi); }
  bool addArgument(const char* name, const char* title, double val);
  RooRealVar* addParameter(const char* name, const char* title, double val, double min, double max, const char* constraintType=0);
  RooRealVar* addParameter(const char* name, const char* title, double min, double max, const char* constraintType=0);
  RooRealVar* addObservable(const char* name, const char* title, double min, double max);
  RooRealVar* addObservable(const char* name, const char* title, int nBins, double min, double max);
  RooRealVar* addObservable(const char* name, const char* title, int nBins, const double* bins);
  TRooHStack* addChannel(const char* name, const char* title, const char* observable, int nBins, double min, double max);
  TRooHStack* addChannel(const char* name, const char* title, const char* observable, int nBins, const double* bins);
  TRooHStack* addChannel(const char* name, const char* title, const char* observable);
  
  bool addSamples(const char* name, const char* title, const char* channels="*", bool allowNegative=false);
  
  TRooHF1* addFactor(const char* name, const char* title, double nomVal=1.);
  bool SetFactorContent(const char* name,double val, const char* parName=0, double parVal=0 );
  
  bool Fill(const char* channel, double x, double w=1.); //fills data!
  Int_t Fill(const char* sample, const char* channel, double x, double w=1., const char* variationName=0, double variationVal=0);
  bool Fill(const char* sampleName, const char* channelNames, TTree* tree, const char* weight="1", const char* variationName=0, double variationVal=0);
  
  
  bool Add(const char* sample, const char* channel,  TH1* h1, const char* variationName=0, double variationVal=0);
  bool Add(const char* sample, const char* channel, RooAbsReal& arg);
  
  void SetBinContent(const char* sampleName, const char* channelName, Int_t bin, double val, const char* variationName=0, double variationVal=0);
  
  
  

  //add a normalization factor to a sample, across all channels
  void Scale(const char* sample, const char* channelNames, RooAbsReal& arg);
  void Scale(const char* sampleName,const char* channelNames, const char* par) { if(var(par)) Scale(sampleName,channelNames,*var(par)); }
  
  //set fill color of sample in given channels (semicolon separated list) (use "*" to specify all channels)
  void SetFillColor(const char* sample, const char* channelNames, Int_t in);
  void SetLineColor(const char* sample, const char* channelNames, Int_t in);
  
  double IntegralAndError(double& err, const char* sampleName, const char* channelName, const TRooFitResult& res) const;
  double IntegralAndError(double& err, const char* sampleName, const char* channelName) const { return IntegralAndError(err,sampleName,channelName,fCurrentFitResult?*fCurrentFitResult:""); }
  
  /*double sampleIntegralAndError(double& err, const char* sampleFullName, const TRooFitResult& res="") const;
  double sampleIntegralAndError(double& err, const char* channelName, unsigned int sampleNumber, const TRooFitResult& fr="") const;*/
  
  TRooH1* sample(const char* sampleName, const char* channelName);
  TRooAbsHStack* channel(const char* name) const;
  TRooHF1* factor(const char* factorName);
  
  void setData(const char* dataName) { 
    if(!data(dataName)) return;
    fCurrentData = dataName; 
  }
  
  using RooWorkspace::data;
  inline RooAbsData* data() { return RooWorkspace::data(fCurrentData); }
  
  void DisableForcedRecommendedOption(bool in) { kDisabledForcedRecommendedOptions=in; } //use to override forcing of the recommended fit options when calling fitTo
  
  double impact(const char* poi, const char* np, bool positive=true);
  void impact(const char* impactPar=0, float correlationThreshold=0);
  
  RooFitResult* fitTo(RooAbsData* theData, const RooArgSet* globalObserables=0, bool doHesse=true);
  RooFitResult* fitTo(const char* dataName=0, bool doHesse=true, const RooArgSet& minosPars=RooArgSet());
  RooFitResult* loadFit(const char* fitName,bool prefit=false);
  RooFitResult* getFit(const char* fitName=0) { return dynamic_cast<RooFitResult*>(obj((fitName==0)?fCurrentFit.Data():fitName)); }
  RooAbsReal* getFitNll(const char* fitName=0);
  
  RooAbsReal* nll(const char* nllName) const { return dynamic_cast<RooAbsReal*>(fNll.find(nllName)); }
  
  void DrawPLL(const char* parName, const char* opt="AL");
  
  double pll(const char* poi, RooAbsData* theData, const RooArgSet* globalObservables=0, bool oneSided=false, bool discovery=false);
  double pll(RooArgSet&& poi, RooAbsData* theData, const RooArgSet* globalObservables=0, bool oneSided=false, bool discovery=false);
  
  void addLabel(const char* label) { fLabels.push_back(label); }
  
  TLegend* GetLegend();
  
  double GetSampleCoefficient(const char* sampleFullName) const;
  
  RooSimultaneous* model(const char* channels="*");
  
  bool generateAsimov(const char* name, const char* title, bool fitToObsData=true);
  std::pair<RooAbsData*,RooArgSet*> generateToy(const char* name, const char* title, bool fitToObsData=true);
  
  
  //controls which channels are visible when drawing things
  //sets the 'hidden' attribute on non-visible channels
  Int_t SetVisibleChannels(const char* filter) { 
    setChannelAttribute("*","hidden",kTRUE); //hide all channels first
    return setChannelAttribute(filter,"hidden",kFALSE); //unhide the selected ones
  }
  Int_t setChannelAttribute(const char* channels,const char* attribute,Bool_t val=kTRUE);
  Int_t setVarAttribute(const char* vars,const char* attribute,Bool_t val=kTRUE);
  
  //draw a channel's stack and overlay the data too
  void channelDraw(const char* channel, Option_t* option="e3005", const TRooFitResult& res = "");
  
  //draws all channels
  virtual void Draw(Option_t* option, const TRooFitResult& res);
  virtual void Draw(Option_t* option="e3005") { 
    //if(fCurrentFit!="") Draw(option,getFit(fCurrentFit));
    if(fCurrentFitResult) Draw(option,*fCurrentFitResult);
    else Draw(option,""); 
  }
  
  //draws all channels, showing how values of channels depend on var
  void DrawDependence(const char* var, Option_t* option="TRI1");
  
  void SetRatioHeight(double in, bool showSignificance=false) { 
    if(fLegend) fLegend->SetTextSize( fLegend->GetTextSize() * (1. - fRatioHeight) / (1. - in) );
    fRatioHeight = in; 
    kShowSignificance=showSignificance;
  }
  
  Bool_t writeToFile(const char* fileName, Bool_t recreate=kTRUE);
  
  virtual void Print(Option_t* opt="") const;
  
  void FindVariations(double relThreshold, Option_t* opt = "samples");
  
  static void setDefaultStyle();
  
private:
  TString fChannelCatName = "channelCat"; //should only change if wrapping a non-conventional workspace
  TString fSimPdfName = "simPdf"; //should only change if wrapping a non-conventional workspace

  std::map<TString,TH1*> fDummyHists;
  
  TString fCurrentData = "obsData";
  
  TString fCurrentFit = "";
  Bool_t fCurrentFitIsPrefit = false; 
  TRooFitResult* fCurrentFitResult = 0;
  
  //RooArgList fStagedChannels; //channels cannot be added until they are frozen (all content filled)
  
  TLegend* fLegend = 0;

  std::vector<TString> fLabels; //plot labels

  Double_t fRatioHeight = 0; //if nonzero, channelDraw will draw ratio plots

  Bool_t fIsHFWorkspace = false; //if true, this is a histfactory workspace
  Bool_t kDisabledForcedRecommendedOptions = false;
  
  Bool_t kShowSignificance = false;
  
  RooArgList fNll; //!
  
  ClassDef(TRooWorkspace,1) // An extended form of a RooWorkspace
};
 
#endif
