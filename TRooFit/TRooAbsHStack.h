/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * Base class for TRooHStack (RooRealSumPdf-like)
 * and TRooHPdfStack (RooAddPdf-like)                                        * 
 * Both implement 'beeston-barlow-lite' acquisition of statFactors:
 * will apply Beeston-Barlow-Lite principle to stat
 * uncerts ... i.e. all subhists will end up sharing
 * the same statFactors                                                      * 
 *****************************************************************************/

#ifndef TROOABSHSTACK
#define TROOABSHSTACK

#include "TRooFit/TRooH1.h"
#include "THStack.h"


//base class for TRooHStack (RooRealSumPdf-like) and TRooHPdfStack (RooAddPdf-like)
class TRooAbsHStack : public TRooAbsH1 {
public:
  using TRooAbsH1::TRooAbsH1;
  ~TRooAbsHStack();
  
  virtual Double_t getBinContent(const RooArgSet* nset) const;
  
  virtual TH1* getNominalHist() const { return 0; } //FIXME ... do we need to return something better ever?
  
  virtual Double_t missingEvents() const;
  

  virtual void Paint(Option_t* option = "");

  void Draw(Option_t* option,const TRooFitResult& r);
  void Draw(Option_t* option = "");
  
  void DrawDependence(const char* _var, Option_t* option=""); //use to show the value of this TRooAbsH1 as a function of the given var
  
  void Add(RooAbsReal* func, bool acquireStatFactors=true);
  THStack* GetStack(const RooFitResult* fr = 0) const;
  THStack* GetStack(const TRooFitResult& fr) const { return GetStack(&fr); }

  THStack* fillStack(THStack* stack, const RooFitResult* fr = 0, bool noRestyle=false) const;

  void SetMinimum(Double_t min) { fMinimum=min; if(fDrawStacks.size()) fDrawStacks.back().stack->SetMinimum(fMinimum); }
  void SetMaximum(Double_t max) { fMaximum=max; if(fDrawStacks.size()) fDrawStacks.back().stack->SetMaximum(fMaximum); }
  
  

  virtual TAxis* GetXaxis() const;
  virtual TAxis* GetYaxis() const;

  virtual void setBlindRange(const char* rangeName) { 
    RooFIter funcIter = compList().fwdIterator() ;
    RooAbsReal* func ;
    while((func=(RooAbsReal*)funcIter.next())) {
      TRooAbsH1* trooFunc = dynamic_cast<TRooAbsH1*>(func);
      if(trooFunc) {
        trooFunc->setBlindRange(rangeName);
      }
    }
    TRooAbsH1::setBlindRange(rangeName);
    
  }

  virtual RooListProxy& compList() = 0;
  virtual const RooListProxy& compList() const = 0;
  virtual const RooListProxy& coeffList() const = 0;

protected:
  
  
  TRooH1* firstHistogram();
  
  struct DrawnStack {
    TVirtualPad* pad = 0;
    THStack* stack = 0;
    TH1* frame = 0; //the histogram we use to draw the axis
    TRooFitResult* fr = 0; //optional
    TString opt;
  };
  std::vector<DrawnStack> fDrawStacks; //!
  
  mutable THStack* fStack = 0; //histograms come from the TRooH1's, not deleted by us

  //this methods are implemented in the concrete classes
  virtual TIterator*& compIter() = 0;
  
  
  virtual void reinit() = 0; //called after adding new components
  //virtual TClass* IsA() const = 0;


public:
  ClassDef(TRooAbsHStack,1)

};
 
#endif
