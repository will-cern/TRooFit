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
  
  virtual TH1* getNominalHist() const { return 0; } //FIXME ... do we need to return something better ever?
  
  virtual Double_t missingEvents() const;
  

  virtual void Paint(Option_t* option = "");

  void Draw(Option_t* option,const TRooFitResult& r);
  void Draw(Option_t* option = "");
  
  void Add(RooAbsReal* func);
  void Add(TRooH1* hist);
  THStack* GetStack(const RooFitResult* fr = 0) const;
  THStack* GetStack(const TRooFitResult& fr) const { return GetStack(&fr); }

  THStack* fillStack(THStack* stack, const RooFitResult* fr = 0, bool noRestyle=false) const;

  void SetMinimum(Double_t min) { fMinimum=min; if(fDrawStacks.size()) fDrawStacks.back().stack->SetMinimum(fMinimum); }
  void SetMaximum(Double_t max) { fMaximum=max; if(fDrawStacks.size()) fDrawStacks.back().stack->SetMaximum(fMaximum);}
  
  Double_t GetMinimum() const { return fMinimum; }
  Double_t GetMaximum() const { return fMaximum; }

  virtual TAxis* GetXaxis() const;


protected:
  
  
  TRooH1* firstHistogram();
  
  struct DrawnStack {
    TVirtualPad* pad = 0;
    THStack* stack = 0;
    TRooFitResult* fr = 0; //optional
    TString opt;
  };
  std::vector<DrawnStack> fDrawStacks; //!
  
  mutable THStack* fStack = 0; //histograms come from the TRooH1's, not deleted by us

  //this methods are implemented in the concrete classes
  virtual TIterator*& compIter() = 0;
  virtual RooListProxy& compList() = 0;virtual const RooListProxy& compList() const = 0;
  virtual void reinit() = 0; //called after adding new components
  //virtual TClass* IsA() const = 0;

  Double_t fMinimum = -1111;
  Double_t fMaximum = -1111;

public:
  ClassDef(TRooAbsHStack,1)

};
 
#endif
