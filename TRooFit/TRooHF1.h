/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOHF1
#define TROOHF1

#include "RooAbsReal.h"
#include "TRooFit/TRooAbsH1Fillable.h"



class TRooHF1 : public RooAbsReal, public TRooAbsH1Fillable {
public:

  TRooHF1() : RooAbsReal(), TRooAbsH1Fillable() { }
  TRooHF1(const char* name, const char* title) : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title) { }
  TRooHF1(const char *name, const char *title,const RooArgList& observables,TH1* hist) : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title,observables,hist) { 
    if(hist) Add( hist );
  }
  TRooHF1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max )  : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title, observables, bins, min, max) { }
  TRooHF1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<double>&& min, std::vector<double>&& max )  : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title, observables, std::move(bins), std::move(min), std::move(max)) { }
  TRooHF1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<const Double_t*>&& binEdges )  : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title, observables, std::move(bins), std::move(binEdges)) { }
  TRooHF1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max, std::vector<const Double_t*>&& binEdges, TH1* hist=0)  : RooAbsReal(name,title),TRooAbsH1Fillable(this,name,title, observables, bins, min, max,std::move(binEdges),hist) { 
    if(hist) Add( hist );
  }
  
  virtual void Paint(Option_t* option = "") { TRooAbsH1Fillable::Paint(option); }
  virtual void Draw(Option_t* option="") { TRooAbsH1Fillable::Draw(option); }
  virtual void Draw(Option_t* option,const TRooFitResult& r) { TRooAbsH1Fillable::Draw(option,r); }
  

  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAbsReal::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAbsReal::getVal(nset); }
  Double_t expectedEvents(const RooArgSet*) const { return 1.; }
  Double_t expectedEvents(const RooArgSet&) const { return 1.; }
  ///----
  
  Double_t getBinError(const RooFitResult& fr) const { return getError(fr); }
  Double_t getBinContent(const RooArgSet* nset=0) const { return getVal(nset); } //FIXME: assumes getVal is flat across the bin
  
  
  TRooHF1(const TRooHF1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TRooHF1(*this,newname); }
  


  RooAbsArg::CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised; } //had to disable cache and track, otherwise fits not work
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlow, Double_t xhi) const { return TRooAbsH1Fillable::binBoundaries(obs,xlow,xhi); }
  Bool_t isBinnedDistribution(const RooArgSet& obs) const { return TRooAbsH1Fillable::isBinnedDistribution(obs); }


protected:
  virtual RooAbsArg* me() { return this; }
  Double_t evaluate() const ;
  
private:

  ClassDef(TRooHF1,1) // The base class for all TRooFit histogram FUNCTIONS (not pdfs) (TRooHFxD)
};


#endif
