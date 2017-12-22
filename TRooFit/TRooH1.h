/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOH1
#define TROOH1


#include "TRooFit/TRooAbsH1Fillable.h"


class TRooH1 : public RooAbsPdf, public TRooAbsH1Fillable {
public:
  friend class TRooHStack;
  
  TRooH1() : RooAbsPdf(), TRooAbsH1Fillable() { }
  TRooH1(const char* name, const char* title) : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title) { }
  TRooH1(const char *name, const char *title,const RooArgList& observables,TH1* hist) : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title,observables,hist) { 
    if(hist) Add( hist ); //adds it, so that stat factors are created ... have to do here so that me() is defined
  }
  TRooH1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max )  : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title, observables, bins, min, max) { }
  TRooH1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<double>&& min, std::vector<double>&& max )  : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title, observables, std::move(bins), std::move(min), std::move(max)) { }
  TRooH1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<const Double_t*>&& binEdges )  : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title, observables, std::move(bins), std::move(binEdges)) { }
  TRooH1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max, std::vector<const Double_t*>&& binEdges, TH1* hist=0)  : RooAbsPdf(name,title),TRooAbsH1Fillable(this,name,title, observables, bins, min, max,std::move(binEdges),hist) { 
    if(hist) Add( hist ); //adds it, so that stat factors are created ... have to do here so that me() is defined
  }
  
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAbsPdf::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAbsPdf::getVal(nset); }
  Double_t expectedEvents(const RooArgSet* nset=0) const;
  Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset); }
  ///----
  
  void resetNormMgr() { _normMgr.sterilize(); _norm=0; } //difference between reset and sterilize?
  
  virtual void Paint(Option_t* option = "") { TRooAbsH1Fillable::Paint(option); }
  virtual void Draw(Option_t* option="") { TRooAbsH1Fillable::Draw(option); }
  virtual void Draw(Option_t* option,const TRooFitResult& r) { TRooAbsH1Fillable::Draw(option,r); }
  
  TRooH1(const TRooH1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TRooH1(*this,newname); }
  

  
  //override getValV so we can suppress warnings about 0 and negative values
  virtual Double_t getValV( const RooArgSet* set = 0 ) const;
  

  ExtendMode extendMode() const { return CanBeExtended ; }
  RooAbsArg::CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised; } //had to disable cache and track, otherwise fits not work
  Bool_t syncNormalization(const RooArgSet* nset, Bool_t adjustProxies=kTRUE) const;


  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlow, Double_t xhi) const { return TRooAbsH1Fillable::binBoundaries(obs,xlow,xhi); }
  Bool_t isBinnedDistribution(const RooArgSet& obs) const { return TRooAbsH1Fillable::isBinnedDistribution(obs); }

  void FillMissing(double w=1.); //adds a missing event ... this wont contribute to expectedEvents (or integral) but will affect getVal(x) results
  void SetMissingContent(double w);
  void SetMissingError(double w);
  Double_t GetMissingContent() const { if(!fMissingBin) return 0; return fMissingBin->GetBinContent(1); }

  //we will need to do some trickery if a NormFactor or ShapeFactor depends on our observables!
  //virtual Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName) const; 
  //virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const;


protected:
  virtual RooAbsArg* me() { return this; }

  Double_t evaluate() const { return TRooAbsH1Fillable::evaluateImpl(true); }

  
private:

  ClassDef(TRooH1,1) // The base class for all TRooFit histograms (TRooHxD)
};


class TRooH0D : public TRooH1 {
public:
  TRooH0D() : TRooH1() { }
  
  TRooH0D(const char *name, const char *title, double val=0, double err=0) : TRooH1(name,title,RooArgList(),{},{}) {
    if(val) SetBinContent(1,val);
    if(err) SetBinError(1,err);
  }

  using TRooAbsH1::Draw;

  using TRooH1::SetBinContent;
  virtual inline void SetBinContent(double val) { SetBinContent(1,val); } //set value without creating a stat uncert

private:
  ClassDef(TRooH0D,1) //  A zero-bin (simple event count) TRooFit histogram [do not use because of RooFit problems]
};
 
#endif
