/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOH1
#define TROOH1


#include "TRooFit/TRooAbsH1.h"


class TRooH1 : public RooAbsPdf, public TRooAbsH1 {
public:
  friend class TRooHStack;
  
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAbsPdf::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAbsPdf::getVal(nset); }
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset) ; }
  virtual TH1* getNominalHist() const { return (fHists.size()) ? fHists[0] : 0; }
  ///----
  
  TRooH1() { RooMsgService::instance().getStream(RooFit::INFO).removeTopic(RooFit::NumIntegration); } ; 
  
  //FIXME need destructor to delete fHists
  
  TRooH1(const char* name, const char* title); //single bin
  TRooH1(const char *name, const char *title,const RooArgList& observables,TH1* hist);
  TRooH1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max );
  TRooH1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<double>&& min, std::vector<double>&& max );
  TRooH1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<const Double_t*>&& binEdges );
  TRooH1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max, std::vector<const Double_t*>&& binEdges, TH1* hist=0);
  
  TRooH1(const TRooH1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TRooH1(*this,newname); }
  

  virtual void Paint(Option_t* option = "") { TRooAbsH1::Paint(option); }
  using TRooAbsH1::Draw;

  bool addParameter( RooAbsArg& arg , int interpCode = 0 ); //current value is taken as value for all existing data
  

  Bool_t Add( const TH1* h1 , Double_t c1 = 1);

  Int_t Fill( double x, RooAbsReal& val );
  Int_t Fill( RooAbsReal& val ); 
  Int_t Fill( const char* name , double w = 1. );
  Int_t Fill( double x , double w =1.); 
  
  virtual void SetBinContent( int bin, RooAbsReal& val ); //turns bin value into a function 
  virtual void SetBinContent( const char* name, double val ); //set category bin
  virtual void SetBinContent( int bin, double val);
  virtual void SetBinError( int bin, double val); //will update the stat uncert

  
  //override getValV so we can suppress warnings about 0 and negative values
  virtual Double_t getValV( const RooArgSet* set = 0 ) const;
  

  
  

  TH1* GetHist(unsigned int paramSet) const { 
    //Returns the raw histogram associated to the given parameter spacepoint (paramSet)
    //You can check the current parameter values put you on spacepoint with TRooH1::getParamSet
    if(paramSet>=fHists.size()) return 0; return fHists[paramSet]; 
    
  }

  
  
  Int_t getParamSet() const; //the index of the current paramset
  

  bool setInterpCode(const char* parName, int code);
  bool setInterpCode(const RooAbsArg& arg, int code);
  
  Int_t getObsInterpCode() const { return fObsInterpCode; }
  
  TRooH1* createTransFactor( TRooH1* transferFrom );
  
  

  virtual ExtendMode extendMode() const { return CanBeExtended ; }
  virtual RooAbsArg::CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised; } //had to disable cache and track, otherwise fits not work
  virtual Bool_t syncNormalization(const RooArgSet* nset, Bool_t adjustProxies=kTRUE) const;

  const std::vector<double>& GetParamSet(int idx) const;
  
  
  void isData(bool forceBinned=false); //activates 'data' mode for histogram ... all fills will go into a RooAbsData (if we already have fills, then RooDataHist, else RooDataSet);
  RooAbsData* data() const { return fData; }


  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlow, Double_t xhi) const;
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const;

  void FillMissing(double w=1.); //adds a missing event ... this wont contribute to expectedEvents (or integral) but will affect getVal(x) results
  void SetMissingContent(double w);
  void SetMissingError(double w);

  //we will need to do some trickery if a NormFactor or ShapeFactor depends on our observables!
  //virtual Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet* normSet, const char* rangeName) const; 
  //virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const;

protected:
  Int_t getOrCreateParamSet();

  RooListProxy fParameters; //parameters are not integrated over
  RooListProxy fValues; //bins can take a functional value .. the functions are kept in this list
  

  Double_t evaluate() const ;

  std::vector<TH1*> fHists; //if 1 or 2D distribution we just use this
  //THnSparseD* fSparseValues = 0;  in future may extend to this for 3D or more
  
  std::map<int,std::vector<int>> fFunctionalBinValues; //the bins that have a function as their value (first int is bin num, second is the fValues entry)
  
  std::vector<std::vector<double>> fParameterSnapshots; //the specific parameter points we have data for
  
  TRooH1* fTransFactor = 0; //for now, only one transfer factor allowed .. in future, may allow for partial transfer factors 
  bool kIsTransNumerator=false;
  
  std::vector<Int_t> fInterpCode; //one for each parameter //0 = piecewise linear, 4 = 6th order polynomal with logarithmic extrapolation
  Int_t fObsInterpCode = 0; //bitwise flag to interpolate for given observable: FIXME: will limit you to 32 observables

  RooAbsData* fData = 0; //gets set up by setData, and then subsequent fills go into this
  RooRealVar* fDataWeightVar = 0; //used when fData is a RooDataSet

  
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
