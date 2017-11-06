/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOHF1
#define TROOHF1

#include "TRooFit/TRooH1.h"



class TRooHF1 : public RooAbsReal, public TRooAbsH1 {
public:
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAbsReal::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAbsReal::getVal(nset); }
  virtual Double_t expectedEvents(const RooArgSet*) const { return 1.; }
  virtual Double_t expectedEvents(const RooArgSet&) const { return 1.; }
  virtual TH1* getNominalHist() const { return (fHists.size()) ? fHists[0] : 0; }
  ///----
  
  virtual Double_t getBinError(const RooFitResult& fr) const { return getError(fr); }
  inline virtual Double_t getBinContent(const RooArgSet* nset=0) const { return getVal(nset); } //FIXME: assumes getVal is flat across the bin
  
  
  TRooHF1() {} ; 
  
  //FIXME need destructor to delete fHists
  
  TRooHF1(const char* name, const char* title); //single bin
  TRooHF1(const char *name, const char *title,const RooArgList& observables,TH1* hist);
  TRooHF1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max );
  TRooHF1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<double>&& min, std::vector<double>&& max );
  TRooHF1(const char *name, const char *title,const RooArgList& observables, std::vector<int>&& bins, std::vector<const Double_t*>&& binEdges );
  TRooHF1(const char *name, const char *title,const RooArgList& observables, const int* bins, const double* min, const double* max, std::vector<const Double_t*>&& binEdges, TH1* hist=0);
  
  TRooHF1(const TRooHF1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TRooHF1(*this,newname); }
  

  virtual void Paint(Option_t* option = "") { TRooAbsH1::Paint(option); }
  virtual void Draw(Option_t* option="") { TRooAbsH1::Draw(option); }
  virtual void Draw(Option_t* option,const TRooFitResult& r) { TRooAbsH1::Draw(option,r); }
  //using TRooAbsH1::Draw; ... seems that this isn't enough, even though works for TRooH1 inside workspaces??

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

  

  
  

  TH1* GetHist(unsigned int paramSet) const { 
    //Returns the raw histogram associated to the given parameter spacepoint (paramSet)
    //You can check the current parameter values put you on spacepoint with TRooH1::getParamSet
    if(paramSet>=fHists.size()) return 0; return fHists[paramSet]; 
    
  }

  
  
  Int_t getParamSet() const; //the index of the current paramset
  

  bool setInterpCode(const char* parName, int code);
  bool setInterpCode(const RooAbsArg& arg, int code);
  
  Int_t getObsInterpCode() const { return fObsInterpCode; }
  
  

  virtual RooAbsArg::CacheMode canNodeBeCached() const { return RooAbsArg::NotAdvised; } //had to disable cache and track, otherwise fits not work

  const std::vector<double>& GetParamSet(int idx) const;
  

  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlow, Double_t xhi) const;
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const;


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
  
  
  std::vector<Int_t> fInterpCode; //one for each parameter //0 = piecewise linear, 4 = 6th order polynomal with logarithmic extrapolation
  Int_t fObsInterpCode = 0; //bitwise flag to interpolate for given observable: FIXME: will limit you to 32 observables


  
private:

  ClassDef(TRooHF1,1) // The base class for all TRooFit histogram FUNCTIONS (not pdfs) (TRooHFxD)
};


#endif
