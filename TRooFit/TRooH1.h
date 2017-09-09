/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOH1
#define TROOH1





//need this horrible hack to unset _owner in cachemanager when calling getObj and setObj
//FIXME: should also protect calls made in getNormObj method (method is virtual, so can override)
#define protected public 
#include "RooObjCacheManager.h" 
#undef protected

#include "RooProdPdf.h"

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooSetProxy.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsCategory.h"
#include "TRooFit/TRooFitResult.h"

#include "RooAbsData.h"

#include "TH1.h"

#include "TVirtualPad.h"


class RooProdPdf;
class TRooH1;

class TRooAbsH1 : public TAttLine, public TAttFill, public TAttMarker {
public:
  friend class TRooAbsHStack;

  TRooAbsH1() { UseCurrentStyle(); }

  TRooAbsH1(const RooArgList& observables, RooAbsArg* me);

  TRooAbsH1(const TRooAbsH1& other, RooAbsArg* me);
  inline virtual ~TRooAbsH1() { SafeDelete(fDrawHistogram); SafeDelete(fThisWithConstraints); }
  
  
  //derived classes must implement this method
  virtual const char* GetName() const = 0; 
  //derived classes must implement this method
  virtual const char* GetTitle() const = 0; 
  virtual Double_t getVal(const RooArgSet* nset = 0) const = 0;
  virtual Double_t getVal(const RooArgSet& nset) const = 0;
  virtual TIterator* clientIterator() const = 0;
  virtual RooAbsArg* cloneTree(const char* newname=0) const = 0;
  virtual RooArgSet* getDependents(const RooArgSet& set) const = 0;
  virtual RooArgSet* getParams(const RooArgSet& set) const = 0;
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const = 0;
  virtual Double_t expectedEvents(const RooArgSet&) const = 0;
  virtual TH1* getNominalHist() const = 0; //retrieves the underlying nominal histogram
  
  virtual Double_t missingEvents() const; //default implementation just checks fMissingBin - override in stacks to combine components
  
  void UseCurrentStyle();
  
  Int_t GetDimension() const { return fObservables.getSize(); }
  
  
  bool addNormFactor( RooAbsReal& factor ); //add a norm factor
  bool addShapeFactor( int bin, RooAbsReal& factor ); //add a shape factor to a bin
  bool addShapeFactor( const char* name, RooAbsReal& factor ); //add shape to a category bin

  inline void SetRangeName(const char* name) { fRangeName = name; }
  virtual const char* GetRangeName(const char* name=0) const { 
    if(name) return name;
    if(fRangeName.Length()) return fRangeName.Data(); 
    return GetName(); 
  } //derived classes will override this
  
  //these functions use the ACTUAL binning of this object
  Int_t FindFixBin( double x ) const;
  std::unique_ptr<RooArgSet> GetShapeFactors(int bin) const; 
  
  
  Int_t getBin(const char* rangeName=0) const; //return the current bin number, for the given binning
  
  //the next lot of functions rely on what is returned by GetRangeName to choose the binning
  Double_t getBinVolume() const; //return volume of current bin
  Double_t getError(const RooFitResult& fr) const;
  Double_t getBinError(const RooFitResult& fr) const { return getBinVolume()*getError(fr); }
  inline Double_t getBinContent(const RooArgSet* nset=0) const { return getBinVolume()*getVal(nset)*(expectedEvents(nset) /*+ missingEvents()*/); } //FIXME: assumes getVal is
  inline Double_t getBinContent(const RooArgSet& nset) const { return getBinContent(&nset); }   //  flat across the bin!
  RooRealVar* getStatFactor(int bin, bool createIf=false);
  
  //the following functions temporarily move the observables to the given bin, then call a method above
  Double_t GetBinContent(int bin,const RooFitResult* fr = 0) const;
  Double_t GetBinContent(int bin,const TRooFitResult& fr) const { return GetBinContent(bin,&fr); }
  Double_t GetBinContent(const char* bin, const RooFitResult* fr = 0) const;
  inline Double_t GetBinContent(const char* bin,const TRooFitResult& fr) const { return GetBinContent(bin,&fr); }
  Double_t GetBinError(int bin, const RooFitResult* fr = 0) const;
  inline Double_t GetBinError(int bin, const TRooFitResult& fr) const { return GetBinError(bin,&fr); }
  TH1* GetHistogram(const RooFitResult* fr = 0, bool includeErrors=false, TH1* histToFill=0) const; //fills a histogram with the values (and errors) corresponding to the fit result
  TH1* GetHistogram(const TRooFitResult& fr) const { return GetHistogram(&fr); }
  void fillHistogram(TH1* histToFill, const RooFitResult* r, bool includeErrors) const;
  void fillGraph(TGraph* graphToFill, const RooFitResult* r, bool includeErrors, int nPoints=100) const;
  
  
  //these functions are used when fitting
  RooAbsPdf&  model(); //returns self, with constraint terms added if necessary
  RooProdPdf* buildConstraints(const RooArgSet& obs, const char* systGroups="", bool addSelf=false) const;
  
  virtual void Draw(Option_t* option = "");
  virtual void Draw(const TRooFitResult& r, Option_t* option = "");
  virtual void Paint(Option_t* option = "");
  
  //other functions
  void setRooFitValV(bool in) { kUseAbsPdfValV = in; }//option to fall back to RooFit's usual evaluation
  virtual void setFloor(bool in) { kMustBePositive = in; } //if true, 'pdf' evaluations cannot go negative
  
protected:
  TH1* createOrAdjustHistogram(TH1* hist, bool noBinLabels=false) const; //rebins and styles the given histogram (creating it if no hist given

  RooListProxy fObservables; //observables define bins ... *must* be defined at construction
  RooListProxy fNormFactors; //these multiply all bins the same
  RooListProxy fShapeFactors; //these multiply a single bin
  
  RooListProxy fStatFactors; //these are a type of shape factor, but here we maintain an owning list
  
  std::map<int,std::vector<int>> fBinsShapeFactors; //maps a bin onto the shape factors it is affected by
  
  struct DrawnHistogram {
    TVirtualPad* pad = 0; //which pad the hist is drawn on
    TObject* hist = 0; //the histogram or graph
    TObject* postHist = 0; //a secondary hist or graph, drawn after the first
    TRooFitResult* fr = 0; //the associated fit result for the parameter snapshot used to fill the histogram (optional)
    TString opt; //the draw option
    TString postHistOpt; //the posthist draw option
  };
  std::vector<DrawnHistogram> fDrawHistograms; //!
  
  mutable TH1* fDrawHistogram = 0; //our histogram that we update and stuff when requests to 'draw' us
  //TRooFitResult* fDrawFitResult = 0;
  
  
  TString fRangeName; //when 'getting' bin contents, will use this binning
  
  TH1* fDummyHist = 0; //currently used to check binning consistencies etc
  
  bool kUseAbsPdfValV=false; //if we should just use the RooAbsPdf getValV method;
  bool kMustBePositive=false; //if true, return in getValV forced to 0 (if integral is <=0 and val is <=0 then we return "1" as value
  
  
  TRooH1* fMissingBin = 0; //when FillMissing call, created a TH0D to hold this contribution to pdf normalization
  RooRealProxy fMissingBinProxy; //hold a proxy to it
  
  RooProdPdf* fThisWithConstraints = 0; //constructed in 'model' method.

private:
    ClassDef(TRooAbsH1,1) // The Abstract Base class for all TRooFit pdfs
};


class TRooH1 : public RooAbsPdf, public TRooAbsH1 {
public:
  friend class TRooHStack;
  
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAbsPdf::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAbsPdf::getVal(nset); }
  TIterator* clientIterator() const { return RooAbsPdf::clientIterator(); }
  RooAbsArg* cloneTree(const char* newname=0) const { return RooAbsPdf::cloneTree(newname); }
  RooArgSet* getDependents(const RooArgSet& set) const { return RooAbsPdf::getDependents(set); }
  RooArgSet* getParams(const RooArgSet& set) const { return RooAbsPdf::getParameters(set); }
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset) ; }
  virtual TH1* getNominalHist() const { return (fHists.size()) ? fHists[0] : 0; }
  ///----
  
  TRooH1() {} ; 
  
  //FIXME need destructor to delete fHists
  
  TRooH1(const char* name, const char* title); //single bin
  
  TRooH1(const char *name, const char *title,const RooArgList& observables,TH1* hist);
  TRooH1(const char *name, const char *title,
	      const RooArgList& observables, const int* bins, const double* min, const double* max );
  
  TRooH1(const char *name, const char *title,
	      const RooArgList& observables, std::vector<int>&& bins, std::vector<double>&& min, std::vector<double>&& max );

  TRooH1(const char *name, const char *title,
	      const RooArgList& observables, std::vector<int>&& bins, std::vector<Double_t*>&& binEdges );
  
  TRooH1(const TRooH1& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new TRooH1(*this,newname); }
  

  virtual void Paint(Option_t* option = "") { TRooAbsH1::Paint(option); }
  using TRooAbsH1::Draw;

  bool addParameter( RooAbsArg& arg ); //current value is taken as value for all existing data
  

  Bool_t Add( const TH1* h1 , Double_t c1 = 1);

  Int_t Fill( RooAbsReal& val ); 
  Int_t Fill( const char* name , double w = 1. );
  Int_t Fill( double x , double w =1.); 
  
  //virtual void SetBinContent( int bin, RooAbsReal& val ); //turns bin value into a function .. not available yet
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
