/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOABSH1
#define TROOABSH1





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
  
  
  //derived classes must implement these methods
  //------
  virtual Double_t getVal(const RooArgSet* nset = 0) const = 0;
  virtual Double_t getVal(const RooArgSet& nset) const = 0;
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const = 0;
  virtual Double_t expectedEvents(const RooArgSet&) const = 0;
  virtual TH1* getNominalHist() const = 0; //retrieves the underlying nominal histogram
  //-------
  
  virtual Double_t missingEvents() const; //default implementation just checks fMissingBin - override in stacks to combine components
  virtual RooAbsReal* createIntegralWM(const RooArgSet& iset,const char* rangeName = 0) const;
  
  void UseCurrentStyle();
  
  Int_t GetDimension() const { return fObservables.getSize(); }
  virtual TAxis* GetXaxis() const;
  virtual TAxis* GetYaxis() const;
  
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
  virtual Double_t getBinError(const RooFitResult& fr) const { return getBinVolume()*getError(fr); }
  inline virtual Double_t getBinContent(const RooArgSet* nset=0) const { return getBinVolume()*getVal(nset)*(expectedEvents(nset) /*+ missingEvents()*/); } //FIXME: assumes getVal is
  inline Double_t getBinContent(const RooArgSet& nset) const { return getBinContent(&nset); }   //  flat across the bin!
  RooRealVar* getStatFactor(int bin, bool createIf=false);
  
  Double_t GetBinVolume(int bin) const;
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
  
  //other methods to mimic TH1 behaviour 
  Double_t Integral(Option_t* opt="") const;
  
  
  //these functions are used when fitting
  RooAbsPdf&  model(); //returns self, with constraint terms added if necessary
  RooProdPdf* buildConstraints(const RooArgSet& obs, const char* systGroups="", bool addSelf=false) const;
  
  virtual void Draw(Option_t* option = "");
  virtual void Draw(Option_t* option,const TRooFitResult& r);
  virtual void Paint(Option_t* option = "");
  
  //other functions
  void setRooFitValV(bool in) { kUseAbsPdfValV = in; }//option to fall back to RooFit's usual evaluation
  virtual void setFloor(bool in, double floorValue=0.) { kMustBePositive = in; fFloorValue=floorValue;  } //if true, 'pdf' evaluations cannot go negative ... floorValue is what will be returned
  
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
  double fFloorValue=0.;
  
  TRooH1* fMissingBin = 0; //when FillMissing call, created a TH0D to hold this contribution to pdf normalization
  RooRealProxy fMissingBinProxy; //hold a proxy to it
  
  RooProdPdf* fThisWithConstraints = 0; //constructed in 'model' method.

private:
     const char* GetName() const { return dynamic_cast<const TNamed*>(this)->GetName(); }
  const char* GetTitle() const { return dynamic_cast<const TNamed*>(this)->GetTitle(); } 
    ClassDef(TRooAbsH1,1) // The Abstract Base class for all TRooFit pdfs
};

#endif
