/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * Extends a RooFitResult object, allowing for visualization
 * the same statFactors                                                      * 
 *****************************************************************************/

#ifndef TROOFITRESULT
#define TROOFITRESULT

#include "RooFitResult.h"
 
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h" 
#include "TGraphErrors.h"

#include "TLine.h"
#include "TBox.h"

#include "RooRealVar.h"
 
class TRooFitResult : public RooFitResult {
public:
  friend class TRooAbsH1; //needed to setConstPars in Draw method
  friend class TRooAbsHStack;
  TRooFitResult(const RooFitResult* res, double errorThreshold=-1)  : RooFitResult((res) ? *res : RooFitResult()) {
    if(!res) return;
    //construct from an existing RooFitResult
    //will discard floating parameters that had errors below the given threshold
    RooFIter itr2( floatParsFinal().fwdIterator() );
    RooAbsArg* arg = 0;
    RooArgList badParams;
    while( (arg = itr2.next()) ) {
      if(dynamic_cast<RooRealVar*>(arg)->getError() < errorThreshold) {
        badParams.add(*arg);
      }
    }
    if(badParams.getSize()) {
      //need to reduce ...
      RooArgList goodParams(*_finalPars); goodParams.remove(badParams,true,true);
      auto covMatrix = reducedCovarianceMatrix(goodParams);
      _finalPars->remove( badParams, true, true );
      _initPars->remove( badParams, true, true );
      setCovarianceMatrix( covMatrix );
    }
    
    adjustCovarianceMatrix();
    
   };
   
  TRooFitResult(const TRooFitResult& fr) : RooFitResult(fr) {
    //don't copy the draw objects though (anything that would get deleted in the destructor ... don't want a double delete)
    fPullFrame = 0;
    fPullGraph = 0;
    fCovHist = 0;
  }

  TRooFitResult(const char* name, const char* title, const RooArgList& finalPars); //will construct a covariance matrix assuming all uncorrelated

  TRooFitResult(const RooArgList& finalPars) : TRooFitResult(0,0,finalPars) { 
    //Construct for a list of parameter values. Will be used as the floatParsFinal
    //The covariance matrix is taken to be uncorrelated
  }; 
  
  TRooFitResult(const RooArgList& finalPars, const RooArgList& _constPars) : TRooFitResult(0,0,finalPars) {
    //constructor that also sets the constPars
    RooArgList tmp;
    tmp.addClone(constPars()); //some of the finalPars may have been const, ending up in here
    tmp.addClone(_constPars);
    setConstParList(tmp);
  }
  
  TRooFitResult(const char* constPars=""); 
  
  virtual void Paint(Option_t* opt) {
    TString sOpt(opt);
    if(sOpt.Contains("cov")||sOpt.Contains("cor")) {
      if(fCovHist) fCovHist->Paint("sameaxis");
    }
    else {
      if(fPullFrame) {
        //fPullFrame->Paint(option);
        for(auto& o : fPullBoxes) o.Paint("3");
        for(auto& o : fPullLines) o.Paint("l");
        fPullFrame->Paint("sameaxis");
      }
      if(fPullGraph) fPullGraph->Paint("p");
    }
    
  }
  
  std::map<TString, TRooFitResult*> GetAssociatedFitResults() const { return fAssociatedFitResults; }
  
  void AddAssociatedFitResult(const char* name, TRooFitResult* result) { fAssociatedFitResults[name] = result; }
  
  
  virtual void Draw(Option_t* option = "pull") { Draw(option,RooArgList()); }
  virtual void Draw(Option_t* option,const RooArgList& args);
  virtual void Draw(Option_t* option,const char* argFilter);
  
  TGraphAsymmErrors* GetPullGraph() { return fPullGraph; }
  
  ~TRooFitResult() {
    if(fPullGraph) delete fPullGraph;
    if(fPullFrame) delete fPullFrame;
  }
  
protected:
  void adjustCovarianceMatrix();
  void resetCovarianceMatrix();

private:
  void init(const RooArgList& pars);


  TH1D* fPullFrame = 0;
  TGraphAsymmErrors* fPullGraph = 0;
  std::vector<TGraph> fPullLines;
  std::vector<TGraphErrors> fPullBoxes;

  TH1* fCovHist = 0; //covariance histogram

  std::map<TString, TRooFitResult*> fAssociatedFitResults; //can associate other fits to this fit

  ClassDef(TRooFitResult,1) // An extended/improved version of a RooFitResult
};
 
#endif
