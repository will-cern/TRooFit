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
#include "TGraphAsymmErrors.h" 

#include "TLine.h"
#include "TBox.h"
 
class TRooFitResult : public RooFitResult {
public:
  friend class TRooAbsH1; //needed to setConstPars in Draw method
  friend class TRooAbsHStack;
  TRooFitResult(const RooFitResult* res)  : RooFitResult(*res) {
    //construct from an existing RooFitResult
   };

  TRooFitResult(const char* name, const char* title, const RooArgList& finalPars); //will construct a covariance matrix assuming all uncorrelated

  TRooFitResult(const RooArgList& finalPars) : TRooFitResult(0,0,finalPars) { 
    //Construct for a list of parameter values. Will be used as the floatParsFinal
    //The covariance matrix is taken to be uncorrelated
  }; 
  
  TRooFitResult(const char* constPars=""); 
  
  virtual void Paint(Option_t* option = "") {
    if(fPullFrame) {
      fPullFrame->Paint(option);
      for(auto& o : fPullBoxes) o.Paint("f");
      for(auto& o : fPullLines) o.Paint("l");
      fPullFrame->Paint("sameaxis");
    }
    if(fPullGraph) fPullGraph->Paint("p");
  }
  
  virtual void Draw(Option_t* option = "pull");
  
  TGraphAsymmErrors* GetPullGraph() { return fPullGraph; }
  
protected:
  

private:
  void init(const RooArgList& pars);


  TH1D* fPullFrame = 0;
  TGraphAsymmErrors* fPullGraph = 0;
  std::vector<TLine> fPullLines;
  std::vector<TBox> fPullBoxes;

  ClassDef(TRooFitResult,1) // An extended/improved version of a RooFitResult
};
 
#endif
