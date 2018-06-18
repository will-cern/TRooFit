
#include "TRooFit/TRooHStack.h"

#include "TMath.h"

#include "RooConstVar.h"

ClassImp(TRooHStack) 

using namespace std;

TRooHStack::TRooHStack(const char* name, const char* title) : RooRealSumPdf(name,title),
  TRooAbsHStack(RooArgList(),this)
{
  //Constructor for a TRooHStack, specify just the name and title 
  //
  //You then add components to the stack with TRooAbsHStack::Add

  _extended = true; _haveLastCoef = true;
}

TRooHStack::TRooHStack(const TRooHStack& other, const char* name) :  
   RooRealSumPdf(other,name),
   TRooAbsHStack(other,this)
{ 
  //Copy constructor

  _extended = true; _haveLastCoef = true;
}
 
TRooHStack::TRooHStack(const RooRealSumPdf& other, const RooArgSet& observables) :
  RooRealSumPdf(other,other.GetName()) {
  //Constructor for making a TRooHStack object out of existing RooRealSumPdf 
  
  RooArgSet* obs = other.getObservables(observables);
  
  fObservables.add(*obs);  
  delete obs;
  
  if(fObservables.getSize()) fObservables.setName("obs");
  
  std::vector<int> bins; std::vector<const Double_t*> binEdges;
    //take binning from the observables .. 
    for(int i=0;i<fObservables.getSize();i++) {
      if(fObservables[i].IsA() == RooRealVar::Class()) {
        RooAbsBinning& myBins = static_cast<RooRealVar&>(fObservables[i]).getBinning();
        bins.push_back( myBins.numBins() );
        binEdges.push_back( myBins.array() );
      }
    }
  
  
  //create the default hist
    if(fObservables.getSize()==0) {
      fDummyHist = new TH1D(other.GetName(),other.GetTitle(),1,-0.5,0.5);
      //specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator");
    } else if(fObservables.getSize()==1) {
      if(dynamic_cast<TObject&>(fObservables[0]).InheritsFrom(RooAbsCategory::Class())) { 
        RooAbsCategory& cat = static_cast<RooAbsCategory&>(fObservables[0]);
        fDummyHist = new TH1D(other.GetName(),other.GetTitle(),cat.numTypes(),-0.5,cat.numTypes()-0.5);
      } else {
        fDummyHist = new TH1D(other.GetName(),other.GetTitle(),bins[0],binEdges[0]);
      }
      //specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator");
    } else if(fObservables.getSize()==2) {
        fDummyHist = new TH2D(other.GetName(),other.GetTitle(),bins[0],binEdges[0],bins[1],binEdges[1]);
    } else {
      std::cout << "not supported!" << std::endl;
    }
 
 
  
  fDummyHist->SetDirectory(0);fDummyHist->Sumw2();
  fDummyHist->GetXaxis()->SetTitle("");fDummyHist->GetYaxis()->SetTitle("");
  
 
  
  
  
}


void TRooHStack::reinit() {
  //internal method that is called when new components are added to the stack 
  //resets the integral managers etc etc 
  
  
  //recreate the func iterator
  delete _funcIter;
  _funcIter = _funcList.createIterator();
  
  //coefList needs extending 
  _coefList.removeAll();
  for(int i=0;i<_funcList.getSize();i++) {
    _coefList.addOwned( *(new RooConstVar(Form("%s_coef%d",GetName(),i),Form("%s_coef%d",GetName(),i),1)) ); 
  }
  //recreate iterator
  delete _coefIter;
  _coefIter = _coefList.createIterator();
  
  //reset the integration manager (what's the diff to below??)
  _normIntMgr.reset();
  
  //reset the RooAbsPdf integration manager
  _normMgr.reset();
  //NOTE: maybe want to set _norm = 0 too ... depends if we see problems?
  
}

void TRooHStack::printMetaArgs(ostream& os) const {
  //internal method that is called when Print method is called 

  //don't print coefficients since we made all of them const = 1
 // Customized printing of arguments of a RooRealSumPdf to more intuitively reflect the contents of the
  // product operator construction

  _funcIter->Reset() ;
  _coefIter->Reset() ;

  Bool_t first(kTRUE) ;
    
  RooAbsArg* func ;
    while((func=(RooAbsArg*)_funcIter->Next())) {
      if (!first) {
	os << " + " ;
      } else {
	first = kFALSE ;
      }
      os << func->GetName() ; 
    }  

  os << " " ;    
}

Double_t TRooHStack::evaluate() const {
  if(fBlindRangeName.Length()) {
    RooFIter obsItr(fObservables.fwdIterator());
    while(auto obs = obsItr.next() ) {
      if(obs->inRange(fBlindRangeName)) return 0;
    }
  }
  double out = RooRealSumPdf::evaluate();
  
  //multiply by all the norm factors
  if(fNormFactors.getSize()) {
    RooFIter itr(fNormFactors.fwdIterator());
    while( RooAbsReal* arg = (RooAbsReal*)itr.next() ) out *= arg->getVal(); //NOTE: should we use _normSet? leads to issues if normfactor is a pdf...
  }
  //and by the shape factors for this bin
  if(fBinsShapeFactors.size()) {
    int bin = getBin(GetName()); //forcefully use OUR binning 
    auto&& fItr = fBinsShapeFactors.find(bin);
    if(fItr!=fBinsShapeFactors.end()) {
      for(auto& sfIdx : fItr->second) {
        out *= ((RooAbsReal&)fShapeFactors[sfIdx]).getVal();
      }
    }
  }
  
  return out;
  
}

//_____________________________________________________________________________
Double_t TRooHStack::getValV(const RooArgSet* nset) const
{
  //Evaluates the TRooHStack's raw value (if nset=0) or pdf value (if nset=observables)
  //The implementation is identical to RooRealSumPdf::getValV but with some minor 
  //changes so suppress warnings about negative values
  //
  //You can force the RooRealSumPdf implementation to be used by using TRooAbsH1::setRooFitValV(true)

  if(kUseAbsPdfValV) return RooRealSumPdf::getValV(nset);
  ///THIS CODE IS COPIED FROM BUT SLIGHTLY MODIFIED VERSION OF RooAbsPdf::getValV
  ///DONE TO SUPPRESS WARNINGS ABOUT NEGATIVE VALUES

  // Return current value, normalizated by integrating over
  // the observables in 'nset'. If 'nset' is 0, the unnormalized value. 
  // is returned. All elements of 'nset' must be lvalues
  //
  // Unnormalized values are not cached
  // Doing so would be complicated as _norm->getVal() could
  // spoil the cache and interfere with returning the cached
  // return value. Since unnormalized calls are typically
  // done in integration calls, there is no performance hit.

  // Fast-track processing of clean-cache objects
  //   if (_operMode==AClean) {
  //     cout << "RooAbsPdf::getValV(" << this << "," << GetName() << ") CLEAN  value = " << _value << endl ;
  //     return _value ;
  //   }

  // Special handling of case without normalization set (used in numeric integration of pdfs)
  if (!nset) {
    RooArgSet* tmp = _normSet ;
    _normSet = 0 ;
    
    Double_t val = evaluate() ;
    _normSet = tmp ;
    Bool_t error = (TMath::IsNaN(val)) ? traceEvalPdf(val) : false; //allows negative and zero values

    //if(kMustBePositive && val<0) { val=0; } ... raw values can be negative

    if (error) {
//       raiseEvalError() ;
      return 0 ;
    }
    return val ;
  }


  // Process change in last data set used
  Bool_t nsetChanged(kFALSE) ;
  if (nset!=_normSet || _norm==0) {
    nsetChanged = syncNormalization(nset) ;
  }

  // Return value of object. Calculated if dirty, otherwise cached value is returned.
  if (isValueDirty() || nsetChanged || _norm->isValueDirty()) {

    // Evaluate numerator
    Double_t rawVal = evaluate() ;
    Bool_t error = (TMath::IsNaN(rawVal)) ? traceEvalPdf(rawVal) : false; //allows negative and zero values // Error checking and printing

    if((kMustBePositive||getFloor()||getFloorGlobal()) && rawVal < fFloorValue) rawVal=fFloorValue;

    // Evaluate denominator
    Double_t normVal(_norm->getVal()) ;
    
    //June2018: should we add this line to the other classes too?
    if((kMustBePositive||getFloor()||getFloorGlobal()) && normVal < fFloorValue) normVal=fFloorValue;
    
//     if (normVal<=0.) {
//       error=kTRUE ;
//       logEvalError("p.d.f normalization integral is zero or negative") ;  
//     }

    // Raise global error flag if problems occur
    if (error) {
//       raiseEvalError() ;
      _value = 0 ;
    } else {
      _value = rawVal / normVal ;
//       cout << "RooAbsPdf::getValV(" << GetName() << ") writing _value = " << rawVal << "/" << normVal << " = " << _value << endl ;
    }
    if(rawVal==0 && normVal==0) _value=1;
    

    clearValueAndShapeDirty() ; //setValueDirty(kFALSE) ;
  } 

  return _value ;
}





