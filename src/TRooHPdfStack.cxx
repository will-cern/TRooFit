
#include "TRooFit/TRooHPdfStack.h"

#include "TMath.h"

//___________________________________
/* BEGIN_HTML
<p>A TRooHPdfStack is the TRooFit version of a RooAddPdf</p>
<p>You should think of the value of this pdf as being the sum of a sub-pdfs, where each pdf
is given a coefficient corresponding to the fraction of expectedEvents contributed by that pdf</p>

<p>i.e. pdfValue = Sum(i) coef_i*pdfValue_i, where coef_i is the fraction of total expected events of that pdf. 
The coef_i are automatically recomputed </p>

<p>Each subpdf is automatically individually normalized to 1, unlike in a TRooHStack of pdfs. </p>

<p>TRooHPdfStack can only combine actual pdfs, whereas a TRooHStack can combine any function (RooAbsReal).</p>

</p>
END_HTML */
//____________________________________

ClassImp(TRooHPdfStack) 

using namespace std;

TRooHPdfStack::TRooHPdfStack(const char* name, const char* title) : RooAddPdf(name,title),
  TRooAbsHStack(RooArgList(),this)
{
  //Constructor for a TRooHPdfStack, specify just the name and title 
  //
  //You then add components to the stack with TRooAbsHStack::Add
  _allExtendable = kTRUE;
}

TRooHPdfStack::TRooHPdfStack(const TRooHPdfStack& other, const char* name) :  
   RooAddPdf(other,name),
   TRooAbsHStack(other,this)
{ 
  //Copy Constructor
  //_allExtendable = kTRUE;
} 


TRooHPdfStack::TRooHPdfStack(const RooAddPdf& other, const RooArgSet& observables) :
  RooAddPdf(other,other.GetName()) {
  //Constructor for making a TRooHPdfStack object out of existing RooAddPdf 
  
  RooArgSet* obs = other.getObservables(observables);
  
  fObservables.add(*obs);  
  delete obs;
  
  if(fObservables.getSize()) fObservables.setName("obs");
  
  fNormFactors.setName("!normFactors");fShapeFactors.setName("!shapeFactors");fStatFactors.setName("!statFactors");
  
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


void TRooHPdfStack::reinit() {
  //internal method that is called when new components are added to the stack 
  //resets the integral managers etc etc
  
  //recreate the pdf iterator
  delete _pdfIter;
  _pdfIter = _pdfList.createIterator();
  
  if(_coefCache) delete[] _coefCache;
  _coefCache = new Double_t[_pdfList.getSize()];
  
  //reset the integration manager
  _normMgr.reset();
  //NOTE: maybe want to set _norm = 0 too ... depends if we see problems?
  
}

Double_t TRooHPdfStack::evaluate() const {
  if(fBlindRangeName.Length()) {
    RooFIter obsItr(fObservables.fwdIterator());
    while(auto obs = obsItr.next() ) {
      if(obs->inRange(fBlindRangeName)) return 0;
    }
  }
  double out = RooAddPdf::evaluate();
  
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
Double_t TRooHPdfStack::getValV(const RooArgSet* nset) const
{
  //Evaluates the TRooHPdfStack's raw value (if nset=0) or pdf value (if nset=observables)
  //The implementation is identical to RooAddPdf::getValV but with some minor 
  //changes so suppress warnings about negative values
  //
  //You can force the RooAddPdf implementation to be used by using TRooAbsH1::setRooFitValV(true)

  if(kUseAbsPdfValV) return RooAddPdf::getValV(nset);
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

    // Evaluate denominator
    Double_t normVal(_norm->getVal()) ;
    
    
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
    if(kMustBePositive && _value < fFloorValue) _value=fFloorValue;

    clearValueAndShapeDirty() ; //setValueDirty(kFALSE) ;
  } 

  return _value ;
}





