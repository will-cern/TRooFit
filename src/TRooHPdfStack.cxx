
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
  _allExtendable = kTRUE;
}

TRooHPdfStack::TRooHPdfStack(const TRooHPdfStack& other, const char* name) :  
   RooAddPdf(other,name),
   TRooAbsHStack(other,this)
{ 
  _allExtendable = kTRUE;
} 


void TRooHPdfStack::reinit() {
  //this method is called when new components are added to the stack 
  
  //recreate the pdf iterator
  delete _pdfIter;
  _pdfIter = _pdfList.createIterator();
  
  if(_coefCache) delete[] _coefCache;
  _coefCache = new Double_t[_pdfList.getSize()];
  
  //reset the integration manager
  _normMgr.reset();
  //NOTE: maybe want to set _norm = 0 too ... depends if we see problems?
  
}


//_____________________________________________________________________________
Double_t TRooHPdfStack::getValV(const RooArgSet* nset) const
{
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
    if(kMustBePositive && _value < 0) _value=0;

    clearValueAndShapeDirty() ; //setValueDirty(kFALSE) ;
  } 

  return _value ;
}





