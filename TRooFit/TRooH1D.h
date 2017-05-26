/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOH1D
#define TROOH1D

#include "TRooFit/TRooH1.h"

#include "RooRealVar.h"

class TRooH1D : public TRooH1 {
public:
  TRooH1D() : TRooH1() { }
  
  TRooH1D(const char *name, const char *title, RooAbsLValue& x, const char* binning=0) : 
    TRooH1(name,title,RooArgList(dynamic_cast<RooAbsArg&>(x)),{x.numBins(binning)},{(x.getBinningPtr(binning)) ? x.getBinningPtr(binning)->array() : 0}) { }
  TRooH1D(const char *name, const char *title, RooRealVar& x, int bins) : TRooH1(name,title,x,{bins},{x.getMin()},{x.getMax()}) { }
  
  virtual void Draw(const TRooFitResult& r, Option_t* option = "") { TRooAbsH1::Draw(r,option); }
  virtual void Draw(Option_t* option = "") { TRooAbsH1::Draw(option); }
  
  //need these for proper streaming
  TRooH1D(const TRooH1D& other, const char* name=0) : TRooH1(other,name) { }
  virtual TObject* clone(const char* newname) const { return new TRooH1D(*this,newname); }
  
private:
  ClassDef(TRooH1D,1)
};

 
#endif
