/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * Modified version of a 
 * TSlider ... controls changing a RooAbsRealLValue                          * 
 *****************************************************************************/

#ifndef TROOSLIDER
#define TROOSLIDER

#include "TSlider.h"
#include "RooAbsRealLValue.h"

#include "TCanvas.h"
 
class TRooSlider : public TSlider {
public:
  TRooSlider() : TSlider() { }
  
  TRooSlider(RooAbsRealLValue* var) :TSlider(var->GetName(),var->GetTitle(),-0.5,-1,0.5,-0.5), fVar(var) { 
    SetObject(this); //calls ourself when needing an update 
  }
  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) {
    double val = GetMaximum()*(fVar->getMax() - fVar->getMin()) + fVar->getMin(); //based on max
    if(fabs(val - fVar->getVal())<1e-9) return;
    fVar->setVal(val);
    gPad->GetCanvas()->Modified();
    gPad->GetCanvas()->Update();
  }
  
private:

  RooAbsRealLValue* fVar = 0;

  ClassDef(TRooSlider,1) // Your description goes here...
};
 
#endif
