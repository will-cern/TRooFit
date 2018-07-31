
#include "TRooFit/TRooExtendedBinding.h"

ClassImp(TRooExtendedBinding)

TRooExtendedBinding::TRooExtendedBinding(const char* name, const char* title, RooAbsPdf& _pdf, const RooArgSet& _nset) :
  RooAbsReal(name,title),
  pdf("pdf","pdf",this,_pdf),
  nset("nset","nset",this)
{
  nset.add(_nset);
}

TRooExtendedBinding::TRooExtendedBinding(const TRooExtendedBinding& other, const char* name) :
  RooAbsReal(other,name),
  pdf("pdf",this,other.pdf),
  nset("nset",this,other.nset)
{

}

Double_t TRooExtendedBinding::evaluate() const {
  return ((RooAbsPdf&)pdf.arg()).expectedEvents(nset);
}