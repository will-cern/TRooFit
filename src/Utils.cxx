
#include "TRooFit/Utils.h"

ClassImp(TRooFit)

TRooFit::TRooFit() : TObject() { }

RooAbsPdf* TRooFit::BuildModel(TRooAbsH1& pdf, RooAbsData& data) {
  return pdf.buildConstraints( *data.get() , "", true );
}