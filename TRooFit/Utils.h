//static methods that are useful for various activities
//they live in a class so that the root dictionary is generated for them
//so you can use them on the command line

#ifndef TROOFIT
#define TROOFIT

#include "RooAbsPdf.h"
#include "TRooH1.h"


class TRooFit : public TObject  {

  public:
    TRooFit();
    
    //takes a pdf and dataset and builds a model from it (i.e. the pdf + constraint terms)
    static RooAbsPdf* BuildModel(TRooAbsH1& pdf, RooAbsData& data);

  
    ClassDef(TRooFit,1);

};

#endif