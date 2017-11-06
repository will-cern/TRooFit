#ifndef TROO_CHI2_VAR
#define TROO_CHI2_VAR

#include "RooAbsReal.h"
#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class TRooChi2Constraint : public RooAbsPdf {
public:

  // Constructors, assignment etc

  enum FuncMode { Function, Pdf, ExtendedPdf } ;

  TRooChi2Constraint(const char *name, const char *title, RooAbsPdf& pdf, RooDataHist& data,
	    Bool_t extended=kFALSE, const char* rangeName=0, const char* addCoefRangeName=0);
	    
  TRooChi2Constraint(const TRooChi2Constraint& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new TRooChi2Constraint(*this,newname); }


  
  virtual ~TRooChi2Constraint();



protected:

  RooRealProxy fChi2;
  
  virtual Bool_t selfNormalized() const { return true; }

  virtual Double_t evaluate() const { return std::exp(-0.5*fChi2); }
  
  ClassDef(TRooChi2Constraint,1) // Chi^2 function of p.d.f w.r.t a binned dataset. pdf value is exp(-0.5*chi^2)
};

// #include "RooGaussian.h"
// 
// class TRooGaussian : public RooGaussian {
// public:
//     using RooGaussian::RooGaussian;
//     
//   virtual Bool_t selfNormalized() const { return true; }
// 
// protected:
//   ClassDef(TRooGaussian,1)
// };


#endif
