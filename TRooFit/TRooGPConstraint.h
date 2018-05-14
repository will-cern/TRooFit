#ifndef TROO_GP_CONSTRAINT
#define TROO_GP_CONSTRAINT

#include "RooAbsReal.h"
#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include "TRooH1.h"
#include "TMatrixD.h"

//pdf value is exp(-Sum_ij[ (y_i-y_ref_i)*Kinv_ij*(y_j-y_ref_j) ])
//where Kinv is the inverse of the kernel
//
//For this to be equivalent to a TRooChi2Constraint, you choose K (the kernel) to be
//a diagonal matrix with diagonal values equal to err^2 for each bin
//I.e. the kernel IS the diagonal covariance matrix!
//
//For gaussian process, the usual kernel would be given by
// K_ij = exp( -(x_i - x_j)^2 / l^2 )


class TRooGPConstraint : public RooAbsPdf {
public:

  // Constructors, assignment etc
  TRooGPConstraint(const char *name, const char *title, TRooH1& pdf, RooDataHist& ref, const TMatrixD& kernel, bool isInverted=false, RooRealVar* strength=0);
	    
  TRooGPConstraint(const TRooGPConstraint& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new TRooGPConstraint(*this,newname); }


  void setKernel(const TMatrixD& kernel, bool isInverted=false);
  
  virtual ~TRooGPConstraint();



protected:

  //RooRealProxy fPdf;
  //RooDataHist* fRef;
  //TMatrixD fKernelInv;
  
  RooRealProxy fGP;
  
  RooRealProxy fStrength; //a scale factor 
  
  virtual Bool_t selfNormalized() const { return true; }

  virtual Double_t evaluate() const { 
    double scaleFactor = (fStrength.absArg()) ? double(fStrength) : 1.;
    return std::exp(-(fGP)); 
  }
  
  ClassDef(TRooGPConstraint,1) // Generalised constraint, inspired by Gaussian Processes
};

#include "RooChi2Var.h"

class TRooGPVar : public RooChi2Var {
public:
  TRooGPVar(const char* name, const char* title, RooAbsPdf& pdf, RooDataHist& data, const TMatrixD& kernel, bool isInverted=false);
  virtual ~TRooGPVar();
  TRooGPVar(const TRooGPVar& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new TRooGPVar(*this,newname); }
  
  void setKernel(const TMatrixD& kernel, bool isInverted=false) { fKernelInv = kernel; if(!isInverted) fKernelInv.Invert(); setValueDirty(); }
  
protected:
  TMatrixD fKernelInv;
  
  virtual Double_t evaluatePartition(Int_t firstEvent, Int_t lastEvent, Int_t stepSize) const;
  
  ClassDef(TRooGPVar,1) //Variation of Chi2, where can specify a kernel too
  
};


#endif
