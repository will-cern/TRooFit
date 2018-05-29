
#include "TRooFit/TRooGPConstraint.h"

#include "TVectorD.h"

ClassImp(TRooGPConstraint) 


using namespace std;



TRooGPConstraint::TRooGPConstraint(const char *name, const char *title, TRooH1& pdf, RooDataHist& ref, const TMatrixD& kernel, bool isInverted, RooRealVar* strengthVar)
	    : RooAbsPdf(name,title), /*fPdf("pdf","pdf",this,pdf), fRef(&ref), fKernelInv(kernel),*/ fGP("gp","gp",this,true,false,false/*owns arg*/), fStrength("gpStrength","strength",this)
{
 
  TRooGPVar* gp = new TRooGPVar("chi2","chi2",pdf,ref,kernel,isInverted);
  fGP.setArg(*gp); //will take ownership of it, so will manage its deletion
  
  if(strengthVar) fStrength.setArg(*strengthVar);
  
  
  //fKernelInv.Invert();
  //need to attach the buffers so that when we loop over fRef it updates pdf var too 
  //fRef->attachBuffers( *pdf.getDependents(*fRef->get()) );
  //FIXME: need to find a way to stop dependents from becoming variables in the fit
 
}
 



//_____________________________________________________________________________
TRooGPConstraint::TRooGPConstraint(const TRooGPConstraint& other, const char* name) : 
  RooAbsPdf(other,name),fGP("chi2",this,other.fGP),fStrength("strength",this,other.fStrength)
{
  // Copy constructor
  //fRef->resetBuffers();
  //fRef->attachBuffers( *fPdf.arg().getDependents(*fRef->get()) );
}





//_____________________________________________________________________________
TRooGPConstraint::~TRooGPConstraint()
{
  // Destructor
}

TMatrixD& TRooGPConstraint::getInvKernel(){ 
  return static_cast<TRooGPVar&>(const_cast<RooAbsReal&>(fGP.arg())).getInvKernel(); 
}

void TRooGPConstraint::setKernel(const TMatrixD& kernel, bool isInverted){ 
  static_cast<TRooGPVar&>(const_cast<RooAbsReal&>(fGP.arg())).setKernel(kernel,isInverted); 
  setValueDirty();
}

// Double_t TRooGPConstraint::evaluate() const {
// 
//   //need to loop over bins of RooDataHist, for each one get the value of the pdf, and compute 
//   //the penality term 
// 
//   double normFactor = ((const RooAbsPdf&)fPdf.arg()).expectedEvents(fRef->get());
//   TVectorD d(fRef->numEntries());
// 
//   for(int i=0;i<fRef->numEntries();i++) {
//     const RooArgSet* a = fRef->get(i);
//     if(!fRef->valid()) continue;
//     fPdf.arg().setValueDirty(); //seems like I have to do this??
//     d(i) = (fPdf.arg().getVal(a) * normFactor * fRef->binVolume()) - fRef->weight();
//   }
// 
// 
//   return std::exp(-(d*(fKernelInv*d)));
// 
// }


ClassImp(TRooGPVar)

TRooGPVar::TRooGPVar(const char *name, const char *title, RooAbsPdf& pdf, RooDataHist& hdata,
	    const TMatrixD& kernel, bool isInverted)
	    : RooChi2Var(name,title,pdf,hdata,true), fKernelInv(kernel)
{
  if(!isInverted) fKernelInv.Invert();
}


//_____________________________________________________________________________
TRooGPVar::TRooGPVar(const TRooGPVar& other, const char* name) : 
  RooChi2Var(other,name),fKernelInv(other.fKernelInv)
{
  // Copy constructor
}



//_____________________________________________________________________________
TRooGPVar::~TRooGPVar()
{
  // Destructor
}

#include "RooAbsDataStore.h"

//_____________________________________________________________________________
Double_t TRooGPVar::evaluatePartition(Int_t firstEvent, Int_t lastEvent, Int_t stepSize) const 
{
  if(stepSize!=1 || firstEvent!=0 || lastEvent!=_dataClone->numEntries()) {
    Warning("TRooGPVar","miisconfigured!!");
  }
  Int_t i ;
  _dataClone->store()->recalculateCache( _projDeps, firstEvent, lastEvent, stepSize, kFALSE) ;
  TVectorD d(lastEvent);
  // Determine normalization factor depending on type of input function
  Double_t normFactor(1) ;
  switch (_funcMode) {
  case Function: normFactor=1 ; break ;
  case Pdf: normFactor = _dataClone->sumEntries() ; break ;
  case ExtendedPdf: normFactor = ((RooAbsPdf*)_funcClone)->expectedEvents(_dataClone->get()) ; break ;
  }

  // Loop over bins of dataset
  RooDataHist* hdata = (RooDataHist*) _dataClone ;
  for (i=firstEvent ; i<lastEvent ; i+=stepSize) {
    // get the data values for this event
    hdata->get(i);
    if (!hdata->valid()) continue;
    const Double_t nData = hdata->weight() ;
    const Double_t nPdf = _funcClone->getVal(_normSet) * normFactor * hdata->binVolume() ;
    d(i) = nPdf-nData ;
  }
  
  //std::cout <<  d*(fKernelInv*d) << std::endl;
    
  return d*(fKernelInv*d) ;
}
