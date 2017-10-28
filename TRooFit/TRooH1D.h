/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           * 
 *****************************************************************************/

#ifndef TROOH1D
#define TROOH1D

#include "TRooFit/TRooH1.h"

#include "RooRealVar.h"

//____________________
//
//                                                            //
// TRooH1D                                                    //
//                                                            //



class TRooH1D : public TRooH1 {
public:
  TRooH1D() : TRooH1() {
   //Default constructor ... generally do not use this
   }
  
  TRooH1D(const char *name, const char *title, RooAbsLValue& x, const char* binning=0) : 
    TRooH1(name,title,RooArgList(dynamic_cast<RooAbsArg&>(x)),{x.numBins(binning)},{(x.getBinningPtr(binning)) ? x.getBinningPtr(binning)->array() : 0}) {
     //Constructor taking either RooCategory or RooRealVar
     //In the latter case, you should specify the name of binning to use
     //which determines the number of bins in this TRooH1D (from x.numBins(binning))
     //
     //binnings can also be non-uniform, just add the non-uniform binning to the x variable
     //and then specify the name of the binning in this constructor
     
     
     
     }
  TRooH1D(const char *name, const char *title, RooRealVar& x, int nbins) : TRooH1(name,title,x,{nbins},{x.getMin()},{x.getMax()}) {
    //Constructor for use with continuous variables only, where you explicitly specify the number of bins
    //The min and max are taken from the RooRealVar 
   
   }
   
   TRooH1D(const char *name, const char *title, RooRealVar& x, int nbins, const Double_t* bins) : 
    TRooH1(name,title,x,{nbins},{bins}) {
      //Constructor for use with continuous variables only, where you explicitly specify the number of bins
      //You also specify the binning
     
     
     }
     
    
   TRooH1D(const char *name, const char *title, RooRealVar& x, TH1* hist) : TRooH1(name,title,x,0,0,0,{},hist) {
    //Constructor from an existing histogram
    //Will use as the nominal value
   
   }
   
  
  virtual void Draw(Option_t* option,const TRooFitResult& r);
  
  //Like Draw method above, but will use current parameter values and errors
  virtual void Draw(Option_t* option = "");
  
  
  
  TRooH1D(const TRooH1D& other, const char* name=0) : TRooH1(other,name) { 
    //need these for proper streaming 
  }
  virtual TObject* clone(const char* newname) const { return new TRooH1D(*this,newname); }
  
private:
  ClassDef(TRooH1D,1) // A TH1D that behaves like a RooFit Pdf
};

 
#endif
