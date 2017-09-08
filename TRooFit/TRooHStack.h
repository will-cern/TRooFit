/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * A RooRealSumPdf that behaves like a THStack                               * 
 * TRooHPdfStack is different in that the components are individually
 * all pdfs, and are individually all normalized  (it is like a RooAddPdf)
 *****************************************************************************/

#ifndef TROOHSTACK
#define TROOHSTACK

#include "RooRealSumPdf.h"
#include "TRooFit/TRooAbsHStack.h"



class TRooHStack : public RooRealSumPdf, public TRooAbsHStack {
public:
  friend class TRooH1;
  using TRooAbsH1::Draw;
  using RooRealSumPdf::cloneTree;
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooRealSumPdf::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooRealSumPdf::getVal(nset); }
  TIterator* clientIterator() const { return RooRealSumPdf::clientIterator(); }
  virtual RooAbsArg* cloneTree(const char* newname=0) const { return RooRealSumPdf::cloneTree(newname); }
  RooArgSet* getDependents(const RooArgSet& set) const { return RooRealSumPdf::getDependents(set); }
  RooArgSet* getParams(const RooArgSet& set) const { return RooRealSumPdf::getParameters(set); }
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const { Double_t out = RooRealSumPdf::getNorm(nset); if(out<0 && kMustBePositive) return 0; return out; }
  virtual Double_t expectedEvents(const RooArgSet& nset) const { Double_t out = RooRealSumPdf::getNorm(nset); if(out<0 && kMustBePositive) return 0; return out; }
  ///----
  virtual Double_t missingEvents() const { return TRooAbsHStack::missingEvents(); }


  TRooHStack() {} ; 
  
  TRooHStack(const char* name, const char* title);
  virtual TObject* clone(const char* newname) const { return new TRooHStack(*this,newname); }
  TRooHStack(const TRooHStack& other, const char* name=0) ;
  
  

  //virtual const char* GetRangeName() const { if(fRooHists.getSize()==0) return 0; return fRooHists[0].GetName(); }
  virtual void Paint(Option_t* option = "") {  TRooAbsHStack::Paint(option); }
  
  virtual void Draw(const TRooFitResult& r, Option_t* option = "") { TRooAbsHStack::Draw(r,option); }
  virtual void Draw(Option_t* option = "") { TRooAbsHStack::Draw(option); }
  
  //using TRooAbsHStack::Draw;
  
 
  //override getValV so we can suppress warnings about 0 and negative values
  virtual Double_t getValV( const RooArgSet* set = 0 ) const;

  virtual void setFloor(Bool_t in) { RooRealSumPdf::setFloor(in); kMustBePositive=in; }

  virtual void printMetaArgs(std::ostream& os) const;

protected:
  ///Methods required by TRooAbsHStack
  virtual TIterator*& compIter() { return _funcIter; }
  virtual RooListProxy& compList() { return _funcList; }
  virtual const RooListProxy& compList() const { return _funcList; }
  virtual void reinit();

  
private:

  ClassDef(TRooHStack,1) // A RooRealSumPdf that behaves like a THStack
};

 
#endif
