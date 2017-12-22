/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * A RooAddPdf that behaves like a THStack                                   * 
 *****************************************************************************/

#ifndef TROOHPDFSTACK
#define TROOHPDFSTACK

#include "RooAddPdf.h"
#include "TRooFit/TRooAbsHStack.h"



class TRooHPdfStack : public RooAddPdf, public TRooAbsHStack {
public:
  friend class TRooH1;
  using TRooAbsH1::Draw;
  using RooAddPdf::cloneTree;
  ///Methods required by TRooAbsH1
  const char* GetName() const { return TNamed::GetName(); }
  const char* GetTitle() const { return TNamed::GetTitle(); }
  Double_t getVal(const RooArgSet* nset = 0) const { return RooAddPdf::getVal(nset); }
  Double_t getVal(const RooArgSet& nset) const { return RooAddPdf::getVal(nset); }
  virtual Double_t expectedEvents(const RooArgSet* nset=0) const { Double_t out = RooAddPdf::expectedEvents(nset); if(out<fFloorValue && kMustBePositive) return fFloorValue; return out; }
  virtual Double_t expectedEvents(const RooArgSet& nset) const { Double_t out = RooAddPdf::expectedEvents(nset); if(out<fFloorValue && kMustBePositive) return fFloorValue; return out; }
  ///----
  virtual Double_t missingEvents() const { return TRooAbsHStack::missingEvents(); }

  void resetNormMgr() { _normMgr.reset(); }
  
  TRooHPdfStack() {} ; 
  
  TRooHPdfStack(const char* name, const char* title);
  virtual TObject* clone(const char* newname) const { return new TRooHPdfStack(*this,newname); }
  TRooHPdfStack(const TRooHPdfStack& other, const char* name=0) ;
  
  

  //virtual const char* GetRangeName() const { if(fRooHists.getSize()==0) return 0; return fRooHists[0].GetName(); }
  virtual void Paint(Option_t* option = "") {  TRooAbsHStack::Paint(option); }
  
  virtual void Draw(Option_t* option,const TRooFitResult& r) { TRooAbsHStack::Draw(option,r); }
  virtual void Draw(Option_t* option = "") { TRooAbsHStack::Draw(option); }
 
  void PrintCoefs() const { for(int i=0;i<_pdfList.getSize();i++) std::cout << _coefCache[i] << std::endl; }


  //override getValV so we can suppress warnings about 0 and negative values
  virtual Double_t getValV( const RooArgSet* set = 0 ) const;

  

protected:
  ///Methods required by TRooAbsHStack
  virtual TIterator*& compIter() { return _pdfIter; }
  virtual RooListProxy& compList() { return _pdfList; }
  virtual const RooListProxy& compList() const { return _pdfList; }
  virtual void reinit();
  virtual double evaluate() const;
  
private:

  ClassDef(TRooHPdfStack,1) // An RooAddPdf that behaves like a THStack
};

 
#endif
