/*****************************************************************************
 * Project: TRooFit   - an extension for RooFit                              *
 *                                                                           *
 * Modified version of a 
 * RooWorkspace ... Helps with model building                          * 
 *****************************************************************************/

#ifndef TROOWORKSPACE
#define TROOWORKSPACE

#include "RooWorkspace.h"

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1.h"

#include "TTree.h"

#include "RooSimultaneous.h"
 
class TRooWorkspace : public RooWorkspace {
public:
  using RooWorkspace::RooWorkspace; 
  
  bool definePoi(const char* poi) { return defineSet("poi",poi); }
  bool addArgument(const char* name, const char* title, double val);
  bool addParameter(const char* name, const char* title, double min, double max, const char* constraintType=0);
  bool addObservable(const char* name, const char* title, double min, double max);
  bool addChannel(const char* name, const char* title, const char* observable, int nBins, double min, double max);
  bool addChannel(const char* name, const char* title, const char* observable, int nBins, const double* bins);
  bool addChannel(const char* name, const char* title, const char* observable);
  bool addSample(const char* name, const char* title, const char* channels="*");
  
  bool dataFill(const char* channel, double x, double w=1.);
  Int_t sampleFill(const char* sample, const char* channel, double x, double w=1.);
  bool sampleAdd(const char* sample, const char* channel,  TH1* h1);
  bool sampleAdd(const char* sample, const char* channel, RooAbsReal& arg);
  
  bool sampleFill(const char* sample, TTree* tree, const char* weight); //fills given sample in all channels where formula have been defined
  
  //add a normalization factor to a sample, across all channels
  void sampleScale(const char* sample,RooAbsReal& arg);
  
  //set fill color of sample in all channels
  void sampleSetFillColor(const char* sample, Int_t in);
  void sampleSetLineColor(const char* sample, Int_t in);
  
  TRooH1* sample(const char* sampleName, const char* channelName);
  TRooHStack* channel(const char* name);
  
  
  RooSimultaneous* model(const char* channels="*");
  
  //draw a channel's stack and overlay the data too
  void channelDraw(const char* channel, Option_t* option="e3005", const TRooFitResult& res = "");
  //draws all channels
  
  virtual void Draw(Option_t* option, const TRooFitResult& res = "");
  virtual void Draw(Option_t* option="e3005") { Draw(option,""); }
  
  Bool_t writeToFile(const char* fileName, Bool_t recreate=kTRUE);
  
private:

  std::map<TString,TH1*> fDummyHists;
  
  RooArgList fStagedChannels; //channels cannot be added until they are frozen (all content filled)
  

  ClassDef(TRooWorkspace,1) // An extended form of a RooWorkspace
};
 
#endif
