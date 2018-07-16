#include "TRooFit/TRooWorkspace.h"

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1D.h"

#include "RooCategory.h"
#include "RooDataSet.h"

#include "TRooFit/Utils.h"

#define protected public
#include "RooGaussian.h"
#undef protected


TRooWorkspace::TRooWorkspace(const RooWorkspace& other) : RooWorkspace(other) {
  RooMsgService::instance().getStream(RooFit::INFO).removeTopic(RooFit::NumIntegration); //stop info message every time

  //dataset names do not get correctly imported by copy constructor!!
  auto otherData = other.allData();
  auto myData = allData();
  
  auto otherItr = otherData.begin();
  auto myItr = myData.begin();
  
  while( otherItr != otherData.end() ) {
    (*myItr)->SetName((*otherItr)->GetName());
    (*myItr)->SetTitle(strlen((*otherItr)->GetTitle())==0 ? (*otherItr)->GetName() : (*otherItr)->GetTitle());
    otherItr++;
    myItr++;
  }
  
  //now see if we are a histfactory workspace (i.e. we don't have TRooHStacks defined for all channels ...
  
  //create TRooHStack for each channel ...
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  while( (c = catIter->Next()) ) {
    TString cName(c->GetName());
    if(channel(cName)==0) {
      fIsHFWorkspace = true; //is a histfactory model!
    } else {
      continue; //no need to make a TRooHStack
    }
    RooRealSumPdf* channelPdf = dynamic_cast<RooRealSumPdf*>(pdf(cName+"_model"));
    TRooHStack* channelStack = 0;
    if(!channelPdf) {
      Error("setHF","Could not find model for channel %s",c->GetName());
      channelStack = new TRooHStack(c->GetName(),c->GetName());
    } else {
      channelStack = new TRooHStack(*channelPdf,*set("ModelConfig_Observables"));
      channelStack->SetName(c->GetName());channelStack->SetTitle(c->GetName());
      channelStack->SetMinimum(0);
    }
    import(*channelStack,RooFit::Silence());
    //create a dummy hist for each channel too
    fDummyHists[c->GetName()] = (channelStack->fDummyHist) ? (TH1*)channelStack->fDummyHist->Clone(c->GetName()) : new TH1D(c->GetName(),c->GetName(),1,0,1);
    fDummyHists[c->GetName()]->SetDirectory(0);
  }
  
  if(fIsHFWorkspace) {
    if(set("ModelConfig_Observables")) defineSet("obs",*set("ModelConfig_Observables"));
    if(set("ModelConfig_POI")) defineSet("poi",*set("ModelConfig_POI"));
    if(set("ModelConfig_GlobalObservables")) saveSnapshot("obsData_globalObservables",*set("ModelConfig_GlobalObservables"));
    
    //go through nuisance parameters and  guess errors for any params that dont have errors ..
    const RooArgSet* np = set("ModelConfig_NuisParams");
    if(np) {
      RooFIter itr = np->fwdIterator();
      while(RooAbsArg* arg = itr.next()) {
        RooRealVar* v = dynamic_cast<RooRealVar*>(arg);
        if(!v) continue;
        if(v->getError()) continue; //already has an error 
        if(v->isConstant()) continue; //shouldn't be the case 
        //HistFactory models have alpha parameters with range -5 to 5 ... error should be 1
        if(TString(v->GetName()).BeginsWith("alpha_") && fabs(v->getMin()+5.)<1e-9 && fabs(v->getMax()-5.)<1e-9) {
          //Info("TRooWorkspace","Guessing nominal error of %s is 1",v->GetName());
          v->setError(1);
        } else if(TString(v->GetName()).BeginsWith("gamma_stat_") && fabs(v->getMin())<1e-9) {
          //stat parameters have their max value set to 5*sigma ...
          v->setError( (v->getMax()-1.)/5. );
        }
      }
    }
  }
  
}


RooRealVar* TRooWorkspace::addParameter(const char* name, const char* title, double val, double min, double max, const char* constraintType) {
  factory(Form("%s[%f,%f]",name,min,max));
  RooRealVar* v = var(name);
  if(!v) return 0;
  v->setVal(val);
  v->SetTitle(title);
  
  if(constraintType) {
    v->setStringAttribute("constraintType",constraintType);
    TString cType(constraintType);
    cType.ToUpper();
    if(cType=="NORMAL") v->setError(1);  
  }
  
  return v;
}

RooRealVar* TRooWorkspace::addParameter(const char* name, const char* title, double min, double max, const char* constraintType) {
  return addParameter(name,title,(max+min)/2.,min,max,constraintType);
}

RooRealVar* TRooWorkspace::addObservable(const char* name, const char* title, double min, double max) {
  RooRealVar* out = addParameter(name,title,min,max);
  if(!out) return 0;
  
  if(!set("obs")) defineSet("obs",name);
  else extendSet("obs",name);
  
  return out;
}

bool TRooWorkspace::addArgument(const char* name, const char* title, double val) {
  factory(Form("%s[%f]",name,val));
  RooRealVar* v = var(name);
  if(!v) return false;
  v->SetTitle(title);
  v->setConstant(true);
  return true;
}

TRooHStack* TRooWorkspace::addChannel(const char* name, const char* title, const char* observable, int nBins, double min, double max) {
  RooCategory* c = cat("channelCat");
  if(!c) {
    factory(Form("channelCat[%s=0]",name));
    c = cat("channelCat");
    if(!set("obs")) defineSet("obs","channelCat");
    else extendSet("obs","channelCat");
    
    if(!c) return 0;
  }
  else {
    c->defineType(name);
  }
  
  
  
  RooRealVar* v = var(observable);
  if(!v) return 0;
  
  //create a TRooHStack and import it ..
  TRooHStack* hs = new TRooHStack(name,title);
  fStagedChannels.add(*hs);
  hs->setFloor(true,1e-9); //this also does the SetMinimim
  //import(*hs);
  
  //need to store a dummyHist for the binning ...
  fDummyHists[name] = new TH1D(observable,"Data",nBins,min,max);
  fDummyHists[name]->SetDirectory(0);
  
  return hs;
}

TRooHStack* TRooWorkspace::addChannel(const char* name, const char* title, const char* observable, int nBins, const double* bins) {
  RooCategory* c = cat("channelCat");
  if(!c) {
    factory(Form("channelCat[%s=0]",name));
    c = cat("channelCat");
    if(!set("obs")) defineSet("obs","channelCat");
    else extendSet("obs","channelCat");
    
    if(!c) return 0;
  }
  else {
    c->defineType(name);
  }
  
  
  
  RooRealVar* v = var(observable);
  if(!v) return 0;
  
  //create a TRooHStack and import it ..
  TRooHStack* hs = new TRooHStack(name,title);
  fStagedChannels.add(*hs);
  hs->SetMinimum(0);
  //import(*hs);
  
  //need to store a dummyHist for the binning ...
  fDummyHists[name] = new TH1D(observable,"Data",nBins,bins);
  fDummyHists[name]->SetDirectory(0);
  
  return hs;
}

TRooHStack* TRooWorkspace::addChannel(const char* name, const char* title, const char* observable) {
  RooRealVar* v = var(observable);
  if(!v) return 0;
  if(v->getBinning(name).isUniform()) {
    return addChannel(name,title,observable,v->getBins(name),v->getMin(name),v->getMax(name));
  } else {
    return addChannel(name,title,observable,v->getBins(name),v->getBinning(name).array());
  }
}

#include "TRegexp.h"

Int_t TRooWorkspace::setChannelAttribute(const char* channels,const char* attribute,Bool_t val) {
  TRegexp pattern(channels,true);
  TIterator* catIter = cat("channelCat")->typeIterator();
  int nCat(0);
  TObject* c;
  while( (c = catIter->Next()) ) {
    TString cName(c->GetName());
    if(cName.Contains(pattern)) {
      channel(c->GetName())->setAttribute(attribute,val);
      nCat++;
    } 
  }
  delete catIter;
  return nCat;
}

Int_t TRooWorkspace::setVarAttribute(const char* channels,const char* attribute,Bool_t val) {
  TRegexp pattern(channels,true);
  RooArgSet _allVars = allVars();
  RooFIter itr = _allVars.fwdIterator();
  RooAbsArg* c;
  int nCat(0);
  while( (c = itr.next()) ) {
    TString cName(c->GetName());
    if(cName.Contains(pattern)) {
      c->setAttribute(attribute,val);
      nCat++;
    } 
  }
  return nCat;
}

bool TRooWorkspace::addSample(const char* name, const char* title, const char* channels, bool allowNegative) {
  
  TRegexp pattern(channels,true);
  
  //loop over channels, create sample for each one that matches pattern
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  while( (c = catIter->Next()) ) {
    TString cName(c->GetName());
    if(cName.Contains(pattern)) {
      fDummyHists[c->GetName()]->Reset(); //ensures it is empty
      TString hName(Form("%s_%s",name,c->GetName()));
      TRooH1D* h = new TRooH1D(hName,title,*var(fDummyHists[c->GetName()]->GetName()),fDummyHists[c->GetName()]);
      h->SetFillColor( TRooFit::GetColorByName(name,true) );
      if(!allowNegative) h->setFloor(true,0); //by default, samples cannot go negative
      //import(*h);
      channel(c->GetName())->Add(h);
    }
  }
  
  delete catIter;
  
  return true;
  
}
  
bool TRooWorkspace::dataFill(const char* channelName, double x, double w) {
  if(!data("obsData")) {
    RooArgSet s(*set("obs"));
    RooRealVar w("weightVar","weightVar",1);
    s.add(w);
    RooDataSet* data = new RooDataSet("obsData","obsData",s,"weightVar");
    import(*data);
  }
  
  cat("channelCat")->setLabel(channelName);
  var( fDummyHists[channelName]->GetName() )->setVal( x );
  
  data("obsData")->add(*set("obs"),w);
  
  return true;
  
}

#include "TCut.h"

bool TRooWorkspace::sampleFill(const char* sampleName, TTree* tree, const char* weight, const char* variationName, double variationVal) {
  TString sWeight(weight);
  //loop over channels, create sample for each one that matches pattern
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  while( (c = catIter->Next()) ) {
    TRooH1* s = sample(sampleName,c->GetName());
    TRooHStack* cc = channel(c->GetName());
    
    TH1* histToFill = (TH1*)fDummyHists[c->GetName()]->Clone("tmpHist");
    histToFill->Reset();
    tree->Draw(Form("%s>>tmpHist",var( fDummyHists[c->GetName()]->GetName() )->getStringAttribute("formula")), TCut("cut",cc->getStringAttribute("formula"))*sWeight);
    if(variationName) {
      sampleAddVariation(sampleName,c->GetName(),variationName,variationVal,histToFill);
    } else {
      s->Add(histToFill);
    }
    delete histToFill;
    
  }
  delete catIter;
  
  return true;
}

Int_t TRooWorkspace::sampleFill(const char* sampleName, const char* channelName,  double x, double w) {
  return sample(sampleName,channelName)->Fill(x,w);
}
bool TRooWorkspace::sampleAdd(const char* sampleName, const char* channelName,  TH1* h1) {
  return sample(sampleName,channelName)->Add(h1);
}
bool TRooWorkspace::sampleAdd(const char* sampleName, const char* channelName,  RooAbsReal& arg) {
  return sample(sampleName,channelName)->Add(arg);
}

bool TRooWorkspace::sampleAddVariation(const char* sampleName, const char* channelName, const char* parName, double parVal, TH1* h1) {
  //get the parameter 
  if(!var(parName)) {
    Info("sampleAddVariation","%s is not defined. Creating a gaussian-constrained parameter",parName);
    addParameter(parName,parName,0,-5,5,"normal");
    return kFALSE;
  }
  return sample(sampleName,channelName)->AddVariation(*var(parName), parVal, h1);
}

void TRooWorkspace::SetVariationBinContent(const char* sampleName, const char* channelName, const char* parName, double parVal, Int_t bin, double val) {
  //get the parameter 
  if(!var(parName)) {
    Info("SetVariationBinContent","%s is not defined. Creating a gaussian-constrained parameter",parName);
    addParameter(parName,parName,0,-5,5,"normal");
  }
  sample(sampleName,channelName)->SetVariationBinContent(*var(parName), parVal, bin,val);
}

double TRooWorkspace::GetSampleCoefficient(const char* sampleFullName) const {
  RooAbsReal* samp = function(sampleFullName);
  if(!samp) {
    Error("GetSampleCoefficient","Unknown sample: %s",sampleFullName);
    return 1.;
  }
  //need to determine which channel this sample belongs to, because coefficient on sample may not be 1 (this is case in histfactory models) 
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  while( (c = catIter->Next()) ) {
    if(channel(c->GetName())->dependsOn(*samp)) break;
  }
  int idx = channel(c->GetName())->funcList().index(samp);
  if(idx!=-1) {
    return dynamic_cast<RooAbsReal*>(channel(c->GetName())->coefList().at(idx))->getVal();
  } 
  Error("GetSampleCoefficient","Failed to find sample coefficient!");
  return 1;
}

double TRooWorkspace::sampleIntegralAndError(double& err, const char* sampleFullName, const TRooFitResult& fr) const {
  RooAbsReal* samp = function(sampleFullName);
  if(!samp) return 0;
  if(samp->InheritsFrom("TRooAbsH1")) {
    return (dynamic_cast<TRooAbsH1*>(samp))->IntegralAndError(err,fr);;
  } else {
    //just a normal function, so we have to manually do integral and error ...
    RooArgSet* myObs = samp->getObservables(*const_cast<TRooWorkspace*>(this)->set("obs"));
   std::unique_ptr<RooAbsReal> inte( samp->createIntegral(*myObs) );
  
    RooArgSet* params = 0;
    RooArgSet* snap = 0;
  
    if(fr.floatParsFinal().getSize()||fr.constPars().getSize()) {
        params = samp->getParameters(*myObs);
        snap = (RooArgSet*)params->snapshot();
        *params = fr.floatParsFinal();
        *params = fr.constPars();
    }
  
    double out = inte->getVal();
    auto res = TRooAbsH1::getError(fr,&*inte,*myObs,0);
    
    err = res.first;
  
    if(fr.floatParsFinal().getSize()||fr.constPars().getSize()) {
      *params = *snap;
      delete snap;
    }
    delete myObs;
  
     double coef = GetSampleCoefficient(sampleFullName);
     err *= coef;
     out *= coef;
  
    return out;
  }
}

TRooH1* TRooWorkspace::sample(const char* sampleName, const char* channelName) {
  RooAbsReal* chan = dynamic_cast<RooAbsReal*>(channel(channelName));
  TRooH1* out = dynamic_cast<TRooH1*>(chan->findServer(Form("%s_%s",sampleName,channelName)));
  return out;
}

TRooHStack* TRooWorkspace::channel(const char* name) const {
  TRooHStack* out =  dynamic_cast<TRooHStack*>(pdf(name));
  if(!out) {
    out = dynamic_cast<TRooHStack*>(fStagedChannels.find(name));
  }
  return out;
}
  

//set fill color of sample in all channels
void TRooWorkspace::sampleSetFillColor(const char* sampleName, Int_t in) {
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  while( (c = catIter->Next()) ) {
    if(sample(sampleName,c->GetName())) sample(sampleName,c->GetName())->SetFillColor(in);
  }
}
void TRooWorkspace::sampleSetLineColor(const char* sampleName, Int_t in) {
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  while( (c = catIter->Next()) ) {
    if(sample(sampleName,c->GetName())) sample(sampleName,c->GetName())->SetLineColor(in);
  }
}

void TRooWorkspace::sampleScale(const char* sampleName,RooAbsReal& arg) {
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  while( (c = catIter->Next()) ) {
    if(sample(sampleName,c->GetName())) sample(sampleName,c->GetName())->Scale(arg);
  }
}

#include "TPRegexp.h"


RooSimultaneous* TRooWorkspace::model(const char* channels) {
  //builds the model for the given channels, putting them in a RooSimultaneous, then imports that and returns
  
  std::vector<TRegexp> patterns;
  TStringToken nameToken(channels,";");
  while(nameToken.NextToken()) {
      TString subName = (TString)nameToken;
      patterns.push_back(TRegexp(subName,true));
  }
  
  TRegexp pattern(channels,true);
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  
  TString simPdfName("simPdf");
  int nComps(0);
  TString factoryString("");
  
  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(pdf("simPdf"));
  
  while( (c = catIter->Next()) ) {
    bool pass=false;
    for(auto& pattern : patterns) if(TString(c->GetName()).Contains(pattern)) pass=true;
    if(pass==false) continue;
    if(channel(c->GetName())->getAttribute("isValidation")) continue; //don't include validation regions when constructing models
    nComps++;
    
    if(simPdf) {
      //take the channel models from the existing simPdf, to handle case where looking at histfactory model
      factoryString += Form(",%s=%s",c->GetName(),simPdf->getPdf(c->GetName())->GetName());
    } else {
      if( !pdf(Form("%s_with_Constraints",c->GetName())) ) {
        import( *channel(c->GetName())->buildConstraints(*set("obs"),"",true), RooFit::RecycleConflictNodes(), RooFit::Silence() );
      }
      factoryString += Form(",%s=%s_with_Constraints",c->GetName(),c->GetName());
      //remove from the stagedChannels list if its there (FIXME: is this a memory leak?)
      fStagedChannels.remove(*channel(c->GetName()),true/*silent*/,true/*match by name*/);
    }
    simPdfName += Form("_%s",c->GetName());
  }
  
  if(nComps>0) {
    if(nComps==cat("channelCat")->numTypes()) { simPdfName="simPdf"; } //all channels available
    if(pdf(simPdfName)) return static_cast<RooSimultaneous*>(pdf(simPdfName));
    factory(Form("SIMUL::%s(channelCat%s)",simPdfName.Data(),factoryString.Data()));
    
   
    
    
    //if got here then need to create the global observables snapshot too ...
    //need list of gobs then ..
    RooAbsPdf* m = pdf(simPdfName);
     //infer the global observables, nuisance parameters, model args (const np) 
    RooArgSet* gobs_and_np = m->getParameters(*set("obs"));

    //remove the poi ...
    //if(set("poi")) gobs_and_np->remove(*set("poi"));

    RooArgSet gobs;
    gobs.add(*gobs_and_np); //will remove the np in a moment

    //now pass this to the getAllConstraints method ... it will modify it to contain only the np!
    RooArgSet* s = m->getAllConstraints(*set("obs"),*gobs_and_np);
     delete s; //don't ever need this 

     //gobs_and_np now only contains the np
     gobs.remove(*gobs_and_np);
     //ensure all global observables are held constant now - important for pll evaluation
     gobs.setAttribAll("Constant",kTRUE);
     
     defineSet(Form("globalObservables%s",TString(simPdfName(6,simPdfName.Length()-6)).Data()),gobs);
     
     saveSnapshot(Form("%s_globalObservables%s",fCurrentData.Data(),TString(simPdfName(6,simPdfName.Length()-6)).Data()),gobs);
     
     //save the list of parameters as a set (poi+np)
     RooArgSet np;RooArgSet args;
      std::unique_ptr<TIterator> itr(gobs_and_np->createIterator());
      RooAbsArg* arg = 0;
      while((arg=(RooAbsArg*)itr->Next())) { 
          if(arg->isConstant()) args.add(*arg);
          else np.add(*arg);
      }
      defineSet(Form("params%s",TString(simPdfName(6,simPdfName.Length()-6)).Data()),np);
      defineSet(Form("args%s",TString(simPdfName(6,simPdfName.Length()-6)).Data()),args);
     
     delete gobs_and_np;
     
  }
  
  return static_cast<RooSimultaneous*>(pdf(simPdfName));
}

bool TRooWorkspace::generateAsimov(const char* name, const char* title, bool fitToObsData) {
  RooAbsPdf* thePdf = model();

  if(fitToObsData && data("obsData")) {
    //fit to data first, holding parameters of interest constant at their current values ...
    RooArgSet* snap = 0;
    if(set("poi")) {
      snap = (RooArgSet*)set("poi")->snapshot();
      const_cast<RooArgSet*>(set("poi"))->setAttribAll("Constant",kTRUE);
    }
    auto nll = TRooFit::createNLL(thePdf,data("obsData"),set("globalObservables"));
    loadSnapshot("obsData_globalObservables");
    auto fitResult = TRooFit::minimize(nll,true,false); //hesse not needed, just need best fit values 
    import(*fitResult,Form("conditionalFit_obsData_for%s",name));
    delete fitResult;
    if(snap) {
       *const_cast<RooArgSet*>(set("poi")) = *snap; //restores constant statuses 
       delete snap;
    }
  }
  
  auto asimov = TRooFit::generateAsimovDataset(thePdf,set("obs"),set("globalObservables"));
  
  saveSnapshot(Form("%s_globalObservables",name),*asimov.second,true);
  asimov.first->SetName(name);asimov.first->SetTitle(title);
  
  import(*asimov.first);
  
  delete asimov.first;
  delete asimov.second;
  
  return kTRUE;
  
}

std::pair<RooAbsData*,RooArgSet*> TRooWorkspace::generateToy(const char* name, const char* title, bool fitToObsData) {
  RooAbsPdf* thePdf = model();

  if(fitToObsData && data("obsData")) {
    //fit to data first, holding parameters of interest constant at their current values ...
    RooArgSet* snap = 0;
    if(set("poi")) {
      snap = (RooArgSet*)set("poi")->snapshot();
      const_cast<RooArgSet*>(set("poi"))->setAttribAll("Constant",kTRUE);
    }
    auto nll = TRooFit::createNLL(thePdf,data("obsData"),set("globalObservables"));
    loadSnapshot("obsData_globalObservables");
    auto fitResult = TRooFit::minimize(nll,true,false); //hesse not needed, just need best fit values 
    import(*fitResult,Form("conditionalFit_obsData_for%s",name));
    delete fitResult;
    if(snap) {
       *const_cast<RooArgSet*>(set("poi")) = *snap; //restores constant statuses 
       delete snap;
    }
  }
  
  auto out = TRooFit::generateToy(thePdf,set("obs"),set("globalObservables"));
  
  if(out.second) out.second->setName(Form("%s_globalObservables",name));
  out.first->SetName(name);out.first->SetTitle(title);
  
  return out;
  
}

double TRooWorkspace::pll(RooAbsData* theData, const RooArgSet* globalObservables, bool oneSided, bool discovery) {
  
  if(set("poi")==0) {
    Error("pll","You must define parameters of interest");
    return 0;
  }
  
  //RooArgSet* allSnap = (RooArgSet*)allVars().snapshot();
  RooArgSet _poi(*set("poi"));
  
  //do an unconditional fit ...
  RooArgSet* snap = (RooArgSet*)allVars().snapshot();
  _poi.setAttribAll("Constant",kFALSE);
  
  RooFitResult* unconditionalFit = fitTo(theData,globalObservables,false /*no need for hesse*/);
  
  allVars() = *snap; //restores values
  
  if(oneSided) {
    bool isGreater = ((RooAbsReal*)unconditionalFit->floatParsFinal().find(_poi.first()->GetName()))->getVal() > ((RooAbsReal*)_poi.first())->getVal();
    if( (!discovery && isGreater) || (discovery && !isGreater) ) {
          
      delete unconditionalFit;
      delete snap;
      return 0;
    }
          
  }
  
  
  //do a conditional fit ...
  _poi.setAttribAll("Constant",kTRUE);
  RooFitResult* conditionalFit = fitTo(theData,globalObservables,false /*no need for hesse*/);
  allVars() = *snap; //restores values
  
  double out = 2.0*(conditionalFit->minNll() - unconditionalFit->minNll());
  
  
  
  delete unconditionalFit;
  delete conditionalFit;
  delete snap;
  
  return out;
  
  
}


#include "TMultiGraph.h"

template <class T> struct greater_abs {
        bool operator() (const T& x, const T& y) const {return fabs(x)>fabs(y);}
      };


RooFitResult* TRooWorkspace::fitTo(RooAbsData* theData, const RooArgSet* globalObservables, bool doHesse) {
  RooSimultaneous* thePdf = model();
  TString simPdfName(thePdf->GetName());
  
  if(globalObservables) allVars() = *globalObservables; //sets the global observables;
  
  //check if we need to reduce the data (remove validation regions)
  std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
  TObject* c;
  TString excludedChannels;
  std::set<TString> excludedChannelsSet;
  while( (c = catIter->Next()) ) {
    if(!thePdf->getPdf(c->GetName())) {
      excludedChannelsSet.insert(c->GetName());
      if(excludedChannels!="") excludedChannels += " || ";
      excludedChannels += Form("channelCat == %d",((RooCatType*)c)->getVal());
    }
  }
  RooAbsData* reducedData=0;
  
  bool internalData = (data(theData->GetName())==theData); //data is already saved into workspace?
  
  if(excludedChannels!="" && theData->sumEntries(excludedChannels)>0) {
    TString dataName = TString::Format("reduced_%s_for_%s",theData->GetName(),thePdf->GetName());
    if(!internalData || !embeddedData(dataName)) {
      Info("fitTo","reducing %s to exclude validation channels",theData->GetName());
      TString cutString="!(";cutString+=excludedChannels;cutString+=")";
      
      //discovered that asimov datasets aren't properly reduced, they have been losing their weights ... so will do a manual reduction ..
      //if dataset is weighted
      if(theData->isWeighted()) {
        RooArgSet obs(*theData->get());RooRealVar weightVar("weightVar","weightVar",1);
        obs.add(weightVar);
        reducedData = new RooDataSet(theData->GetName(),theData->GetTitle(),obs,"weightVar");
        for(int i=0;i<theData->numEntries();i++) {
          if(excludedChannelsSet.find(theData->get(i)->getCatLabel("channelCat"))!=excludedChannelsSet.end()) continue; //skip excluded channels
          
          reducedData->add(*theData->get(),theData->weight());
        }
      } else {
        reducedData = theData->reduce(cutString);
      }
      
      reducedData->SetName(dataName);
      if(internalData) {
        import(*reducedData,RooFit::Embedded());
        delete reducedData;
      }
    }
    theData = (internalData) ? embeddedData(dataName) : reducedData;
  }
  
  //load or otherwise create the nll
  RooAbsReal* nll = (!internalData) ? 0 : function(Form("nll_%s_%s",theData->GetName(),thePdf->GetName()));
  if(!nll) {
    //it seems that if any parameters are const and we later make them unconst then the nll doesnt update properly (the const optimization is borked) 
    //so unconst all poi to be safe ...
    RooArgSet* snap = 0;
    if(set("poi")) {
      snap = (RooArgSet*)set("poi")->snapshot();
      const_cast<RooArgSet*>(set("poi"))->setAttribAll("Constant",kFALSE);
    }
    nll = TRooFit::createNLL(thePdf,theData,set(Form("globalObservables%s",TString(simPdfName(6,simPdfName.Length()-6)).Data())));
    if(snap) {
      *const_cast<RooArgSet*>(set("poi")) = *snap;
      delete snap;
    }
    nll->SetName(Form("nll_%s_%s",theData->GetName(),thePdf->GetName()));
    if(internalData) {
      import(*nll,RooFit::Silence());
      delete nll;
      nll = function(Form("nll_%s_%s",theData->GetName(),thePdf->GetName()));
    }
  }
  
  //any floating parameter that does not feature in the pdf or dataset should be held constant ... 
  //if this isn't done then saw getError method in TRooAbsH1 ends up putting variables like stat gamma's of validation regions into nset, which shouldn't be in the nset
  RooArgSet* allPars = thePdf->getParameters(RooArgSet());
  RooArgSet _allVars = allVars();
  RooArgSet* allFloatVars = static_cast<RooArgSet*>(_allVars.selectByAttrib("Constant",kFALSE));
  allFloatVars->remove(*allPars);allFloatVars->remove(*theData->get(),true,true);
  if(allFloatVars->getSize()) {
    Info("fitTo","Setting the following parameters constant: %s", allFloatVars->contentsString().c_str());
    allFloatVars->setAttribAll("Constant",kTRUE);
  }
  delete allPars; delete allFloatVars;
  
  if(!kDisabledForcedRecommendedOptions) TRooFit::setRecommendedDefaultOptions();
  RooFitResult* result = TRooFit::minimize(nll, true, doHesse);
  if(internalData) result->Print();
  result->SetName(Form("fitTo_%s",theData->GetName()));
  result->SetTitle(result->GetTitle());
  
  if(internalData) nll->setAttribute(result->GetName()); //used to track which nll corresponds to which fit
  
  if(reducedData && !internalData) delete reducedData;
  
  //before returning we will override _minLL with the actual NLL value ... offsetting could have messed up the value
  result->_minNLL = nll->getVal();
  
  if(!internalData) delete nll;
  
  return result;
  
}

RooAbsReal* TRooWorkspace::getFitNll(const char* fitName) {
  //This method looks for the nll function corresponding to the given fit name
  
  TString sFitName = (fitName) ? fitName : fCurrentFit.Data();
  
  if(sFitName=="") return 0;
  
  std::unique_ptr<RooAbsCollection> possibleNll( allFunctions().selectByAttrib(sFitName,true) );

  if(possibleNll->getSize()==0) return 0;
  else if(possibleNll->getSize()==1) return static_cast<RooAbsReal*>(possibleNll->first());
  else {
    Error("getFitNll","Found %d possible nll functions for %s ?? ",possibleNll->getSize(),sFitName.Data());
    return 0;
  }
}

void TRooWorkspace::impact(const char* impactPar) {
  //computes impact on given parameter, or parameter of interest
  
  TString parName;
  
  if(impactPar==0 && set("poi")->getSize()==1)  parName=set("poi")->first()->GetName();
  else if(impactPar) parName = impactPar;
  else {
    return;
  }
  
  RooAbsReal* nll = getFitNll();
  if(!nll) return;
  
  RooFitResult* out = getFit();
  if(!out) return;
  
  if(impactPar && out->floatParsInit().find(impactPar)) {
     RooArgSet* globalFloatPars = nll->getObservables(out->floatParsInit()); 
     
     TGraph posImpact; posImpact.SetFillColor(kCyan); //will hold impacts from positive errors
     TGraph negImpact; negImpact.SetFillColor(38); //hold impacts from negative errrors
     std::vector<TString> argNames;
     std::map<TString,std::pair<RooFitResult*,RooFitResult*>> impactResults;//save the fit restults in a map
    TMultiGraph impactGraph;
    
      
    
      //get impact for top X (most correlated) variables
      int myIdx = out->floatParsFinal().index(parName);
      std::vector<float> cors;
      for(int i=0;i<out->floatParsFinal().getSize();i++) {
        if(myIdx==i) continue;
        cors.push_back(out->correlation(myIdx,i));
      }
      std::sort(cors.begin(),cors.end(),greater_abs<float>());
    
      RooFIter parItr = out->floatParsInit().fwdIterator();
    Info("fitTo","%d variables to compute impact of ...",out->floatParsInit().getSize()-1);
    while( RooAbsArg* arg = parItr.next() ) {
      if(parName==arg->GetName()) continue; //dont compute self-impact!
      if(cors.size()>10 && fabs(out->correlation(myIdx,out->floatParsFinal().index(arg->GetName()))) < fabs(cors[10])) continue;
      //to compute e.g. post-fit impact due to lumi uncert, shift and fix par before redoing global fit ...
    std::cout << "correlation = " << out->correlation(myIdx,out->floatParsFinal().index(arg->GetName())) << std::endl;
      for(int j=-1;j<2;j+=2) { //j=0 is positive error, j=1 is negative errr
    
        *globalFloatPars = out->floatParsInit(); //restore initial state for global fit (includes restoring non-constant state)
        RooRealVar* parInNLL = ((RooRealVar*)nll->getParameters(RooArgSet())->find(arg->GetName()));
        RooRealVar* parInGlobalFit = ((RooRealVar*)out->floatParsFinal().find(arg->GetName()));
        parInNLL->setConstant(1);
        parInNLL->setVal( parInGlobalFit->getVal() + ((j<0) ? parInGlobalFit->getErrorLo() : parInGlobalFit->getErrorHi()) ); //set to desired value
        //rerun the fit ..
        Info("fitTo","Computing impact on %s due to %s",parName.Data(),arg->GetName());
        RooFitResult* impactFit = TRooFit::minimize(nll);
        double postfitImpact = ((RooRealVar*)impactFit->floatParsFinal().find(parName))->getVal() - ((RooRealVar*)out->floatParsFinal().find(parName))->getVal();

        if(j>0) {
          posImpact.SetPoint(posImpact.GetN(),posImpact.GetN()+1,postfitImpact);
          impactResults[arg->GetName()].second=impactFit;
        } else {
          negImpact.SetPoint(negImpact.GetN(),negImpact.GetN()+1,postfitImpact);
          impactResults[arg->GetName()].first=impactFit;
        }
      
      }
      argNames.push_back(arg->GetTitle());
   }
  
    impactGraph.Add(&negImpact);impactGraph.Add(&posImpact);
    impactGraph.SetName("impact");
    import(impactGraph);
    /*impactGraph.Draw("AB");
    for(int i=0;i<argNames_mu5.size();i++) {
     impactGraph.GetHistogram()->GetXaxis()->SetBinLabel(impactGraph.GetHistogram()->GetNbinsX()*(i+0.5)/(nPars-1),argNames[i]);
    }*/

  }
  
  
}

RooFitResult* TRooWorkspace::fitTo(const char* dataName, bool doHesse, const RooArgSet& minosPars) {
  RooAbsData* theData = data((dataName)?dataName:fCurrentData.Data());
  
  if(!theData) {
    Error("fitTo","Data %s is not in the workspace",dataName);
    return 0;
  }
  
  fCurrentData = theData->GetName(); //switch to this data
  
  model(); //call this method just the bring globalObservables snapshot into existence if necessary
  RooFitResult* result = fitTo(theData,getSnapshot(Form("%s_globalObservables",theData->GetName())),doHesse);
  import(*result,true/*overwrites existing*/); //FIXME: the statusHistory is not imported because RooFitResult Clone doesnt copy

  
  RooFitResult* out = loadFit(result->GetName());
  out->_statusHistory = result->_statusHistory; //HACK!!

  
  delete result;
  
  
  
  
  
  //delete nll;
  
  //if(reducedData) delete reducedData;
   
  return out;
}

RooFitResult* TRooWorkspace::loadFit(const char* fitName, bool prefit) {
  RooFitResult* result = dynamic_cast<RooFitResult*>(obj(fitName));
  if(!result) {
    Error("loadFit","Could not find fit %s",fitName);
    return 0;
  }
  
  allVars() = result->constPars();
  if(prefit) {
    allVars() = result->floatParsInit();
    fCurrentFitIsPrefit = true;
  } else {
    allVars() = result->floatParsFinal();
    fCurrentFitIsPrefit = false;
  }
  
  fCurrentFit = fitName;
  if(fCurrentFitResult) delete fCurrentFitResult;
  fCurrentFitResult = new TRooFitResult( result );
  
  return result;
  
  
}


#include "TCanvas.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TStyle.h"

TLegend* TRooWorkspace::GetLegend() {
  if(!fLegend) { 
    fLegend = new TLegend(0.5,1.-gStyle->GetPadTopMargin()-0.2,1.-gStyle->GetPadRightMargin(),1.-gStyle->GetPadTopMargin()); 
    fLegend->SetLineWidth(0);
    fLegend->SetFillStyle(0);
    fLegend->SetTextSize(gStyle->GetTextSize()*0.75);
  }
  return fLegend;
}

//draw a channel's stack and overlay the data too
void TRooWorkspace::channelDraw(const char* channelName, const char* opt, const TRooFitResult& res) {

  RooMsgService::instance().getStream(RooFit::INFO).removeTopic(RooFit::NumIntegration); //stop info message every time


  TVirtualPad* pad = gPad;
  if(!pad) {
    gROOT->MakeDefCanvas();
    pad = gPad;
  }

  TString sOpt(opt);
  sOpt.ToLower();
  
  if(fCurrentFitIsPrefit && !sOpt.Contains("init")) sOpt += "init";

  TLegend* myLegend = (TLegend*)GetLegend()->Clone("legend");
  myLegend->SetBit(kCanDelete); //so that it is deleted when pad cleared

  channel(channelName)->Draw(sOpt,res);
  
  pad->SetName(channel(channelName)->GetName());pad->SetTitle(channel(channelName)->GetTitle());
  
  if(channel(channelName)->getAttribute("isValidation")) gPad->SetFillColor(kGray);
  fDummyHists[channelName]->Reset();
  if(data(fCurrentData)) {
    //determine observables for this channel ...
    RooArgSet* chanObs = channel(channelName)->getObservables(data(fCurrentData));
    RooArgList l(*chanObs); delete chanObs;
    data(fCurrentData)->fillHistogram(fDummyHists[channelName], l, Form("channelCat==channelCat::%s",channelName));
    //update all errors to usual poisson
    for(int i=1;i<=fDummyHists[channelName]->GetNbinsX();i++) fDummyHists[channelName]->SetBinError(i,sqrt(fDummyHists[channelName]->GetBinContent(i)));
    fDummyHists[channelName]->SetMarkerStyle(20);
    fDummyHists[channelName]->SetLineColor(kBlack);
    fDummyHists[channelName]->Draw("same");
    myLegend->AddEntry(fDummyHists[channelName],data(fCurrentData)->GetTitle(),"lpe");
  }
  
  THStack* myStack = (THStack*)gPad->GetListOfPrimitives()->FindObject(Form("%s_stack",channel(channelName)->GetName()));
  if(myStack) {
    //add stack entries to legend
    TList* l = myStack->GetHists();
    for(int i=l->GetEntries()-1;i>=0;i--) { //go in reverse order
      myLegend->AddEntry(l->At(i),l->At(i)->GetTitle(),"f");
    }
  }
  
  //adjust legend size based on # of rows 
  
  myLegend->SetY1( myLegend->GetY2() - myLegend->GetNRows()*myLegend->GetTextSize() );
  myLegend->ConvertNDCtoPad();
  myLegend->Draw();

  
  if(!sOpt.Contains("same")) {
    TLatex latex;latex.SetNDC();latex.SetTextSize(gStyle->GetTextSize()*0.75);
    Double_t xPos = gPad->GetLeftMargin()+12./gPad->GetCanvas()->GetWw();
    Double_t yPos = 1. - gPad->GetTopMargin() -12./gPad->GetCanvas()->GetWh() - latex.GetTextSize();
    for(auto label : fLabels) {
      latex.DrawLatex(xPos,yPos,label);
      yPos -= latex.GetTextSize();
    }
    latex.DrawLatex(xPos,yPos,channel(channelName)->GetTitle());
    yPos -= latex.GetTextSize();
    if(fCurrentFit!="" && !fCurrentFitIsPrefit) {
      if(getFit(fCurrentFit)->status()%1000!=0) latex.SetTextColor(kRed);
      latex.DrawLatex(xPos,yPos,"Post-Fit");
      latex.SetTextColor(kBlack);
    }
    else latex.DrawLatex(xPos,yPos,"Pre-Fit");
  }

  pad->Update();

  //create space for the legend by adjusting the range on the stack ...
  double currMax = pad->GetUymax();
  double currMin = pad->GetUymin();
  
  if(currMax < fDummyHists[channelName]->GetMaximum()) currMax=fDummyHists[channelName]->GetMaximum();
  
  //double yScale = (currMax-currMin)/(1.-gPad->GetTopMargin()-gPad->GetBottomMargin());
  
  //want new max to give a newYscale where currMax is at shifted location 
  double newYscale = (currMax-currMin)/(myLegend->GetY1NDC()-pad->GetBottomMargin());
  
  currMax = currMin + newYscale*(1.-pad->GetTopMargin()-pad->GetBottomMargin());
  
  //currMax += yScale*(myLegend->GetY2NDC()-myLegend->GetY1NDC());
  
  channel(channelName)->GetYaxis()->SetRangeUser(currMin,currMax);
  pad->Modified(1);
  pad->Update();
  
}

Bool_t TRooWorkspace::writeToFile(const char* fileName, Bool_t recreate) {
  //ensure all staged channels are imported ...
  std::unique_ptr<TIterator> itr( fStagedChannels.createIterator() );
  TObject* obj;
  while( (obj = itr->Next()) ) {
    import(*static_cast<RooAbsPdf*>(obj));
  }
  fStagedChannels.removeAll();
  
  return RooWorkspace::writeToFile(fileName,recreate);
  
}


void TRooWorkspace::Draw(Option_t* option, const TRooFitResult& res) {
  //draws each channel
  TVirtualPad* pad = gPad;
  if(!pad) {
    gROOT->MakeDefCanvas();
    pad = gPad;
  }
  
  TString sOpt(option);
  sOpt.ToLower();
  
  if(sOpt.Contains("pull") || sOpt.Contains("cor")) {
    //just draw the pull plot of the fit result ..
    RooArgSet* visibleArgs = static_cast<RooArgSet*>(allVars().selectByAttrib("hidden",kFALSE));
    RooArgList* l = static_cast<RooArgList*>(fCurrentFitResult->floatParsFinal().selectCommon(*visibleArgs));
    delete visibleArgs;
    if(fCurrentFitResult) fCurrentFitResult->Draw(sOpt,*l);
    delete l;
    return;
  }
  

  
  
  if(!sOpt.Contains("same")) {
    pad->Clear();
  
    int nCat=0;
    TIterator* catIter = cat("channelCat")->typeIterator();
    TObject* c;
    while( (c = catIter->Next()) ) if(!channel(c->GetName())->getAttribute("hidden")) nCat++;
    delete catIter;
    
    if(nCat>1) {
      int nRows = nCat/sqrt(nCat);
      int nCols = nCat/nRows;
      if(nRows*nCols < nCat) nCols++;
      pad->Divide(nCols,nRows);
    }
  }
  
  
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  int i=1;
  while( (c = catIter->Next()) ) {
    if(channel(c->GetName())->getAttribute("hidden")) continue;
    pad->cd(i++);
    channelDraw(c->GetName(),option,res);
  }
  pad->cd(0);

  delete catIter;

}

//draws all channels, showing how values of channels depend on var
void TRooWorkspace::DrawDependence(const char* _var, Option_t* option) {
  RooRealVar* theVar = var(_var);
  if(!theVar) {
    Error("DrawDependence","%s not found",_var);
    return;
  }

  //draws each channel
  TVirtualPad* pad = gPad;
  if(!pad) {
    gROOT->MakeDefCanvas();
    pad = gPad;
  }
  
  //always draw to the whole canvas ...
  pad->GetCanvas()->cd();
  pad = gPad;
  
  TString sOpt(option);
  sOpt.ToLower();
  if(!sOpt.Contains("same")) {
    pad->Clear();
    
    int nCat=0;
    TIterator* catIter = cat("channelCat")->typeIterator();
    TObject* c;
    while( (c = catIter->Next()) ) if(!channel(c->GetName())->getAttribute("hidden")) nCat++;
    delete catIter;
    
    if(nCat>1) {
      int nRows = nCat/sqrt(nCat);
      int nCols = nCat/nRows;
      if(nRows*nCols < nCat) nCols++;
      pad->Divide(nCols,nRows);
    }
  }

  bool doSlice = sOpt.Contains("slice");
  sOpt.ReplaceAll("slice","");
  bool doPull = sOpt.Contains("pull");
  sOpt.ReplaceAll("pull","");
  
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  int i=1;
  while( (c = catIter->Next()) ) {
    TRooHStack* chan = channel(c->GetName());
    if(chan->getAttribute("hidden")) continue;
    pad->cd(i++);
    if(chan->getAttribute("isValidation")) gPad->SetFillColor(kGray);
    if(!chan->dependsOn(*theVar)) continue; //no dependence
    
    if(doSlice) {
      //perform 3 slices, around current value of parameter, based on range ...
      double curVal = theVar->getVal();
      
      //double maxVal = theVar->getMax();
      //double minVal = theVar->getMin();
      //double minDiff = std::min(maxVal-curVal,curVal-minVal);
      //double step = minDiff/2.;
      
      int cols[5] = {kOrange,kRed,kBlack,kBlue,kCyan};
      std::vector<TH1*> hists;
      for(int j=-1;j<2;j++) {
        //double val = curVal + j*step;
        double val = curVal; if(j<0) val += theVar->getErrorLo(); else if(j>0) val += theVar->getErrorHi();
        TRooFitResult r(TString::Format("%s=%f",_var,val));
        TH1* hist = (TH1*)chan->GetHistogram(&r,false)->Clone(chan->GetName());
        hist->SetDirectory(0);
        hist->SetTitle(Form("%s = %g",_var,val));
        hist->SetLineWidth(1);
        hist->SetLineColor(cols[j+2]);
        hist->SetBit(kCanDelete);
        hists.push_back(hist);
      }
      //hists[0]->Divide(hists[2]);hists[1]->Divide(hists[2]);hists[3]->Divide(hists[2]);hists[4]->Divide(hists[2]);
      //hists[2]->Divide(hists[2]);
      for(unsigned int j=0;j<hists.size();j++) {
        hists[j]->Draw((j==0)?"":"same");
      }
      
    } else if(doPull) {
     
     //draw the current and (if a fit is loaded, the prefit) pulls for the given variable ... defined here as (prediction-data)/sigma_data
     TH1* dataHist = (TH1*)fDummyHists[c->GetName()]->Clone("dataHist");
     RooArgSet* chanObs = chan->getObservables(data(fCurrentData));
     RooArgList l(*chanObs); delete chanObs;
     data(fCurrentData)->fillHistogram(dataHist, l, Form("channelCat==channelCat::%s",c->GetName()));
     
     double varVals[6] = {theVar->getVal(),theVar->getVal()+theVar->getErrorHi(),theVar->getVal()+theVar->getErrorLo(),0,0,0};
     RooRealVar* varInit = 0;
     if(getFit()) {
      varInit = dynamic_cast<RooRealVar*>(getFit()->floatParsInit().find(theVar->GetName()));
      if(varInit) {
        varVals[3] = varInit->getVal(); varVals[4] = varInit->getVal()+varInit->getErrorHi(); varVals[5] = varInit->getVal()+varInit->getErrorLo();
      }
     }
     int cols[6] = {kBlack,kBlue,kRed,kBlack,kBlue,kRed};
     for(int j=0;j<6;j++) {
      if(j>2 && !varInit) continue; //no prefit values available, so do not draw ..
      
      TRooFitResult r(TString::Format("%s=%f",_var,varVals[j]));
      TH1* myHist = (TH1*)chan->GetHistogram(&r,false)->Clone(chan->GetName());
      
      myHist->SetDirectory(0);
      myHist->SetTitle(Form("%s = %g",_var,varVals[j]));
      myHist->SetLineWidth(1);
      myHist->SetLineColor(cols[j]);
      myHist->SetLineStyle(1+(j>2)); //dashed lines for prefit
      myHist->SetBit(kCanDelete);
      myHist->Add( dataHist , -1. ); //subtract the data
      
      //scale each bin based on dataHist error
      for(int k=1;k<=myHist->GetNbinsX();k++) {
        double dataError = (myHist->GetBinContent(k)<0) ? (dataHist->GetBinContent(k)-0.5*TMath::ChisquareQuantile(TMath::Prob(1,1)/2.,2.*(dataHist->GetBinContent(k)))) : ((0.5*TMath::ChisquareQuantile(1.-TMath::Prob(1,1)/2.,2.*(dataHist->GetBinContent(k)+1))) - dataHist->GetBinContent(k));
        myHist->SetBinContent(k,myHist->GetBinContent(k)/dataError);
        myHist->SetBinError(k,0);
      }
      
      myHist->Draw((j==0)?"":"same");
      myHist->GetYaxis()->SetRangeUser(-2,2);
      gPad->Modified();
      
      
     }
     TLine ll;ll.SetLineStyle(2);ll.DrawLine(dataHist->GetXaxis()->GetXmin(),-1,dataHist->GetXaxis()->GetXmax(),-1);ll.DrawLine(dataHist->GetXaxis()->GetXmin(),1,dataHist->GetXaxis()->GetXmax(),1);
     delete dataHist;
     
     
    } else {
    
      TMultiGraph* allGraphs = new TMultiGraph; allGraphs->SetTitle(Form("%s;%s;%s",chan->GetTitle(),theVar->GetTitle(),( var(chan->GetObservableName(0))->numBins(chan->GetRangeName())==1 || !sOpt.Contains("bins"))?"Integral":"Bin Content - Current Content"));
    
      std::vector<TRooAbsH1*> comps;
      std::vector<TRooAbsH1*> myComps; //created here
      if(!sOpt.Contains("samples")) comps.push_back(chan);
      else {
        //break down by sample ...
        RooFIter fItr = chan->funcList().fwdIterator();
        RooAbsArg* arg;
        while( (arg = fItr.next()) ) {
          if(arg->InheritsFrom("TRooAbsH1")) {
            comps.push_back(dynamic_cast<TRooAbsH1*>(arg));
            break;
          }
          
          //if got here ... we will need to create a temporary TRooH1 for the sample and use that 
          TRooH1D* myComp = new TRooH1D(arg->GetName(),arg->GetTitle(),*var(chan->GetObservableName(0)),chan->GetRangeName());
          double coef = GetSampleCoefficient(arg->GetName());
          for(int i=1;i<=myComp->GetXaxis()->GetNbins();i++) {
            myComp->SetBinContent(i,coef);
          }
          myComp->Scale(*static_cast<RooAbsReal*>(arg));
          
          TString myTitle(myComp->GetTitle());
          if(myTitle.Contains(TString("_")+chan->GetTitle())) myTitle = myTitle.ReplaceAll(TString("_")+chan->GetTitle(),"");
          myTitle.ReplaceAll("L_x_",""); 
          myTitle.ReplaceAll("_overallSyst",""); 
          myTitle.ReplaceAll("_x_StatUncert","");
          myTitle.ReplaceAll("_x_HistSyst","");
          myTitle.ReplaceAll("_x_Exp","");
          myComp->SetTitle(myTitle); 
          myComp->SetFillColor(TRooFit::GetColorByName(myTitle,true));
          
          
          myComps.push_back(myComp);
          
          comps.push_back(myComp);
          
        }
      }
      
      for(auto& comp : comps) {
    
        TGraph2D* myGraph = new TGraph2D;
        
        RooRealVar* xVar = var(comp->GetObservableName(0));
        
        //graph needs at least 2 points in x-axis to render correctly delauny triangles
        comp->fillGraph(myGraph,RooArgList(*xVar,*theVar),(xVar->numBins()==1)?2:-1,21);
        myGraph->SetName(comp->GetName());
        myGraph->SetTitle(comp->GetTitle());
        /*myGraph->SetBit(kCanDelete);
        myGraph->Draw(option);*/
        
        //graph2D is now a series of points, split this up into 1D graphs ...
        int nBins = xVar->numBins(comp->GetRangeName());
        
        
        
        if(nBins==1 || myGraph->GetN() == nBins*21) {
          if(!sOpt.Contains("bins")) {
            //integrate the values ... across the bins 
            TGraph* gg = new TGraph;
            
            std::vector<double> pointVals(21,0.);
            for(int k=0;k<nBins;k++) {
              
              for(int j=k*21;j<(k+1)*21;j++) {
                pointVals[j%21] += myGraph->GetZ()[j]*comp->GetBinVolume(k+1) - (sOpt.Contains("samples"))*comp->GetBinContent(k+1); //shift to 0 if doing over samples
                gg->SetPoint(j%21,myGraph->GetY()[j],pointVals[j%21]); 
              }
            }
            allGraphs->Add(gg);
            
          } else {
            if(sOpt.Contains("samples") && !dynamic_cast<RooAbsReal*>(comp)->dependsOn(*theVar)) continue; //don't show non-dependent samples
            for(int k=0;k<nBins;k++) {
              TGraph* gg = new TGraph;
              if(sOpt.Contains("samples")) {
                gg->SetTitle(Form("%s bin %d",comp->GetTitle(),k+1));
              } else {
                gg->SetTitle(Form("Bin %d",k+1));
              }
              for(int j=k*21;j<(k+1)*21;j++) {
                double val = myGraph->GetZ()[j]*comp->GetBinVolume(k+1) - (nBins!=1)*comp->GetBinContent(k+1); //shift so that 0 = nominal if looking at multiple bins
                if(fabs(val)<1e-12) val = 0; //seems to be a slight discrepency between values ... possibly difference in getBinVolume (used in GetBinContent) vs GetBinVolume)
                gg->SetPoint(j%21,myGraph->GetY()[j],val); 
              }
              allGraphs->Add(gg);
            }
          }
        } else {
          Error("DrawDependence","Wrong number of points :-( ");
        }
        
        delete myGraph;
      } //loop over comps
      
      for(int j=0;j<allGraphs->GetListOfGraphs()->GetEntries();j++) {
        TGraph* gg = (TGraph*)allGraphs->GetListOfGraphs()->At(j);
        if(sOpt.Contains("bins")) {
          gg->SetLineColor(j+2);
        } else if(sOpt.Contains("samples")) {
          gg->SetTitle(comps[j]->GetTitle());
          gg->SetLineWidth(2);
          //gg->SetLineColor(j+2);
          if(!dynamic_cast<RooAbsReal*>(comps[j])->dependsOn(*theVar)) gg->SetLineStyle(2); //should we even show non-dependent samples?
          gg->SetLineColor(comps[j]->GetFillColor());
        } else {
          gg->SetLineWidth(2);
        }
      }
      
      
      allGraphs->SetBit(kCanDelete);
      allGraphs->Draw("AL");
      
      for(auto& comp : myComps) delete comp;
      
      
    }
    gPad->Update();
    
  }
  pad->cd(0);

  delete catIter;
}


void TRooWorkspace::Print(Option_t* opt) const {
  TString sOpt(opt);
  sOpt.ToLower();
  
  if(sOpt.Contains("yields")) RooMsgService::instance().getStream(RooFit::INFO).removeTopic(RooFit::NumIntegration); //stop info message every time
  
  if(sOpt.Contains("channels")) {
    //list the channels of this workspace ...
    std::unique_ptr<TIterator> catIter(cat("channelCat")->typeIterator());
    TObject* c;
    while( (c = catIter->Next()) ) {
      TRooHStack* chan = channel(c->GetName());
      std::cout << c->GetName();
      if(!chan) {
        std::cout << " - INVALID CHANNEL!!!" << std::endl;
        continue;
      }
      if(chan->getAttribute("hidden")) std::cout << " [HIDDEN]";
      if(chan->getAttribute("isValidation")) std::cout << " [BLINDED]";
      if(sOpt.Contains("yields")) {
        double err; 
        double inte = chan->IntegralAndError(err,(fCurrentFitResult)?*fCurrentFitResult:"");
        std::cout << Form(" [ %g +/- %g ]",inte,err);
      }
      std::cout << std::endl;
      if(sOpt.Contains("samples")) {
        RooFIter fItr = chan->funcList().fwdIterator();
        RooAbsArg* arg;
        while( (arg = fItr.next()) ) {
          std::cout << "   " << arg->GetName();
          if(sOpt.Contains("yields")) {
            double err; 
            double inte = sampleIntegralAndError(err,arg->GetName(),(fCurrentFitResult)?*fCurrentFitResult:"");
            std::cout << Form(" [ %g +/- %g ]",inte,err);
          }
          std::cout << std::endl;
        }
      }
    }
  } else if(sOpt.Contains("fits")) {
    //print RooFitResult's saved in the genericObjects list 
    auto genObjs = allGenericObjects();
    for(TObject* obj : genObjs) {
      if(!obj->InheritsFrom(RooFitResult::Class())) continue;
      std::cout << obj->GetName();
      if(fCurrentFit==obj->GetName()) std::cout << " <-- CURRENT FIT";
      std::cout << std::endl;
    }
  } else if(sOpt.Contains("params")) {
    std::map<TString,std::vector<TString>> paramStrings;
    //go through float parameters that are not observables (not in the 'obs' set)
    RooArgSet* _allVars = static_cast<RooArgSet*>(allVars().selectByAttrib("Constant",kFALSE));
    _allVars->remove(*(const_cast<TRooWorkspace*>(this)->set("obs")));
    RooArgSet _allPdfs = allPdfs();
    RooFIter itr = _allVars->fwdIterator();
    RooRealVar* _var = 0;
    while( (_var = dynamic_cast<RooRealVar*>(itr.next())) ) {
      //if(var->isConstant()) continue;if(set("obs")->find(*var)) continue;
      
      //check for constraint term .. which is a pdf featuring just this parameter 
      RooArgSet otherPars(*_allVars); otherPars.remove(*_var);
      RooFIter pitr = _allPdfs.fwdIterator();
      RooAbsPdf* constraintPdf = 0;
      while( RooAbsPdf* _pdf = dynamic_cast<RooAbsPdf*>(pitr.next()) ) {
        if(!_pdf->dependsOn(*_var)) continue;
        if(_pdf->dependsOn(otherPars)) continue;
        constraintPdf = _pdf; break;
      }
      if(constraintPdf) {
        if(constraintPdf->InheritsFrom("RooGaussian")) {
          //determine the "mean" (usually a global observable) and the standard deviation ...
          RooGaussian* gPdf = static_cast<RooGaussian*>(constraintPdf);
          double mean=0;
          if(strcmp(gPdf->x.arg().GetName(),_var->GetName())==0) mean = gPdf->mean;
          else mean = gPdf->x;
          double sigma= gPdf->sigma;
          paramStrings["gaussian"].push_back( Form("%s [%g,%g]",_var->GetName(),mean,sigma) );
        } else if(constraintPdf->InheritsFrom("RooPoisson")) {
          paramStrings["poisson"].push_back(_var->GetName());
        } else {
          paramStrings["other"].push_back(_var->GetName());
        }
      } else {
        paramStrings["unconstrained"].push_back(_var->GetName());
      }
      
    }
    delete _allVars;
    std::cout << std::endl;
    std::cout << "Unconstrained parameters (usually floating normalizations and parameters of interest)" << std::endl;
    std::cout << "------------------------" << std::endl;
    for(auto& s : paramStrings["unconstrained"]) {
      std::cout << s;
      if(const_cast<TRooWorkspace*>(this)->set("poi") && const_cast<TRooWorkspace*>(this)->set("poi")->find(s)) std::cout << " [PARAMETER OF INTEREST]";
      std::cout << std::endl;
    }
    if(paramStrings["gaussian"].size()) {
      std::cout << std::endl;
      std::cout << "Gaussian-constrained parameters (usually systematics)" << std::endl;
      std::cout << "------------------------" << std::endl;
      for(auto& s : paramStrings["gaussian"]) {
        std::cout << s << std::endl;
      }
    }
    if(paramStrings["poisson"].size()) {
      std::cout << std::endl;
      std::cout << "Poisson-constrained parameters (usually statistical uncertainties)" << std::endl;
      std::cout << "------------------------" << std::endl;
      for(auto& s : paramStrings["poisson"]) {
        std::cout << s << std::endl;
      }
    }
    if(paramStrings["other"].size()) {
      std::cout << std::endl;
      std::cout << "Unknown-constrained parameters" << std::endl;
      std::cout << "------------------------" << std::endl;
      for(auto& s : paramStrings["unconstrained"]) {
        std::cout << s << std::endl;
      }
    }
  } else if(sOpt.Contains("data")) {
    auto _allData = allData();
    for(auto d : _allData) {
      std::cout << d->GetName() << " (" << d->numEntries() << " entries)";
      if(d->isWeighted()) std::cout << " (" << d->sumEntries() << " weighted)";
      if(fCurrentData==d->GetName()) std::cout << " <-- CURRENT DATA";
      std::cout << std::endl;
      
    }
  
  } else {
    RooWorkspace::Print(opt);
  }
  
}



void TRooWorkspace::setDefaultStyle() {
  
  gStyle->SetOptStat(0);

  Int_t icol=0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetFrameFillColor(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetStatColor(icol);
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetTitleYOffset(1.0);
  Int_t font=42;//43; changed to relative because on multiplots with panels in vertical direction, the offset wasn't right for absolute font sizes
  Double_t tsize=0.045;//18;
   Double_t lsize = 0.04;//15;
  gStyle->SetTextFont(font);
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");
  gStyle->SetLabelSize(lsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(lsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(lsize,"z");
  gStyle->SetTitleSize(tsize,"z");
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.7);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetEndErrorSize(0.);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);

  gStyle->SetLegendFont(font);

   const Int_t NRGBs = 5;
   const Int_t NCont = 255;

   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}
