#include "TRooFit/TRooWorkspace.h"

#include "TRooFit/TRooHStack.h"
#include "TRooFit/TRooH1D.h"

#include "RooCategory.h"
#include "RooDataSet.h"


TRooWorkspace::TRooWorkspace(const RooWorkspace& other) : RooWorkspace(other) {
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
  
  
}


RooRealVar* TRooWorkspace::addParameter(const char* name, const char* title, double val, double min, double max, const char* constraintType) {
  factory(Form("%s[%f,%f]",name,min,max));
  RooRealVar* v = var(name);
  if(!v) return 0;
  v->setVal(val);
  v->SetTitle(title);
  
  if(constraintType) v->setStringAttribute("constraintType",constraintType);
  
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
  hs->SetMinimum(0);
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

bool TRooWorkspace::addSample(const char* name, const char* title, const char* channels) {
  
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

bool TRooWorkspace::sampleFill(const char* sampleName, TTree* tree, const char* weight) {
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
    s->Add(histToFill);
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
    Error("sampleAddVariation","%s is not defined. Please call addParameter method first to define the parameter",parName);
    return kFALSE;
  }
  return sample(sampleName,channelName)->AddVariation(*var(parName), parVal, h1);
}

TRooH1* TRooWorkspace::sample(const char* sampleName, const char* channelName) {
  RooAbsReal* chan = dynamic_cast<RooAbsReal*>(channel(channelName));
  TRooH1* out = dynamic_cast<TRooH1*>(chan->findServer(Form("%s_%s",sampleName,channelName)));
  return out;
}

TRooHStack* TRooWorkspace::channel(const char* name) {
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
#include "TRooFit/Utils.h"

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
  TString factoryString("");
  int nComps(0);
  
  while( (c = catIter->Next()) ) {
    bool pass=false;
    for(auto& pattern : patterns) if(TString(c->GetName()).Contains(pattern)) pass=true;
    if(pass==false) continue;
    nComps++;
    if(channel(c->GetName())->getAttribute("isValidation")) continue; //don't include validation regions when constructing models
    if( !pdf(Form("%s_with_Constraints",c->GetName())) ) {
      import( *channel(c->GetName())->buildConstraints(*set("obs"),"",true), RooFit::RecycleConflictNodes(), RooFit::Silence() );
    }
    
    factoryString += Form(",%s=%s_with_Constraints",c->GetName(),c->GetName());
    simPdfName += Form("_%s",c->GetName());
    
    //remove from the stagedChannels list if its there (FIXME: is this a memory leak?)
    fStagedChannels.remove(*channel(c->GetName()),true/*silent*/,true/*match by name*/);
    
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
     
     saveSnapshot(Form("obsData_globalObservables%s",TString(simPdfName(6,simPdfName.Length()-6)).Data()),gobs);
     
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

RooFitResult* TRooWorkspace::fitTo(const char* dataName) {
  RooAbsData* theData = data((dataName)?dataName:fCurrentData.Data());
  
  if(!theData) {
    Error("fitTo","Data %s is not in the workspace",dataName);
    return 0;
  }
  
  fCurrentData = dataName; //switch to this data
  
  
  RooAbsPdf* thePdf = model();
  auto nll = TRooFit::createNLL(thePdf,theData,set("globalObservables"));
  
  //load the globalObservables that go with the data ...
  loadSnapshot(Form("%s_globalObservables",theData->GetName()));
  
  
  RooFitResult* result = TRooFit::minimize(nll);
  import(*result,Form("fitTo_%s",theData->GetName()));
  delete result;
  delete nll;
  
  result = loadFit(Form("fitTo_%s",theData->GetName()));
  
  if(result->status()%1000!=0) {
    Error("fitTo","Fit status is %d",result->status());
    loadFit(Form("fitTo_%s",theData->GetName()),true);
  }
   
  return result;
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

  TVirtualPad* pad = gPad;
  if(!pad) {
    gROOT->MakeDefCanvas();
    pad = gPad;
  }

  

  TLegend* myLegend = (TLegend*)GetLegend()->Clone("legend");
  myLegend->SetBit(kCanDelete); //so that it is deleted when pad cleared

  channel(channelName)->Draw(opt,res);
  
  
  
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

  TString sOpt(opt);
  sOpt.ToLower();
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
    if(fCurrentFit!="" && !fCurrentFitIsPrefit) latex.DrawLatex(xPos,yPos,"Post-Fit");
    else latex.DrawLatex(xPos,yPos,"Pre-Fit");
  }

  gPad->Update();

  //create space for the legend by adjusting the range on the stack ...
  double currMax = gPad->GetUymax();
  double currMin = gPad->GetUymin();
  
  //double yScale = (currMax-currMin)/(1.-gPad->GetTopMargin()-gPad->GetBottomMargin());
  
  //want new max to give a newYscale where currMax is at shifted location 
  double newYscale = (currMax-currMin)/(myLegend->GetY1NDC()-gPad->GetBottomMargin());
  
  currMax = currMin + newYscale*(1.-gPad->GetTopMargin()-gPad->GetBottomMargin());
  
  //currMax += yScale*(myLegend->GetY2NDC()-myLegend->GetY1NDC());
  
  channel(channelName)->GetYaxis()->SetRangeUser(currMin,currMax);
  gPad->Modified(1);
  gPad->Update();
  
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
  if(!sOpt.Contains("same")) {
    pad->Clear();
  
    int nCat = cat("channelCat")->numTypes();
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
  
    int nCat = cat("channelCat")->numTypes();
    if(nCat>1) {
      int nRows = nCat/sqrt(nCat);
      int nCols = nCat/nRows;
      if(nRows*nCols < nCat) nCols++;
      pad->Divide(nCols,nRows);
    }
  }

  bool doSlice = sOpt.Contains("slice");
  sOpt.ReplaceAll("slice","");
  
  TIterator* catIter = cat("channelCat")->typeIterator();
  TObject* c;
  int i=1;
  while( (c = catIter->Next()) ) {
    pad->cd(i++);
    TRooHStack* chan = channel(c->GetName());
    if(!chan->dependsOn(*theVar)) continue; //no dependence
    
    if(doSlice) {
      //perform 5 slices, around current value of parameter, based on range ...
      double curVal = theVar->getVal();
      double maxVal = theVar->getMax();
      double minVal = theVar->getMin();
      
      double minDiff = std::min(maxVal-curVal,curVal-minVal);
      double step = minDiff/2.;
      
      int cols[5] = {kOrange,kRed,kBlack,kBlue,kCyan};
      std::vector<TH1*> hists;
      for(int j=-2;j<3;j++) {
        double val = curVal + j*step;
        TRooFitResult r(TString::Format("%s=%f",_var,val));
        TH1* hist = (TH1*)chan->GetHistogram(&r,false)->Clone(chan->GetName());
        hist->SetTitle(Form("%s = %g",_var,val));
        hist->SetLineWidth(1);
        hist->SetLineColor(cols[j+2]);
        hist->SetBit(kCanDelete);
        hists.push_back(hist);
      }
      hists[0]->Divide(hists[2]);hists[1]->Divide(hists[2]);hists[3]->Divide(hists[2]);hists[4]->Divide(hists[2]);
      hists[2]->Divide(hists[2]);
      for(int j=0;j<5;j++) {
        hists[j]->Draw((j==0)?"":"same");
      }
      
    } else {
    
      TGraph2D* myGraph = new TGraph2D;
      //graph needs at least 2 points in x-axis to render correctly delauny triangles
      chan->fillGraph(myGraph,RooArgList(*var(chan->GetObservableName(0)),*theVar),(var(chan->GetObservableName(0))->numBins()==1)?2:-1,10);
      myGraph->SetName(chan->GetName());
      myGraph->SetTitle(chan->GetTitle());
      myGraph->SetBit(kCanDelete);
      myGraph->Draw(option);
    }
    gPad->Update();
    
  }
  pad->cd(0);

  delete catIter;
}



/*
std::pair<RooDataSet*,RooArgSet*> RooHypoModel::createToyDataSet(bool doBinned,const char* channels) {

   RooAbsPdf* thePdf = model(channels);

   //determine global observables for this model ...
   RooArgSet* gobs_and_np = model->getParameters(*set("obs"));

   //remove the poi ...
   gobs_and_np->remove(*set("poi"));

   RooArgSet gobs;
   gobs.add(*gobs_and_np); //will remove the np in a moment


   //now pass this to the getAllConstraints method ... it will modify it to contain only the np!
   RooArgSet* s = model->getAllConstraints(*set("obs"),*gobs_and_np);
   delete s; //don't ever need this 

   //gobs_and_np now only contains the np
   gobs.remove(*gobs_and_np);

   RooArgSet* toy_gobs = new RooArgSet;gobs.snapshot(*toy_gobs);
   //and then generate global and non-global observables 

   RooDataSet* globals = thePdf->generateSimGlobal(gobs,1);
   *toy_gobs = *globals->get(0);
   delete globals;

   RooDataSet* ddata;
   if(!thePdf->canBeExtended()) {
      //must be a counting experiment ... generate 1 dataset entry per category ... so must recurse through simultaneous pdf
      //usually only one simpdf though ... because user can use a RooSuperCategory to combine several categories!
      ddata = new RooDataSet("toyData","toyData",*set("obs"));
      buildDataset(thePdf,*ddata, false);
   } else {
      //just generate as normal 
      //have to first check if expected events are zero ... if it is, then the dataset is just empty ('generate' method fails when expectedEvents = 0)
      if(thePdf->expectedEvents(*set("obs"))==0) {
         ddata = new RooDataSet("toyData","toyData",*set("obs"));
      } else {
         //std::cout << " expected = " << thePdf->expectedEvents(*data()->dataset()->get()) << std::endl;
         ddata = (doBinned) ? thePdf->generate(*set("obs"),RooFit::Extended(),RooFit::AllBinned()) : thePdf->generate(*set("obs"),RooFit::Extended());
      }
   }

  std::pair<RooDataSet*,RooArgSet*> out = std::make_pair(ddata,toy_gobs);

   return out;

}


void RooHypoModel::buildDataset(RooAbsPdf* thePdf, RooDataSet& out, RooArgSet& gobs, bool doAsimov ) {
   //check if this pdf is a simultaneous 
   if(thePdf->IsA()==RooSimultaneous::Class()) {
      //need to get the observable set for each category 
      RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(thePdf);
      //loop over index category 
      std::unique_ptr<TIterator> itr(simPdf->indexCat().typeIterator());
      RooCatType* tt = NULL;
      while((tt=(RooCatType*) itr->Next())) {
         std::cout << " in category " << tt->GetName() << " of " << thePdf->GetName() << std::endl;
         const_cast<RooArgSet*>(out.get())->setCatIndex(simPdf->indexCat().GetName(),tt->getVal());
         if(simPdf->getPdf(tt->GetName())) buildDataset(simPdf->getPdf(tt->GetName()), out, doAsimov);
      }
      return;
   }

   //got here, so not in a simultaneous 


   if(doAsimov) {
      //... so get the terms of the pdf that constrain the observables:
      RooArgSet params(poi()); params.add(np()); //the poi and np
      RooArgSet* constraints = thePdf->getAllConstraints(gobs,params);
      //constraints are now all the parts of the pdf that do not feature global observables ... so they must feature the observables!
      //loop over constraint pdfs, recognise gaussian, poisson, lognormal
      TIterator* citer = constraints->createIterator();
      RooAbsPdf* pdf = 0;
      while((pdf=(RooAbsPdf*)citer->Next())) { 
         //determine which obs this pdf constrains. There should only be one!
         std::unique_ptr<RooArgSet> cgobs(pdf->getObservables(*out.get()));
         if(cgobs->getSize()!=1) { std::cout << "constraint " << pdf->GetName() << " constrains " << cgobs->getSize() << " observables. skipping..." << std::endl; continue; }
         RooAbsArg* gobs_arg = cgobs->first();
   
         //now iterate over the servers ... the first server to depend on a nuisance parameter or the poi we assume is the correct server to evaluate ..
         std::unique_ptr<TIterator> itr(pdf->serverIterator());
         for(RooAbsArg* arg = static_cast<RooAbsArg*>(itr->Next()); arg != 0; arg = static_cast<RooAbsArg*>(itr->Next())) {
            RooAbsReal * rar = dynamic_cast<RooAbsReal *>(arg); 
            if( rar && ( rar->dependsOn(np()) || rar->dependsOn(poi()) ) ) {
//               std::cout << " setting " << gobs_arg->GetName() << " equal to value of " << arg->GetName() << " (" << rar->getVal() << ")" << std::endl;
               const_cast<RooArgSet*>(out.get())->setRealValue(gobs_arg->GetName(),rar->getVal());
            }
         }
      }
      delete citer;
      //now add the current values to the dataset 
      out.add(*out.get());
   } else {
      //just generate a single toy with the list of observables for all the things are actually depend on
      RooArgSet* obs = thePdf->getObservables(*out.get());
      RooDataSet* toy = thePdf->generate(*obs,1);
      const_cast<RooArgSet*>(out.get())->assignValueOnly(*toy->get(0));
      out.add(*out.get());
      delete toy;delete obs;
   }

}

*/



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