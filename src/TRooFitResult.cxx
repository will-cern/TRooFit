
#include "TRooFit/TRooFitResult.h"
#include "RooRealVar.h"

#include "TPRegexp.h"

#include "TROOT.h"
#include "TPad.h"

ClassImp(TRooFitResult) 

using namespace std;


void TRooFitResult::init(const RooArgList& pars) {
  //check all are RooRealVar, remove any that are not
  RooFIter itr( pars.fwdIterator() );
  RooAbsArg* arg = 0;
  RooArgList realpars;
  RooArgList constpars;
  while( (arg = itr.next()) ) {
    if(arg->isConstant()) constpars.add(*arg);
    else if(arg->IsA() != RooRealVar::Class()) {
      Warning("TRooFitResult","%s is not a RooRealVar, ignoring",arg->GetName());
    } 
    else realpars.add(*arg);
  }
  
  setConstParList(constpars);
  setInitParList(realpars);
  setFinalParList(realpars);
  
  TMatrixDSym cov(realpars.getSize());
  
  RooFIter itr2( realpars.fwdIterator() );
  
  int i=0;
  while( (arg = itr2.next()) ) {
    cov(i,i) = pow(dynamic_cast<RooRealVar*>(arg)->getError(),2);
    i++;
  }
  setCovarianceMatrix(cov);
}


TRooFitResult::TRooFitResult(const char* name, const char* title, const RooArgList& pars) 
  : RooFitResult(name,title) {
  
  init(pars);
  
}

TRooFitResult::TRooFitResult(const char* constPars) : RooFitResult() {
  //constructor that takes a string "x=y,a=b" etc ... copies those into constPars
  //This can be used for quickly checking what a TRooFit histogram looks like 
  //at a given parameter value ... e.g. 
  // h.Draw("param=2.0")

  //parse the finalPars expression for "x=y" terms
  RooArgList pars;
  
  TStringToken nameToken(constPars,",");
  while(nameToken.NextToken()) {
      TString subName = (TString)nameToken;
      //split around "=" sign
      TString parName = subName(0,subName.Index("="));
      TString parVal = subName(subName.Index("=")+1,subName.Length());
      
      RooRealVar* v = new RooRealVar(parName,parName,parVal.Atof());
      pars.addOwned(*v); //so that will delete when done

  }

  //init(pars);
  setConstParList(pars);
  setInitParList(RooArgList());
  setFinalParList(RooArgList());

}

void TRooFitResult::Draw(Option_t* option) {
  //Option: pull - draws pull plot for floating variables (provided their initial error was nonzero)

  TString opt(option);
  opt.ToLower();
  if(opt.Contains("pull")) {
    //only show float pars that had initial errors
    std::vector<int> indices;
    for(int i=0;i<_initPars->getSize();i++) {
      RooRealVar* rinit = static_cast<RooRealVar*>(_initPars->at(i));
      if(rinit->getError()) {indices.push_back(i);continue;}
      //check for constraintType ... can use that to get the initial error 
      if(rinit->getStringAttribute("constraintType")) {
        TString s(rinit->getStringAttribute("constraintType"));
        s.ToUpper();
        if(s=="STATPOISSON") {
          double tau = pow(TString(rinit->getStringAttribute("sumw")).Atof(),2)/TString(rinit->getStringAttribute("sumw2")).Atof();
          rinit->setError(1./sqrt(tau)); indices.push_back(i); continue;
        }
        else if(s=="NORMAL") { rinit->setError(1); indices.push_back(i); continue; }
        else if(s.BeginsWith("GAUSSIAN(")) {
          TString stddevStr = TString(s(s.Index(",")+1,s.Index(")")-(s.Index(","))-1));
          rinit->setError(stddevStr.Atof()); indices.push_back(i); continue;
        }
      }
      //before giving up, see if we can infer from name
      //HistFactory models have alpha parameters with range -5 to 5 ... error should be 1
      if(TString(rinit->GetName()).BeginsWith("alpha_") && fabs(rinit->getMin()+5.)<1e-9 && fabs(rinit->getMax()-5.)<1e-9) {
        rinit->setError(1); indices.push_back(i);
      }
    }
    if(indices.size()) {
      if(!fPullFrame) {
        fPullFrame = new TH1D("frame",GetTitle(),indices.size(),0,indices.size());fPullFrame->SetLineColor(0);
        fPullFrame->SetStats(0);
        for(uint i=0;i<indices.size();i++) {
          fPullFrame->GetXaxis()->SetBinLabel(i+1,floatParsFinal()[indices[i]].GetTitle());
        }
        fPullFrame->SetMaximum(3);fPullFrame->SetMinimum(-3);
        fPullFrame->GetYaxis()->SetTitle("(f-i)/#sigma_{i}");
        fPullLines.push_back(TLine(0,0,indices.size(),0));fPullLines.back().SetLineStyle(2);
        fPullBoxes.push_back(TBox(0,-2,indices.size(),2));fPullBoxes.back().SetFillColor(kYellow);
        fPullBoxes.push_back(TBox(0,-1,indices.size(),1));fPullBoxes.back().SetFillColor(kGreen);
      }
      if(!fPullGraph) {
        fPullGraph = new TGraphAsymmErrors;fPullGraph->SetMarkerStyle(20);
        for(uint i=0;i<indices.size();i++) {
          RooRealVar* r = static_cast<RooRealVar*>(_finalPars->at(indices[i]));
          RooRealVar* rinit = static_cast<RooRealVar*>(_initPars->at(indices[i]));
          double pull = (r->getVal()-rinit->getVal())/rinit->getError();
          double pullUp = (r->getVal()+r->getErrorHi()-rinit->getVal())/rinit->getError();
          double pullDown = (r->getVal()+r->getErrorLo()-rinit->getVal())/rinit->getError();
          fPullGraph->SetPoint(fPullGraph->GetN(),i+0.5,pull);
          fPullGraph->SetPointError(fPullGraph->GetN()-1,0,0,pull-pullDown,pullUp-pull);
          if(pullDown < fPullFrame->GetMinimum()) fPullFrame->SetMinimum(pullDown-1);
          if(pullUp > fPullFrame->GetMaximum()) fPullFrame->SetMaximum(pullUp+1);
        }
      }
    }
  }
  
  if (gPad) {
      if (!gPad->IsEditable()) gROOT->MakeDefCanvas();
      if (!opt.Contains("same")) {
         //the following statement is necessary in case one attempts to draw
         //a temporary histogram already in the current pad
         if (TestBit(kCanDelete)) gPad->GetListOfPrimitives()->Remove(this);
         gPad->Clear();
      }
  }
  
  TObject::Draw(option);
}
