
#include "TRooFit/TRooFitResult.h"
#include "RooRealVar.h"

#include "TPRegexp.h"

#include "TROOT.h"
#include "TPad.h"

#include <algorithm>

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
      TH1D* hh = 0; hh->Draw();
    } 
    else realpars.add(*arg);
  }
  
  setFinalParList(realpars);
  setConstParList(constpars);
  resetCovarianceMatrix(); //fills covariance matrix from float pars final
  setInitParList(floatParsFinal());
  
  removeFromDir(_dir); //don't put in the TDirectory
}

void TRooFitResult::adjustCovarianceMatrix() {
  //when the fit has put a parameter near the limit, the covariance matrix 
  //ends up being underestimated ... so we will compare the error to the diagonal in the covariance 
  //matrix and derive scale factors to be applied to the rows and columns ...
  
  
  RooFIter itr2( floatParsFinal().fwdIterator() );
  RooAbsArg* arg = 0;
  int i=0;
  
  std::vector<double> error_sf(floatParsFinal().getSize(),1.);
  
  TMatrixDSym cov( covarianceMatrix() );
  
  while( (arg = itr2.next()) ) {
    
    double bigErr = std::max(dynamic_cast<RooRealVar*>(arg)->getErrorHi(),-dynamic_cast<RooRealVar*>(arg)->getErrorLo());
    if(cov(i,i)==0) { cov(i,i)=pow(bigErr,2); } //when hessian was bad, we will recover the errors into the covariance matrix
    else {error_sf[i] = bigErr / sqrt(cov(i,i)); }
    i++;
  }
  
  for(int i=0;i<floatParsFinal().getSize();i++) {
    if(error_sf[i]>1.01) { //1.01 to cover rounding errors
      //inflating entries in ith row and column
      Warning("adjustCovarianceMatrix","Inflating %s covariance matrix entries by %g",floatParsFinal().at(i)->GetName(),error_sf[i]);
      for(int j=0;j<cov.GetNcols();j++) {
        cov(i,j) = cov(i,j)*error_sf[i];
        cov(j,i) = cov(j,i)*error_sf[i];
      }
    }
  }
  
  int covQualBackup = _covQual;
  setCovarianceMatrix( cov );
  _covQual = covQualBackup;
  
}


void TRooFitResult::resetCovarianceMatrix() {
  TMatrixDSym cov(floatParsFinal().getSize());
  
  RooFIter itr2( floatParsFinal().fwdIterator() );
  RooAbsArg* arg = 0;
  int i=0;
  while( (arg = itr2.next()) ) {
    cov(i,i) = pow(dynamic_cast<RooRealVar*>(arg)->getError(),2);
    if(cov(i,i)<1e-9) {
      //try to get error by checking for constraintType 
      //check for constraintType ... 
      if(arg->getStringAttribute("constraintType")) {
        TString s(arg->getStringAttribute("constraintType"));
        s.ToUpper();
        if(s=="STATPOISSON") {
          double tau = pow(TString(arg->getStringAttribute("sumw")).Atof(),2)/TString(arg->getStringAttribute("sumw2")).Atof();
          dynamic_cast<RooRealVar*>(arg)->setError(1./sqrt(tau));
          if(std::isnan(dynamic_cast<RooRealVar*>(arg)->getError())||std::isinf(dynamic_cast<RooRealVar*>(arg)->getError())) dynamic_cast<RooRealVar*>(arg)->setError(1e9);
        }
        else if(s=="NORMAL") { dynamic_cast<RooRealVar*>(arg)->setError(1);  }
        else if(s.BeginsWith("GAUSSIAN(")) {
          TString stddevStr = TString(s(s.Index(",")+1,s.Index(")")-(s.Index(","))-1));
          dynamic_cast<RooRealVar*>(arg)->setError(stddevStr.Atof());
        }
        cov(i,i) = pow(dynamic_cast<RooRealVar*>(arg)->getError(),2);
      }
      
    }
    i++;
  }
  int covQualBackup = _covQual;
  setCovarianceMatrix(cov);
  _covQual = covQualBackup;
}

TRooFitResult::TRooFitResult(const char* name, const char* title, const RooArgList& pars) 
  : RooFitResult(name,title) {
  
  init(pars);
  
}

TRooFitResult::TRooFitResult(const char* constPars) : RooFitResult() {
  //constructor that takes a string "x=y;a=b" etc ... copies those into constPars
  //Anything without an = sign in it will be interpreted as a floating parameter ... as long as it has "~" in it
  //This can be used for quickly checking what a TRooFit histogram looks like 
  //at a given parameter value ... e.g. 
  // h.Draw("param=2.0")

  //parse the finalPars expression for "x=y" terms
  RooArgList pars;
  RooArgList floats;
  
  TStringToken nameToken(constPars,";");
  while(nameToken.NextToken()) {
      TString subName = (TString)nameToken;
      //split around "=" sign
      if(subName.Contains("=")) {
        TString parName = subName(0,subName.Index("="));
        TString parVal = subName(subName.Index("=")+1,subName.Length());
        RooRealVar* v = new RooRealVar(parName,parName,parVal.Atof());
        pars.addOwned(*v); //so that will delete when done
      } else if(subName.Contains("~")) {
        subName.ReplaceAll("~","");
        //assume parameter is actually to belong to final pars
        RooRealVar* v =  new RooRealVar(subName,subName,0);
        v->setAttribute("injectValueAndError"); //used in TRooAbsH1::GetBinError to inject parameter values and errors
        floats.addOwned( *v );
      }
      

  }

  //init(pars);
  setConstParList(pars);
  setInitParList(RooArgList());
  setFinalParList(floats);

}

#include "TRegexp.h"
#include "TPRegexp.h"

void TRooFitResult::Draw(Option_t* option, const char* argFilter) {
  TRegexp pattern(argFilter,true);
  TPRegexp pattern2(argFilter);
  RooFIter itr = floatParsFinal().fwdIterator();
  RooArgList args;
  while(RooAbsArg* arg = itr.next() ) {
    TString aName(arg->GetName());
    if(aName.Contains(pattern) || pattern2.MatchB(aName)) args.add(*arg);
  }
  Draw(option,args);
}

double extractNumber(const TString& sWhat, const TString& pat) {
      double out = 0.05;
      int startPos = sWhat.Index(pat)+pat.Length();
      int endPos = (sWhat.Index(" ",sWhat.Index(pat))==-1) ? sWhat.Length() :   //goes to end of sWhat
                                sWhat.Index(" ",sWhat.Index(pat)); //goes to next occurance of a space 
      TString s(sWhat(startPos, endPos-startPos));
      if(s=="") return out; //defaults to 0.05;
      out = s.Atof();
      if(out==0||out>=1) return 0.05;
      return out;
}



void TRooFitResult::Draw(Option_t* option, const RooArgList& args) {
  //Option: pull - draws pull plot for floating variables (provided their initial error was nonzero)

  TString opt(option);
  opt.ToLower();
  if(opt.Contains("pull")) {
    //only show float pars that had initial errors
    std::vector<int> indices;
    for(int i=0;i<_initPars->getSize();i++) {
      RooRealVar* rinit = static_cast<RooRealVar*>(_initPars->at(i));
      if(args.getSize()&&!args.find(*rinit)) continue; //if args are specified, then must be in the list
      if(rinit->getError()) {indices.push_back(i);continue;}
      //check for constraintType ... can use that to get the initial error 
      if(rinit->getStringAttribute("constraintType")) {
        TString s(rinit->getStringAttribute("constraintType"));
        s.ToUpper();
        if(s=="STATPOISSON") {
          double tau = pow(TString(rinit->getStringAttribute("sumw")).Atof(),2)/TString(rinit->getStringAttribute("sumw2")).Atof();
          rinit->setError(1./sqrt(tau)); 
          if(std::isnan(rinit->getError())||std::isinf(rinit->getError())) rinit->setError(1e9);
          indices.push_back(i); continue;
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
        fPullFrame = new TH1D("frame",GetTitle(),indices.size(),0,indices.size());fPullFrame->SetLineColor(0);fPullFrame->SetDirectory(0);
        fPullFrame->SetStats(0);
        for(uint i=0;i<indices.size();i++) {
          fPullFrame->GetXaxis()->SetBinLabel(i+1,floatParsFinal()[indices[i]].GetTitle());
        }
        fPullFrame->SetMaximum(3);fPullFrame->SetMinimum(-3);
        fPullFrame->GetYaxis()->SetTitle("(f-i)/#sigma_{i}");
        /*fPullLines.push_back(TLine(0,0,indices.size(),0));*/fPullLines.push_back(TGraph(2));fPullLines.back().SetPoint(0,0,0);fPullLines.back().SetPoint(1,indices.size(),0);fPullLines.back().SetLineStyle(2);
        fPullBoxes.push_back(TGraphErrors(2));fPullBoxes.back().SetPoint(0,0,0);fPullBoxes.back().SetPoint(1,indices.size(),0);fPullBoxes.back().SetPointError(0,0,2);fPullBoxes.back().SetPointError(1,0,2);fPullBoxes.back().SetFillColor(kYellow);
        fPullBoxes.push_back(TGraphErrors(2));fPullBoxes.back().SetPoint(0,0,0);fPullBoxes.back().SetPoint(1,indices.size(),0);fPullBoxes.back().SetPointError(0,0,1);fPullBoxes.back().SetPointError(1,0,1);fPullBoxes.back().SetFillColor(kGreen);
        //fPullBoxes.push_back(TBox(0,-1,indices.size(),1));fPullBoxes.back().SetFillColor(kGreen);
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
  } else if(opt.Contains("cov")||opt.Contains("cor")) {
    //drawing a covariance matrix ..
    //obtain the reduced matrix ..
    RooArgList l;
    const RooArgList* a = &args;
    if(!args.getSize()) a = &floatParsFinal();
    
    RooAbsArg* singleArg = 0;
    if(a->getSize()==1) {
      //this means user wants to see the correlations or variance wrt to a single parameter ...
      singleArg = a->at(0);
      a = &floatParsFinal();
    }
    
    auto cov = reducedCovarianceMatrix(*a);
    TString myTitle = opt.Contains("cov")?"Covariance":"Correlation";
    if(opt.Contains("threshold=")) {
      //going to filter variables so only ones with at least one non-diagonal entry with magnitude greater than threshold are shown ...
      double thresh = extractNumber(opt,"threshold=");
      for(int i=1;i<=a->getSize();i++) {
        bool keep(false);
        for(int j=1;j<=a->getSize();j++) {
          if(i==j) continue; //skip diagonals
          if(singleArg && strcmp(a->at(j-1)->GetName(),singleArg->GetName())!=0) continue; //when doing a single arg, only look at that arg
          double val = opt.Contains("cor") ? cov(i-1,j-1)/sqrt(cov(i-1,i-1)*cov(j-1,j-1)) : cov(i-1,j-1);
          if(fabs(val) >= thresh) { keep = true; break; }
        }
        if(keep) l.add(*a->at(i-1));
      }
      if(l.getSize()!=a->getSize()) {
        if(l.getSize()==0) {
          Info("Draw","No values above threshold %g",thresh);
          return;
        }
        if(singleArg) l.add(*singleArg);
        a = &l;
        cov.ResizeTo(a->getSize(),a->getSize());
        cov = reducedCovarianceMatrix(*a);
        myTitle += Form(" (Threshold = %g)",thresh);
      }
    }
    
    TH1* covHist = 0;
    if(singleArg) {
      myTitle = Form("%s %s",singleArg->GetTitle(),myTitle.Data());
      covHist = new TH1D(opt.Contains("cov") ? "covariance" : "correlation",myTitle,a->getSize()-1,0,a->getSize()-1);
      covHist->SetDirectory(0);
      int j = a->index(singleArg)+1;
      std::vector<std::pair<std::string,double>> vals;
      int count(0);
      for(int i=1;i<=a->getSize();i++) {
        if(i==j) continue; //don't plot self correlation/covariance
        vals.push_back( std::make_pair(a->at(i-1)->GetTitle(), (opt.Contains("cor") ? cov(i-1,j-1)/sqrt(cov(i-1,i-1)*cov(j-1,j-1)) : cov(i-1,j-1) ) ) );
        count++;
        covHist->GetXaxis()->SetBinLabel(count,a->at(i-1)->GetTitle());
        if(opt.Contains("cor")) {
          covHist->SetBinContent(count,j,cov(i-1,j-1)/sqrt(cov(i-1,i-1)*cov(j-1,j-1)));
        } else {
          covHist->SetBinContent(count,j,cov(i-1,j-1));
        }
      }
      if(opt.Contains("cor")) { covHist->SetAxisRange(-1,1,"Y"); }
      covHist->SetBarWidth(0.9);covHist->SetBarOffset(0.05);
      covHist->SetFillColor(4);
      
      std::sort(vals.begin(),vals.end(), [](auto &left, auto &right) { return fabs(left.second) > fabs(right.second); });
      for(unsigned int i=0;i<vals.size();i++) {
        covHist->GetXaxis()->SetBinLabel(i+1,vals[i].first.c_str());
        covHist->SetBinContent(i+1,vals[i].second);
      }
      
      //create a copy histogram and add it ...
      TH1* copyHist = static_cast<TH1*>(covHist->Clone(covHist->GetName())); copyHist->SetDirectory(0);
      for(int i=1;i<=copyHist->GetNbinsX();i++) copyHist->SetBinContent(i,-copyHist->GetBinContent(i));
      copyHist->SetFillColor(kCyan);
      covHist->GetListOfFunctions()->Add(copyHist,"b same");
      
    } else {
      covHist = new TH2D(opt.Contains("cov") ? "covariance" : "correlation",myTitle,a->getSize(),0,a->getSize(),a->getSize(),0,a->getSize());
      covHist->SetDirectory(0);
      for(int i=1;i<=a->getSize();i++) {
        covHist->GetXaxis()->SetBinLabel(i, a->at(i-1)->GetTitle());
        for(int j=1;j<=a->getSize();j++) {
          if(i==1) covHist->GetYaxis()->SetBinLabel(j, a->at(j-1)->GetTitle());
          if(opt.Contains("cor")) {
            covHist->SetBinContent(i,j,cov(i-1,j-1)/sqrt(cov(i-1,i-1)*cov(j-1,j-1)));
          } else {
            covHist->SetBinContent(i,j,cov(i-1,j-1));
          }
        }
      }
      if(opt.Contains("cor")) { covHist->SetAxisRange(-1,1,"Z"); }
    }
    
    if(fCovHist) delete fCovHist;
    fCovHist = covHist;
    fCovHist->SetStats(0);
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
  if(!opt.Contains("same") && opt.Contains("pull") && fPullFrame) fPullFrame->Draw();
  if(!opt.Contains("same") && (opt.Contains("cov")||opt.Contains("cor"))) {
    if(fCovHist->GetDimension()==2) {
      if(opt.Contains("text")) fCovHist->Draw("COLZ TEXT");
      else fCovHist->Draw("COLZ");
    } else if(fCovHist->GetDimension()==1) {
      if(opt.Contains("text")) fCovHist->Draw("B TEXT");
      else fCovHist->Draw("B");
    }
  }
  TObject::Draw(option);
}
