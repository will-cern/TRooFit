
#include "TRooFit/TRooAbsHStack.h"

#include "TMath.h"

#include "TRooFit/TRooHStack.h"

using namespace std;


TRooAbsHStack::~TRooAbsHStack() {
 SafeDelete(fStack); //FIXME: should delete the contained hists too!
}

TRooH1* TRooAbsHStack::firstHistogram() {
  compIter()->Reset();
  TObject* obj = 0;
  while( (obj = compIter()->Next()) ) {
    if(obj->InheritsFrom(TRooH1::Class())) return static_cast<TRooH1*>(obj);
  }
  return 0;
}

Double_t TRooAbsHStack::missingEvents() const {
  double out = 0;
  RooFIter funcIter = compList().fwdIterator() ;
  RooAbsReal* func ;
  while((func=(RooAbsReal*)funcIter.next())) {
    if(func->InheritsFrom(TRooAbsH1::Class())) {
      out += (dynamic_cast<TRooAbsH1*>(func))->missingEvents();
    }
  }
  return out;
}

void TRooAbsHStack::Add(RooAbsReal* func) {
  if(this->IsA() == TRooHStack::Class()) {
    //must be extendable to be added
    if(func->InheritsFrom(RooAbsPdf::Class()) && !static_cast<RooAbsPdf*>(func)->canBeExtended()) {
      Error("Add","Cannot add an unextended pdf to a TRooHStack. Please wrap it in a RooExtendPdf first");
      return;
    }
  }
  compList().add(*func);
  reinit();
}

void TRooAbsHStack::Add(TRooH1* hist) {
  //if this is the first one, take observables from it 
  if(!firstHistogram()) {
    fObservables.add(hist->fObservables);  
    fDummyHist = static_cast<TH1*>(hist->fDummyHist->Clone(GetName()));
    fDummyHist->SetTitle(GetTitle());
    //FIXME .. only needed this for GetHistogram method .. should remove that dependency
    //fHists.push_back((TH1*)hist->fHists[0]->Clone(GetName()));fHists[0]->Reset();fHists[0]->SetDirectory(0);
    SetRangeName( hist->GetName() );
  }

  compList().add(*hist);
  compList().setName("hists"); //so that it now shows up in Print
  //claim all statFactors as my own, if I can
  //binning must be consistent
  
  bool consistent(true);
  try {
    //TH1::CheckConsistency(hist->fHists[0],fHists[0]);
  } catch(...) {
    consistent = false;
  }
  
  
  if(hist->fStatFactors.isOwning() && consistent) {
    hist->fStatFactors.releaseOwnership();
    hist->fStatFactors.setName("!statFactors"); //hide from prints
    
    //need to combine common stat factors before adding these new ones
    //this includes updating hist->fShapeFactors to reference the new, common, statFactor
    for(int i=0;i<fStatFactors.getSize();i++) {
      RooAbsArg* ss = fStatFactors.at(i);
      int binNum = TString(ss->getStringAttribute("statBinNumber")).Atoi();
      for(int j=0;j<hist->fStatFactors.getSize();j++) {
        RooAbsArg* s = hist->fStatFactors.at(j);
        if(binNum != TString(s->getStringAttribute("statBinNumber")).Atoi()) continue;
        //found a match, need to replace!
        //combine sumw and sumw2 attributes
        Info("Add","Combining StatFactors: Replacing %s with %s in %s",s->GetName(),ss->GetName(),hist->GetName());
        hist->fStatFactors.replace(*s,*ss);
        hist->fShapeFactors.replace(*s,*ss);
        ss->setStringAttribute("sumw",Form("%f",(TString(ss->getStringAttribute("sumw")).Atof() + TString(s->getStringAttribute("sumw")).Atof())));
        ss->setStringAttribute("sumw2",Form("%f",(TString(ss->getStringAttribute("sumw2")).Atof() + TString(s->getStringAttribute("sumw2")).Atof())));
        delete s;
        break;
      }
    }
    //rename the vars to match my name, to avoid confusion, and add them to my stat factors
    RooFIter iter = hist->fStatFactors.fwdIterator();
    while( RooAbsArg* arg = iter.next() ) {
      if(fStatFactors.find(*arg)) continue;
      TString oldName(arg->GetName()); //dont correct the ones we already replaced above
      TString newName(oldName.Replace(0,strlen(hist->GetName()),GetName()));
      arg->SetName(newName);
      fStatFactors.addOwned(*arg);
    }
    if(fStatFactors.getSize()) fStatFactors.setName("statFactors");
  } else {
    Info("Add","Cannot acquire statFactors of %s",hist->GetName());
  }
  
  reinit();
}

THStack* TRooAbsHStack::GetStack(const RooFitResult* fr) const {
  if(!fStack) { 
    fStack = new THStack(Form("%s_stack",GetName()),GetTitle()); 
    //transfer min and max to stack 
    fStack->SetMinimum(GetMinimum());
    fStack->SetMaximum(GetMaximum());
  }
  fillStack( fStack, fr, false );
  return fStack;
}


THStack* TRooAbsHStack::fillStack(THStack* stack, const RooFitResult* fr, bool noRestyle) const {
  TList existingHists;
  
  
  //clears the existing stack ...but we will reuse hists if we can!
  if(stack->GetHists()) {
    existingHists.AddAll(stack->GetHists());
    stack->GetHists()->Clear(); 
  }
   

  RooFIter funcIter = compList().fwdIterator() ;
  RooAbsReal* func ;
  int i=0;

  while((func=(RooAbsReal*)funcIter.next())) {
    TH1* hist = (existingHists.GetSize()>i) ? (TH1*)(existingHists.At(i)) : 0;
    if(func->InheritsFrom(TRooH1::Class())) {
      static_cast<TRooH1*>(func)->SetRangeName(GetRangeName()); //FIXME: perhaps should set this when we acquire the hist
      if(!(hist&&noRestyle)) hist = static_cast<TRooH1*>(func)->createOrAdjustHistogram(hist);
      (static_cast<TRooH1*>(func))->TRooAbsH1::fillHistogram(hist,fr,false); //fills the hist
    } else if(func->InheritsFrom(RooAbsPdf::Class())) {
      //need to fill a histogram ... use binning of GetRangeName() 
       
       
       
       if(!(hist&&noRestyle)) hist = createOrAdjustHistogram(hist); 
      
      //FIXME: need to use fit result to move parameters!
      
      RooAbsPdf* pdf = (RooAbsPdf*)func;
      //use cdf to fill bins ... i.e. correctly integrate over pdfs 
      double expectedEvents = pdf->expectedEvents(fObservables); //will scale to this
      
      RooAbsReal* cdf = pdf->createCdf(fObservables);
      
      std::vector<double> prevIntegral(hist->GetNbinsY(),0.); //need one running total for each 'row' if doing up to 2D
      for(int i=1;i<=hist->GetNbinsX();i++) {
        if(fObservables.getSize()>0) {
          if(fObservables[0].InheritsFrom(RooAbsRealLValue::Class())) {
            ((RooAbsRealLValue&)fObservables[0]).setVal( hist->GetXaxis()->GetBinUpEdge(i) );
          } else {
            //discrete ... FIXME
          }
        }
        for(int j=1;j<=hist->GetNbinsY();j++) { //max 2D support for now 
          if(fObservables.getSize()>1) {
            if(fObservables[1].InheritsFrom(RooAbsRealLValue::Class())) {
              ((RooAbsRealLValue&)fObservables[1]).setVal( hist->GetYaxis()->GetBinUpEdge(j) );
            } else {
              //discrete ... FIXME
            }
          }
          
          double nextIntegral = cdf->getVal();
          hist->SetBinContent(hist->GetBin(i,j,0),(nextIntegral-prevIntegral[j-1])*expectedEvents);
          prevIntegral[j-1]=nextIntegral; 
        }
      }
      //func->fillHistogram(hist,fObservables); //this method assumes function is flat across each bin too (like getBinContent of TRooH1)
      
      
      delete cdf;
      
    } else {
      std::cout << "NOT SUPPORTED!!" << std::endl;
      
    }
    if(hist) stack->Add(hist);
    i++;
  }
  
  return stack;
}


#include "TROOT.h"
#include "TPad.h"
void TRooAbsHStack::Draw(Option_t *option)
{
   TRooAbsHStack::Draw(TRooFitResult(option),option);
}


void TRooAbsHStack::Draw(const TRooFitResult& r, Option_t* option) {
  TString opt = option;
  opt.ToLower();
  

  TRooFitResult* r2 = 0;
  bool hadInit=false;
  
  if(r.floatParsFinal().getSize()>0 || r.constPars().getSize()>0) { //wont have any pars if was doing a straightforward draw
    if(opt.Contains("init")) {
      //request to draw initial parameters instead of final
      r2 = new TRooFitResult(r.floatParsInit());
      opt.ReplaceAll("init","");hadInit=true;
    } else {
      r2 = new TRooFitResult(r.floatParsFinal());
    }
    r2->setConstParList(r.constPars());
  }
  
      // Draw this stack ... potentially as a hist!
  
   bool found(false);
   TObject* me = dynamic_cast<TObject*>(this);
   
   if (gPad) {
      if (!gPad->IsEditable()) gROOT->MakeDefCanvas();
      if (!opt.Contains("same")) {
         //the following statement is necessary in case one attempts to draw
         //a temporary histogram already in the current pad
         if (me->TestBit(kCanDelete)) gPad->GetListOfPrimitives()->Remove(me);
         gPad->Clear();
         
         //also delete any DrawnStacks that matched this pad:
         auto itr = fDrawStacks.begin();
         while( itr != fDrawStacks.end() ) {
          if(itr->pad==gPad) {
            SafeDelete(itr->stack);
            SafeDelete(itr->fr);
            fDrawStacks.erase(itr);
          } else {
            ++itr;
          }
         }
      } else {
        //check if I'm already in the list of primitives ... if so, we wont add me a second time
        if(gPad->GetListOfPrimitives()->FindObject(me)) {
          gPad->Modified(true);
          found = true;
        }
      }
   }
   if(!found) me->AppendPad(opt.Data()); //will create gPad
   
   if(!gPad->IsEditable()) return;
   
   if(opt.Contains("pdf")) {
    //FIXME: ideally would draw a stack of pdfs (in a TMultiGraph?)
    TRooAbsH1::Draw(r,opt);
   
   } else if(!opt.Contains("hist")) {
    fDrawStacks.emplace_back( DrawnStack() );
    fDrawStacks.back().pad = gPad;
    fDrawStacks.back().stack = new THStack(Form("%s_stack",GetName()),GetTitle());
    //transfer min and max to stack 
    fDrawStacks.back().stack->SetMinimum(GetMinimum());
    fDrawStacks.back().stack->SetMaximum(GetMaximum());
    int fillType=0;
    if(opt.Contains("e3")) {
        fillType = TString(opt(opt.Index("e3")+1,opt.Length())).Atoi();
        if(fillType>=3000 && fillType<=3999) {
          opt.ReplaceAll(TString::Format("e%d",fillType),"");
        } else {
          fillType=0;
        }
    }
    std::cout << "opt = " << opt << "," << fillType << std::endl;
    
    fDrawStacks.back().opt = opt;
    fillStack(fDrawStacks.back().stack,r2,false/* noRestyle*/);
    fDrawStacks.back().fr = r2;
    
    gPad->GetListOfPrimitives()->Add( fDrawStacks.back().stack , opt ); //adding stack directly because cant figure out how to get clickable axis without it
    
    //if drawing with option "e3XXX" then need to also draw as a histogram
    if(fillType) {
        TRooAbsH1::Draw(r,TString::Format("%s e%dsame",(hadInit)?"init":"",fillType));
        //FIXME: would like to have stack's maximum match up to error bar maximum
    }
    
   } else {
    //not drawing the stack, pass onto parent class to draw as a hist instead
    opt.ReplaceAll("hist","");
    TRooAbsH1::Draw(r,opt);
   }
  

}

void TRooAbsHStack::Paint(Option_t* option) {
 
    //update the min and max on the last stack
    if(fDrawStacks.size()) {
       //fDrawStacks.back().stack->SetMinimum(GetMinimum());
       //fDrawStacks.back().stack->SetMaximum(GetMaximum()); //NOTE: maybe move this min max setting into 'styleStack' function
    }
 
    //fill the stacks and paint them
    
//     uint i = 0;
//     for(auto& stack : fDrawStacks) {
//       if(stack.pad == gPad) {  
//         fillStack(stack.stack,stack.fr,!(i==fDrawStacks.size()-1) /* restyles the last stack only */); //updates with current values
//         //stack.stack->Paint(stack.opt); ..don't paint here because added histogram to primitives directly (above)
//       }
//       i++;
//     }
 
    //GetStack(0,true); fStack->Paint(option); 
    
    
    //draw any histograms we painted here too!
    TRooAbsH1::Paint(option);
    
}
