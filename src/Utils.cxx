
#include "TRooFit/Utils.h"

#include "RooMinimizer.h"

TFile* TRooFit::m_debugFile = 0;
TObject* TRooFit::m_msgObj = 0;

std::map<TString,Int_t> TRooFit::m_colorsByName;

ClassImp(TRooFit)

TRooFit::TRooFit() : TObject() { }


double TRooFit::significance(double obs, double exp, double relUncert, double relUncertDown) {
  double r = (obs<exp) ? relUncertDown : relUncert;
  double v_b = (r*exp)*(r*exp); //variance of b i.e. sigma_b^2
  
  if(fabs(exp)<1e-9 && obs==0) return 0;
  
  if(obs >= exp) {
    //Cowan's formula based on Wilks theorem approximation of one-sided likelihood ratio test statistic
    return sqrt( 2.0*( obs*((obs!=0)?log(obs*(exp+v_b)/(exp*exp + obs*v_b)):0) - (exp*exp/v_b)*log(1.0 + v_b*(obs-exp)/(exp*(exp+v_b)))) );
  } else {
    //Uncapped test statistic case, I believe just leads to a negative sign for the significance
    return -sqrt( 2.0*( obs*((obs!=0)?log(obs*(exp+v_b)/(exp*exp + obs*v_b)):0) - (exp*exp/v_b)*log(1.0 + v_b*(obs-exp)/(exp*(exp+v_b)))) );
  }
}

RooAbsPdf* TRooFit::BuildModel(TRooAbsH1& pdf, RooAbsData& data) {
  return pdf.buildConstraints( *data.get() , "", true );
}

Int_t TRooFit::GetColorByName(const char* name, bool generateColor) {
  if(m_colorsByName.find(name)!=m_colorsByName.end()) return m_colorsByName[name];
  if(!generateColor) return 0;
  
  //get list of all colors 
  std::set<int> allCols;
  for(auto c : m_colorsByName) allCols.insert(c.second);
  
  
  int nextColor=2; //start at 2 so that never return 1 (which is kBlack)
  if(allCols.size()) {
    auto colItr = allCols.begin();
    int lastCol = (*colItr);
    colItr++;
    while( colItr!=allCols.end() && lastCol+1 == (*colItr) ) {lastCol=(*colItr); colItr++;}
    nextColor=lastCol+1;
  }
  
  m_colorsByName[name] = nextColor;
  
  return nextColor;
  
}

void TRooFit::RemoveErrors(TH1* hist, double relUncertThreshold) {
  for(int i=1;i<=hist->GetXaxis()->GetNbins();i++) {
    for(int j=1;j<=hist->GetYaxis()->GetNbins();j++) {
      for(int k=1;k<=hist->GetZaxis()->GetNbins();k++) {
        if(hist->GetBinError(i,j,k)/(hist->GetBinContent(i,j,k)+1e-15)<relUncertThreshold) hist->SetBinError(i,j,k,0.);
      }
    }
  }

}


RooStats::ModelConfig* TRooFit::CreateModelConfig(RooWorkspace& w, const char* modelName, const char* dataName, const char* poiName ) {

  RooAbsData* data = w.data(dataName);
  if(!data) {
    std::cout << "unable to find data" << std::endl;
    return 0;
  }
  RooAbsPdf* model = w.pdf(modelName);
  if(!model) {
    std::cout << "unable to find model" << std::endl;
    return 0;
  }
  RooRealVar* poi = w.var(poiName);
  if(!poi) {
    std::cout << "unable to find parameter of interest" << std::endl;
    return 0;
  }
  

  RooStats::ModelConfig* mc = new RooStats::ModelConfig("ModelConfig",&w);
  mc->SetPdf(*model);
  
  
   RooArgSet* obs = model->getObservables(*data);
   mc->SetObservables(*obs);
   delete obs;


   RooArgSet* a = model->getObservables(*poi);
   mc->SetParametersOfInterest(*a);
   delete a;


   
   //infer the global observables, nuisance parameters, model args (const np) 
   RooArgSet* gobs_and_np = model->getParameters(*data);

   //remove the poi ...
   gobs_and_np->remove(*mc->GetParametersOfInterest());

   RooArgSet gobs;
   gobs.add(*gobs_and_np); //will remove the np in a moment


   //now pass this to the getAllConstraints method ... it will modify it to contain only the np!
   RooArgSet* s = model->getAllConstraints(*mc->GetObservables(),*gobs_and_np);
   delete s; //don't ever need this 

   //gobs_and_np now only contains the np
   gobs.remove(*gobs_and_np);
   //ensure all global observables are held constant now - important for pll evaluation
   gobs.setAttribAll("Constant",kTRUE);

   mc->SetGlobalObservables(gobs);

   //finally, split out the constant np ... these are the 'model arguments' 
   RooArgSet np;RooArgSet args;
   std::unique_ptr<TIterator> itr(gobs_and_np->createIterator());
   RooAbsArg* arg = 0;
   while((arg=(RooAbsArg*)itr->Next())) { 
      if(arg->isConstant()) args.add(*arg);
      else np.add(*arg);
   }
   delete gobs_and_np;
   
   mc->SetNuisanceParameters(np);

   return mc;

}

void TRooFit::setRecommendedDefaultOptions() {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
}

RooAbsReal* TRooFit::createNLL(RooAbsPdf* pdf, RooAbsData* data, const RooArgSet* gobs) {
  //need to specify global observables so that they are included in the normalization set 
  //Use the offset option because it improves fitting convergence, but consequence is that the 
  //value of the FCN stored in the FitResult will not correspond to the actual NLL, it corresponds to the offset value
  //the offset is determined by the initial (first) value of the NLL function
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); 
  RooAbsReal* out = pdf->createNLL( *data, RooFit::GlobalObservables(*gobs), RooFit::Offset(1) );
  RooMsgService::instance().setGlobalKillBelow(msglevel); 
  return out;
}

std::pair<RooAbsData*,RooArgSet*> TRooFit::generateAsimovDataset(RooAbsPdf* thePdf, RooAbsData* data, const RooArgSet* gobs) { //should only depend on poi_prime value and args 
  RooArgSet* obs = thePdf->getObservables(data);
  auto out = generateAsimovDataset(thePdf, obs, gobs);
  delete obs;
  return out;
}

//should perform a conditional fit before computing dataset
std::pair<RooAbsData*,RooArgSet*> TRooFit::generateAsimovDataset(RooAbsPdf* thePdf, const RooArgSet* obs, const RooArgSet* gobs) { //should only depend on poi_prime value and args

  std::pair<RooAbsData*,RooArgSet*> out = std::make_pair(nullptr,nullptr);
  
  

   //Get the Asimov Dataset!
   out.first = thePdf->generate(*obs,RooFit::ExpectedData(),RooFit::Extended()); //must specify we want extended ... because if it's a RooSimultaneous, RooSimSplitGenContext wont generate non-integer weights unless 'extended'
   /*
   if(thePdf->InheritsFrom(RooSimultaneous::Class())) {
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(thePdf);
    TIterator* iter = simPdf->indexCat().typeIterator();
    RooCatType* tt;
    while( (tt=(RooCatType*)iter->Next()) ) {
      RooAbsPdf* subPdf = simPdf->getPdf(tt->GetName());
      
      //adjust all binnings to match boundaries if needed ... 
      RooFIter obsItr = obs->fwdIterator();
      while( RooAbsArg* obsArg = obsItr.next() ) {
        if(!subPdf->isBinnedDistribution(*obsArg)) continue;
        
      }
      
    }
    delete iter;
    //iterate over pdfs ...
    
   } else {
    out.first = thePdf->generate(*obs,RooFit::ExpectedData(),RooFit::Extended()); //must specify we want extended ... because if it's a RooSimultaneous, RooSimSplitGenContext wont generate non-integer weights unless 'extended'
   }*/
   
   
   if(gobs) {
      RooArgSet* asimov_gobs = new RooArgSet("asimov_gobs");
      gobs->snapshot(*asimov_gobs);
      out.second = asimov_gobs;
      
      //get all floating parameters ... 
      RooArgSet* allPars = thePdf->getParameters(obs);
      RooAbsCollection* floatPars = allPars->selectByAttrib("Constant",kFALSE);
      
      
      RooArgSet temp; temp.add(*floatPars);/*use a copy since don't want to modify the np set*/
      RooArgSet* constraints = thePdf->getAllConstraints(*obs,temp);
    
      //loop over constraint pdfs, recognise gaussian, poisson, lognormal
      TIterator* citer = constraints->createIterator();
      RooAbsPdf* pdf = 0;
      while((pdf=(RooAbsPdf*)citer->Next())) { 
          //determine which gobs this pdf constrains. There should only be one!
          std::unique_ptr<RooArgSet> cgobs(pdf->getObservables(out.second));
          if(cgobs->getSize()!=1) { std::cout << "constraint " << pdf->GetName() << " constrains " << cgobs->getSize() << " global observables. skipping..." << std::endl; continue; }
          RooAbsArg* gobs_arg = cgobs->first();
    
          //now iterate over the servers ... the first server to depend on a nuisance parameter we assume is the correct server to evaluate ..
          std::unique_ptr<TIterator> itr(pdf->serverIterator());
          for(RooAbsArg* arg = static_cast<RooAbsArg*>(itr->Next()); arg != 0; arg = static_cast<RooAbsArg*>(itr->Next())) {
            RooAbsReal * rar = dynamic_cast<RooAbsReal *>(arg); 
            if( rar && ( rar->dependsOn(*floatPars) ) ) {
                out.second->setRealValue(gobs_arg->GetName(),rar->getVal()); //NOTE: We could add some checks that all gobs actually get set
            }
          }
      }
      delete citer;
    
      delete constraints;
      
      delete floatPars;
      delete allPars;
      
      
      
   }
   
   return out;

}

std::pair<RooAbsData*,RooArgSet*> TRooFit::generateToy(RooAbsPdf* thePdf, RooAbsData* data, const RooArgSet* gobs, bool doBinned) {
  RooArgSet* obs = thePdf->getObservables(data);
  auto out = generateToy(thePdf, obs, gobs, doBinned);
  delete obs;
  return out;
}

std::pair<RooAbsData*,RooArgSet*> TRooFit::generateToy(RooAbsPdf* thePdf, const RooArgSet* obs, const RooArgSet* gobs, bool doBinned) {
   
   
   std::pair<RooAbsData*,RooArgSet*> out = std::make_pair(nullptr,nullptr);
   
   if(gobs) {
      RooArgSet* toy_gobs = new RooArgSet("toy_gobs");
      gobs->snapshot(*toy_gobs);

      //ensure we use the gobs from the model ...
      RooArgSet* model_gobs = thePdf->getObservables(*gobs);
      RooDataSet* globals = thePdf->generateSimGlobal(*model_gobs,1);
      *toy_gobs = *globals->get(0);
      delete globals; delete model_gobs;
      out.second = toy_gobs;
    }
   
    //have to first check if expected events are zero ... if it is, then the dataset is just empty ('generate' method fails when expectedEvents = 0)
    
    if(thePdf->expectedEvents(*obs)==0) {
        out.first = new RooDataSet("toyData","toyData",*obs);
    } else {
      
      out.first = (doBinned) ? thePdf->generate(*obs,RooFit::Extended(),RooFit::AllBinned()) : thePdf->generate(*obs,RooFit::Extended());
    }
    
    
   return out;

}

RooFitResult* TRooFit::minimize(RooAbsReal* nll, bool save, bool hesse) {
  
  
  int printLevel  =   ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  TString minim=      ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
  TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
  int strategy =      ::ROOT::Math::MinimizerOptions::DefaultStrategy();
  double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
  
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if(printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); 
  
  RooMinimizer _minimizer(*nll);
  _minimizer.optimizeConst(2);
  _minimizer.setPrintLevel(printLevel);
  _minimizer.setMinimizerType(minim);
  _minimizer.setStrategy(strategy);
  _minimizer.setEps( tol );
  
  int status = 0;
  for(int tries=1,maxtries=4; tries <= maxtries; ++tries) {
    status = _minimizer.minimize(minim.Data(),algorithm.Data());
    //int status = _minimizer->migrad();
    if(status%1000 == 0) break; //fit was good
    if(tries >= maxtries) break; //giving up

    //NOTE: minuit2 seems to distort the tolerance in a weird way, so that tol becomes 100 times smaller than specified
    //Also note that if fits are failing because of edm over max, it can be a good idea to activate the Offset option when building nll
    msg().Warning("minimize","Fit Status=%d (edm=%f, tol=%f), Rescanning #%d...",status,_minimizer.fitter()->Result().Edm(),tol,tries);
    _minimizer.minimize(minim,"Scan");
    if(tries == 2) _minimizer.setStrategy(strategy+1); //up the strategy
    if(tries == 3) { minim = "Minuit"; algorithm = "migradImproved"; }
  }
  
  if(hesse) _minimizer.hesse(); //note: I have seen that you can get 'full covariance quality' without running hesse ... is that expected?
  
  RooFitResult* out = 0;
  if(save) out = _minimizer.save("fitResult","fitResult");
  
  if(status%1000 != 0) {
    msg().Warning("minimize"," fit status %d",status);
  }
  
  if(printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel); 
  
  return out;
}

//The reason for this minos method is because the builtin ROOT method seems to be both inefficient and also fails to find solutions (e.g. runs out of iterations)
RooFitResult* TRooFit::minos(RooAbsReal* nll, const RooArgSet& pars, RooFitResult* unconditionalFitResult, bool respectBoundaries) {
  
  int printLevel  =   ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  
  const double minosPrecision = 0.002;
  
  //for each parameter will search the 2*pll distribution for where it equals +1/-1 
  //search assumes parabolic shape
  
  RooFIter itr = pars.fwdIterator();
  
  if(unconditionalFitResult==0) unconditionalFitResult = TRooFit::minimize(nll,true);
  
  RooArgSet* floatPars = nll->getObservables(unconditionalFitResult->floatParsFinal());
  
  //ensure all other parameters are held constant (at the fit result's constPar values if necesary)
  RooArgSet* otherPars = nll->getParameters(unconditionalFitResult->floatParsFinal());
  *otherPars = unconditionalFitResult->constPars();
  otherPars->setAttribAll("Constant",kTRUE);
  delete otherPars;
  
  //start by moving all float pars to best fit values and ensuring they are floating for the fit
  *floatPars = unconditionalFitResult->floatParsFinal();
  floatPars->setAttribAll("Constant",kFALSE);
  
  double nll_absMin = nll->getVal();//unconditionalFitResult->minNll(); -- cannot use minNll because might be offset corrected
  
  int status(0);
  
  while( RooAbsArg* arg = itr.next() ) {
    RooRealVar* par = dynamic_cast<RooRealVar*>(arg);
    if(!par) continue; //only do minos errors for real vars 
    par = dynamic_cast<RooRealVar*>(floatPars->find(*par)); //ensures we are using the actual version of the variable in nll
    if(!par) {
      std::cout << Form("Warning: %s is not floating in the fit, so no minos error will be computed",par->GetName()) << std::endl;
      continue; //only do minos errors on pars that were floating in unconditional fit
    }

    
    
    //starting guess will be existing error on parameters
    //can use par since above assignment operation should have copied over all the errors
    double bestFitVal = par->getVal();
    double parMax = par->getMax(); double parMin = par->getMin();
    double bestFitErrHi = par->getErrorHi();
    double bestFitErrLo = par->getErrorLo();
   
    double sigma_hi = findSigma(nll, nll_absMin, par, bestFitVal+bestFitErrHi, bestFitVal, 1, minosPrecision);
    double sigma_lo = findSigma(nll, nll_absMin, par, bestFitVal+bestFitErrLo, bestFitVal, -1, minosPrecision);
    
    if(printLevel>0) { msg().Info("TRooFit::minos","%s +/- errors are %f / %f",arg->GetName(),sigma_hi,sigma_lo); }
    
    if(std::isnan(sigma_hi)) {
      sigma_hi = 0;
      unconditionalFitResult->floatParsFinal().find(*par)->setAttribute("minosErrorHi_isnan");
      status=3;
    }
    if(std::isnan(sigma_lo)) {
      sigma_lo = 0;
      unconditionalFitResult->floatParsFinal().find(*par)->setAttribute("minosErrorLo_isnan");
      status=3;
    }
    
    if(respectBoundaries) {
      if(sigma_hi+bestFitVal>parMax) sigma_hi=parMax-bestFitVal;
      if(bestFitVal-sigma_lo<parMin) sigma_lo=bestFitVal-parMin;
    }
    
    //copy the errors into the fit result
    ((RooRealVar*)unconditionalFitResult->floatParsFinal().find(*par))->setAsymError(-sigma_lo,sigma_hi);
    
    //finish by restoring floatPars to unconditional values
    *floatPars = unconditionalFitResult->floatParsFinal();
    
    //and ensure range is restored 
    par->setMin(parMin); par->setMax(parMax);
    
  }
  unconditionalFitResult->_statusHistory.push_back( std::make_pair("TF_MINOS", status ) );
  
  delete floatPars;
  
  return unconditionalFitResult;
}

//this function makes a series of minos calls, each time it promotes parameters in the given group from floating to const
//parameters are flagged in groups according to attributes (the group name is an attribute name)
//the first element of the return is always the unconditionalFitResult
std::vector<RooFitResult*> TRooFit::minos_series(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult, bool runInitialMinos) {

  if(unconditionalFitResult==0) runInitialMinos=true; //calling minos will create the unconditional fit result for us

  if(runInitialMinos) {
    std::cout << Form("minos_series: running unconditional fit and (global) minos errors for %d parameters",pars.getSize()) << std::endl;
    unconditionalFitResult = TRooFit::minos(nll,pars,unconditionalFitResult);
  }

  std::vector<RooFitResult*> out;
  
  out.push_back(unconditionalFitResult);
  
  RooArgSet* floatPars = nll->getObservables(unconditionalFitResult->floatParsFinal());
  
  RooFitResult* previousResult = unconditionalFitResult;
  
  for(auto& group : groups) {
    //create our RooFitResult, which is a copy of the last fit result
    RooFitResult* nextResult = new RooFitResult(*previousResult);
    
  
    RooAbsCollection* parsInGroup = floatPars->selectByAttrib(group,kTRUE);
    
    RooArgSet floatResults;
    RooFIter itr = nextResult->floatParsFinal().fwdIterator();
    while( RooAbsArg* arg = itr.next() ) {
      if(parsInGroup->find(arg->GetName())) floatResults.add(*arg); //track the actual args here so that we can call "remove" later
    }
    
    if(floatResults.getSize()==0) {
      //no floating left
      std::cout << Form("minos_series: no float parameters in group %s (may have been covered by previous groups in series)",group.Data()) << std::endl;
    } else {
      if(parsInGroup->getSize()==0) {
        msg().Warning("minos_series","No parameters found in group %s",group.Data());
      }
      //make selected parameters const, and move them from floatPars to constPars in the fit result
      parsInGroup->setAttribAll("Constant",kTRUE);
      //must create a copy since remove below will delete (because floatParsFinal is an owning collection)
      RooArgSet floatCopy; floatCopy.addClone(floatResults);
      const_cast<RooArgList&>(nextResult->floatParsFinal()).remove(floatResults);
      floatCopy.setAttribAll("Constant",kTRUE); //just for pedantic sake
      const_cast<RooArgList&>(nextResult->constPars()).addClone(floatCopy);
      
      //now run the minos 
      std::cout << Form("minos_series: running minos for group %s. Parameters in group: %s",group.Data(),floatCopy.contentsString().c_str()) << std::endl;
      
      TRooFit::minos(nll,pars,nextResult);
      
      
    }
    delete parsInGroup;
    
    out.push_back(nextResult);
    previousResult = nextResult;
  }
  
  //finish by restoring floating status of all parameters
  floatPars->setAttribAll("Constant",kFALSE);
  delete floatPars;
  
  return out;

}

//uses minos_series to produce a set of RooArgList where the asymmetric error on each arg is obtained from a minos series
//The map always contains a 'STAT' entry (or 'STATCORR' and 'STATUNCORR', which is the uncert on the remainder
//You also always get "TOTAL" entry which is all the uncerts together
//if doSTATCORR is true, and more than one pars is given, then STAT will be factorised into a STATCORR and STATUNCORR component
//STATCORR is estimated from pll after fixing ... it's an estimate of the uncertainty due to correlations with other pars
const std::map<TString,RooArgList*> TRooFit::breakdown(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult, bool runInitialMinos, bool doSTATCORR) {

  auto results = minos_series(nll,pars,groups,unconditionalFitResult,runInitialMinos);
  
  //for each fit result, extract the asymm errors and take the quadrature difference to the previous fit result as the error due to that group
  RooArgList allPars(pars);
  
  std::map<TString,RooArgList*> out;
  
  //add the total uncert (which comes from the first result in the vector, which is unconditional fit)
  out["TOTAL"] = new RooArgList;
  RooFIter itr = results[0]->floatParsFinal().fwdIterator();
  while( RooAbsArg* arg = itr.next() ) {
    if(allPars.find(arg->GetName())) out["TOTAL"]->addClone(*arg);
  }

  
  
  for(uint i=1;i<results.size();i++) {
    RooArgList myList; myList.addClone(allPars);
    RooFIter itr = myList.fwdIterator();
    while( RooRealVar* v = (RooRealVar*)itr.next() ) {
      RooRealVar* currV = (RooRealVar*)results[i]->floatParsFinal().find(*v);
      RooRealVar* prevV = (RooRealVar*)results[i-1]->floatParsFinal().find(*v);
      double errHi = 0;
      if(fabs(prevV->getErrorHi()) < fabs(currV->getErrorHi())) {
        std::cout << Form("Warning %s has zero error for group %s",v->GetName(),groups[i-1].Data()) << std::endl;
      } else {
        errHi = sqrt( pow(prevV->getErrorHi(),2) - pow(currV->getErrorHi(),2) );
      }
      double errLo = 0;
      if(fabs(prevV->getErrorLo()) < fabs(currV->getErrorLo())) {
        std::cout << Form("Warning %s has zero error for group %s",v->GetName(),groups[i-1].Data()) << std::endl;
      } else {
        errLo = sqrt( pow(prevV->getErrorLo(),2) - pow(currV->getErrorLo(),2) );
      }
      
      v->setAsymError(-errLo,errHi);
      
    }
    out[groups[i-1]] = new RooArgList;
    out[groups[i-1]]->addClone(myList);
  }
  
  if(!doSTATCORR) {
    //add one more set for the remainder, it's just the uncerts on the last fit result
    RooArgList myList;
    RooFIter itr = results.back()->floatParsFinal().fwdIterator();
    while( RooAbsArg* arg = itr.next() ) {
      if(allPars.find(arg->GetName())) myList.add(*arg);
    }
    out["STAT"] = new RooArgList;
    out["STAT"]->addClone(myList);
  } else {
    //run an additional minos fit for each parameter, fixing all but that parameter to work out the correlated stat uncert ...
    std::cout << Form("breakdown: calculating correlated and uncorrelated remainder uncerts for %d parameters ...",pars.getSize()) << std::endl;
    RooArgList myList;RooArgList myList2;
    RooFIter itr = results.back()->floatParsFinal().fwdIterator();
    while( RooAbsArg* arg = itr.next() ) {
      if(!allPars.find(arg->GetName())) continue;
      
      RooFitResult* decorrFit = new RooFitResult(*results.back());
      //need to move all pars except this par from float to const
      
      RooArgSet floatCopy; floatCopy.addClone(decorrFit->floatParsFinal());
      floatCopy.remove(*arg,true,true);
      floatCopy.setAttribAll("Constant",kTRUE); //just for pedantic sake
      const_cast<RooArgList&>(decorrFit->constPars()).addClone(floatCopy);
      const_cast<RooArgList&>(decorrFit->floatParsFinal()).removeAll();
      const_cast<RooArgList&>(decorrFit->floatParsFinal()).addClone(*arg); //float pars is now the par on its own
      
      
      TRooFit::minos(nll,*pars.find(*arg),decorrFit);
      
      //decorrFit now holds the uncorrelated uncert .. 
      
      RooRealVar* currV = (RooRealVar*)decorrFit->floatParsFinal().find(*arg);
      RooRealVar* prevV = (RooRealVar*)arg;
      double errHi = 0;
      if(fabs(prevV->getErrorHi()) < fabs(currV->getErrorHi())) {
        std::cout << Form("Warning %s has zero error for STATCORR",arg->GetName()) << std::endl;
      } else {
        errHi = sqrt( pow(prevV->getErrorHi(),2) - pow(currV->getErrorHi(),2) );
      }
      double errLo = 0;
      if(fabs(prevV->getErrorLo()) < fabs(currV->getErrorLo())) {
        std::cout << Form("Warning %s has zero error for STATCORR",arg->GetName()) << std::endl;
      } else {
        errLo = sqrt( pow(prevV->getErrorLo(),2) - pow(currV->getErrorLo(),2) );
      }
      
      //copy arg and update the asymm error
      myList.addClone(*arg);
      RooRealVar* v = (RooRealVar*)myList.find(*arg);
      v->setAsymError(-errLo,errHi);
      
      
      myList2.addClone(*currV); //the remaining uncert (just due to that var);
      
      delete decorrFit;
      
    }
    out["STATCORR"] = new RooArgList;out["STATUNCORR"] = new RooArgList;
    out["STATCORR"]->addClone(myList);
    out["STATUNCORR"]->addClone(myList2);
    
    //finish by restoring float status of all parameters
    allPars.setAttribAll("Constant",kFALSE);
    
  }
  
  //finish up by setting name and sorting all lists so that they are consistent
  for(auto& l : out) {
    l.second->setName(l.first);
    l.second->sort();
  }
  
  //delete all fits (except input fit) ...
  for(auto fit : results) {
    if(fit!=unconditionalFitResult) delete fit;
  }
  
  return out;

}



double getSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val, double val_best, double& tmu, bool refit);
void setVal(RooRealVar* par, double val);
double getNLL(RooAbsReal* nll, bool refit);
double getTpar(RooAbsReal* nll, double nll_min, RooRealVar* par, double val, bool refit);

using namespace std;

// ____________________________________________________________________________|__________
// get sigma assuming nll -2logLR is parabolic in par
double TRooFit::findSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val_guess, double val_best, double N_sigma, double precision, bool mode) {
  int printLevel  =   ::ROOT::Math::MinimizerOptions::DefaultPrintLevel();

  TGraph* debugGraph = 0;
  if(GetDebugFile()) {
    TFile* f = GetDebugFile();
    debugGraph = (TGraph*)f->Get(par->GetName());
    if(!debugGraph) {
      debugGraph = new TGraph;
      debugGraph->SetName(par->GetName());
      debugGraph->SetMarkerStyle(5);
      //start by adding the nll_min value for parameter 
      debugGraph->SetPoint(0,par->getVal(),0);
    }
  }


  if (mode == 0)
  {
    bool isConst = par->isConstant();
    par->setConstant(true);
    
    //check if we need to minimize or not
    RooArgSet* pars = nll->getParameters(RooArgSet());
    auto floatPars = pars->selectByAttrib("Constant",kFALSE);
    int nPars = floatPars->getSize();
    delete floatPars; delete pars;
    
    
    
    
    //int direction = int(N_sigma/fabs(N_sigma));
    double tmu;
    //int nrDamping = 1;
    //double damping_factor = 1.0;
    //double damping_factor_pre = damping_factor;
    int nrItr = 0;
    double sigma_guess = fabs((val_guess-val_best)/N_sigma);
    double val_pre = val_guess-10*precision*sigma_guess;
    while (fabs(val_pre-val_guess) > precision*sigma_guess) {
        if(printLevel>1) cout << "----------------------" << endl;
        if(printLevel>1) cout << "Starting iteration " << nrItr << " of " << nll->GetName() << " and parameter " << par->GetName() << endl;
        val_pre = val_guess;
        //damping_factor_pre = damping_factor;
        sigma_guess = getSigma(nll, nll_min, par, val_guess, val_best, tmu,nPars);

        if(debugGraph) {
          //save the pll value (delta(NLL))
          debugGraph->SetPoint(debugGraph->GetN(),val_guess,tmu);
        }

        double corr = /*damping_factor**/(val_pre - val_best - N_sigma*sigma_guess);

        

        // subtract off the difference in the new and damped correction
        val_guess -= corr;

        if(printLevel>1) {
          cout << "nPars:          " << nPars << std::endl;
          cout << "NLL:            " << nll->GetName() << " = " << nll->getVal() << endl;
          cout << "delta(NLL):     " << nll->getVal()-nll_min << endl;
         cout << "NLL min: " << nll_min << std::endl;
          cout << "N_sigma*sigma(pre):   " << fabs(val_pre-val_best) << endl;
          cout << "sigma(guess):   " << sigma_guess << endl;
          cout << "par(guess):     " << val_guess+corr << endl;
          cout << "true val:       " << val_best << endl;
          cout << "tmu:            " << tmu << endl;
          cout << "Precision:      " << sigma_guess*precision << endl;
          cout << "Correction:     " << (-corr<0?" ":"") << -corr << endl;
          cout << "N_sigma*sigma(guess): " << fabs(val_guess-val_best) << endl;
          cout << endl;
        }
        
        

        nrItr++;
        if (nrItr > 25) {
            cout << "Infinite loop detected in getSigma(). Please intervene." << endl;
            break;
        }
    }

    if(printLevel>1) cout << "Found sigma for nll " << nll->GetName() << ": " << (val_guess-val_best)/N_sigma << endl;
    if(printLevel>1) cout << "Finished in " << nrItr << " iterations." << endl;
    if(printLevel>1) cout << endl;
    par->setConstant(isConst);
    
    if(debugGraph) {
      //add a line to indicate the final return value ...
      TLine* l = new TLine(val_guess,0,val_guess,tmu);
      l->SetLineStyle(2);l->SetLineColor(kRed);
      debugGraph->GetListOfFunctions()->Add(l);
      TDirectory* tmpDir= gDirectory;
      GetDebugFile()->cd();
      debugGraph->Write();
      tmpDir->cd();
    }
    
    return (val_guess-val_best)/N_sigma;
  }
  else
  {
    bool isConst = par->isConstant();
    par->setVal(0);
    par->setConstant(1);
    minimize(nll);
    double nll_val = nll->getVal();
    double q0 = 2*(nll_val - nll_min);
    par->setVal(val_best);
    par->setConstant(isConst);

    cout << "q0:       " << q0 << endl;
    cout << "true val: " << val_best << endl;
    cout << "sigma:    " << val_best/sqrt(q0) << endl;

    
    return val_best/sqrt(q0);
  }
}

// ____________________________________________________________________________|__________
double getSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val, double val_best, double& tmu, bool refit=true) {
    tmu = getTpar(nll, nll_min, par, val,refit);
    return fabs(val-val_best)/sqrt(tmu);
}

// ____________________________________________________________________________|__________
double getTpar(RooAbsReal* nll, double nll_min, RooRealVar* par, double val, bool refit=true) {
    setVal(par, val);
    double nll_val = getNLL(nll,refit);
    return 2*(nll_val-nll_min);
}

// ____________________________________________________________________________|__________
double getNLL(RooAbsReal* nll, bool refit=true) {
    
    
    if(!refit) {
      return nll->getVal(); //no fit to do, just return nll;
    } else {
      auto result = TRooFit::minimize(nll,true,false);
      double val = nll->getVal();//result->minNll(); --cannot use minNll because it might have been offset corrected
      delete result;
      return val;
    }
}

// ____________________________________________________________________________|__________
void setVal(RooRealVar* par, double val) {
    if (val > 0 && par->getMax() < val) par->setMax(2*val);
    if (val < 0 && par->getMin() > val) par->setMin(2*val);
    par->setVal(val);
}



//----------------- goodness of fit tests ------------

#include "RooCategory.h"
#include "RooSimultaneous.h"

double TRooFit::GetBakerCousins(RooAbsPdf* model, RooAbsData* data) {
  //log likelihood ratio 'test statistic' as studied by Baker-Cousins
  //The denominator of the ratio is the likelihood for the 'perfect' binned model
  
  //NOTE: 2*LLR is approximately chi2 distribution, ndof ~= nBins (comes out a little higher) -- the constraint terms shift the unsaturated nll (saturated nll doesnt include)
  
  
  //we must start by determing the binning for each observable in each channel of the data
  RooArgSet* obs = model->getObservables(data);
  
  //do we have a category var? If so, assume that is the channel category ...
  RooCategory* cat = 0;
  RooFIter obsItr = obs->fwdIterator();
  while( RooAbsArg* arg = obsItr.next() ) {
    if(arg->IsA() == RooCategory::Class() ) {
      if(cat!=0) {
        msg().Error("Get2LLR","Multiple RooCategory found - unable to infer channel category");
        return 0;
      }
      cat = dynamic_cast<RooCategory*>(arg);
    }
  }
  
  std::vector<RooRealVar*> channelObs; //observable used in each channel
  
  if(cat) {
    //find the RooSimultaneous in the model (might not be the top node) and use each category pdf to try to infer observable
    
    RooSimultaneous* simPdf = 0;
    if(model->InheritsFrom(RooSimultaneous::Class())) {
      simPdf = dynamic_cast<RooSimultaneous*>(model);
    } else {
      //find in the server list ..
      RooArgSet s; model->treeNodeServerList(&s);
      RooFIter itr = s.fwdIterator();
      while( RooAbsArg* arg = itr.next() ) {
        if(arg->InheritsFrom(RooSimultaneous::Class())) {
          if(simPdf!=0) {
            msg().Error("Get2LLR","Multiple RooSimultaneous found - unable to infer main channel pdf");
            return 0;
          }
          simPdf = dynamic_cast<RooSimultaneous*>(arg);
        }
      }
    }
    
    if(!simPdf)  {
      msg().Error("Get2LLR","Unable to infer RooSimultaneous"); 
      return 0;
    }
    
    for(int i=0;i<cat->numTypes();i++) {
      RooAbsPdf* m = simPdf->getPdf( cat->lookupType(i)->GetName() );
      if(!m) {
        msg().Error("Get2LLR","Could not find channel %s in pdf %s",cat->lookupType(i)->GetName(),simPdf->GetName());
        return 0;
      }
      RooArgSet* mobs = m->getObservables(data);
      RooRealVar* var = 0;
      RooFIter obsItr = mobs->fwdIterator();
      while( RooAbsArg* arg = obsItr.next() ) {
        if(arg->IsA() == RooRealVar::Class() ) {
          if(var!=0) {
            msg().Error("Get2LLR","Multiple RooRealVar found - unable to observable");
            return 0;
          }
          var = dynamic_cast<RooRealVar*>(arg);
        }
      }
      delete mobs;
      if(!var) {
        msg().Error("Get2LLR","Could not determine observable for channel %s",cat->lookupType(i)->GetName());
        return 0;
      }
      channelObs.push_back( (RooRealVar*)data->get()->find(*var) ); //save the version of the var present in the dataset
    }
    
    
    
  } else {
    //no channel category
    //there should only be one observable in the dataset
    RooRealVar* var = 0;
    RooFIter obsItr = obs->fwdIterator();
    while( RooAbsArg* arg = obsItr.next() ) {
      if(arg->IsA() == RooRealVar::Class() ) {
        if(var!=0) {
          msg().Error("Get2LLR","Multiple RooRealVar found - unable to observable");
          return 0;
        }
        var = dynamic_cast<RooRealVar*>(arg);
      }
    }
    if(!var) {
      msg().Error("Get2LLR","Could not determine observable");
      return 0;
    }
    channelObs.push_back( (RooRealVar*)data->get()->find(*var) ); //save the version of the var present in the dataset
  }
  
  

  double nll = 0;
  double nll_saturated = 0; //denominator
  
  
  //evaluate the NLL of the data
  
  for(int i=0;i<data->numEntries();i++) {
    *obs = *data->get(i); //move observables values of ith entry
    nll -= model->getLogVal( obs )*data->weight(); //get probability of ith entry (scale by weight)
  }
  //if model is extended, add extended term ..
  if(model->canBeExtended()) nll += model->extendedTerm(data->sumEntries(),obs);
  
  //nll_saturated is the likelihood from the 'saturated' model where each bin prediction equals the data
  //use data histograms to estimate this ...
  //formula is:
  // nll2 = - sum_i [ n_i * log(n_i/(w_i*N)) ] + N - NlogN = -sum_i [ n_i*log(n_i/w_i) ] + N
  //last equality comes from sum_i [n_i * log(N)] = Nlog N 
  
  for(uint i=0;i<channelObs.size();i++) {
    TString channelName = (cat) ? cat->lookupType(i)->GetName() : "";
    for(int j=0;j<channelObs[i]->numBins(channelName);j++) {
      double nObs(0);
      if(cat) nObs = data->sumEntries(Form("%s==%d&&%s>=%f&&%s<%f",cat->GetName(),i,channelObs[i]->GetName(),channelObs[i]->getBinning(channelName).binLow(j),channelObs[i]->GetName(),channelObs[i]->getBinning(channelName).binHigh(j)));
      else nObs =  data->sumEntries(Form("%s>=%f&&%s<%f",channelObs[i]->GetName(),channelObs[i]->getBinning(channelName).binLow(j),channelObs[i]->GetName(),channelObs[i]->getBinning(channelName).binHigh(j)));
      if(nObs<=0) continue;
      nll_saturated -= nObs*log(nObs/channelObs[i]->getBinning(channelName).binWidth(j));
      nll_saturated += nObs; //adds the 'N' term
    }
  }
  
  delete obs;

  
  return 2.*(nll - nll_saturated);
  
}


#include "Math/ProbFuncMathCore.h"

Double_t TRooFit::Phi_m(double mu, double mu_prime, double a, double sigma, int compatCode) {
   //implemented only for special cases of compatibility function
   if(compatCode==0) {
      return 0;
   } 

   //will need sigma unless mu = mu_prime ...
   if(mu!=mu_prime && !sigma) {
    msg().Error("Phi_m","Estimate of mu_hat sigma not provided but required");
   } else if(mu==mu_prime) {
    sigma = 1; //don't need to worry about sigma value when mu= mu_prime
   }
   
   if(compatCode==1) {
      //ignore region below x*sigma+mu_prime = mu ... x = (mu - mu_prime)/sigma 
      if(a < (mu-mu_prime)/sigma) return 0;
      return ROOT::Math::gaussian_cdf(a) - ROOT::Math::gaussian_cdf((mu-mu_prime)/sigma);
   } else if(compatCode==2) {
      //ignore region above x*sigma+mu_prime = mu ... 
      if(a > (mu-mu_prime)/sigma) return ROOT::Math::gaussian_cdf((mu-mu_prime)/sigma);
      return ROOT::Math::gaussian_cdf(a);
   } else if(compatCode==3) {
      //cutoff at x*sigma+mu_prime = mu for mu<0 , and turns back on at x*sigma+mu_prime=mu for mu>0
      //ignore region between x = (-mu-mu_prime)/sigma and x = (mu-mu_prime)/sigma
      double edge1 = (-mu - mu_prime)/sigma;
      double edge2 = (mu-mu_prime)/sigma;

      if(edge1>edge2) { double tmp=edge1; edge1=edge2; edge2=tmp; }

      if(a<edge1) return ROOT::Math::gaussian_cdf(a);
      if(a<edge2) return ROOT::Math::gaussian_cdf(edge1);
      return ROOT::Math::gaussian_cdf(a) - (ROOT::Math::gaussian_cdf(edge2) - ROOT::Math::gaussian_cdf(edge1));
   }
   msg().Error("Phi_m","Unknown compatibility function ... can't evaluate");
   return 0;
}

//must use asimov dataset corresponding to mu_prime, evaluating pll at mu

Double_t TRooFit::asymptoticPValue(double k, RooRealVar* mu, double mu_prime, double sigma, int compatCode) {
  if(k<0) return 1.;
  if(k==0) {
    if(compatCode==2) return 0.5; //when doing discovery (one-sided negative) use a 0.5 pValue
    return 1.; //case to catch the delta function that ends up at exactly 0 for the one-sided tests 
  }
      //get the poi value that defines the test statistic, and the poi_prime hypothesis we are testing 
   //when setting limits, these are often the same value 

   Double_t poiVal = mu->getVal();
   //Double_t k = getVal();

   //if(fOutputLevel<=1) Info("asymptoticPValue","poiVal = %f poiPrimeVal=%f k=%f", poiVal, poi_primeVal, k);
   Double_t Lambda_y = 0;

   double lowBound = mu->getMin();
   double upBound = mu->getMax();

   Lambda_y = (poiVal-mu_prime)/sigma;


   Double_t k_low = (lowBound == -std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((poiVal - lowBound)/sigma,2);
   Double_t k_high = (upBound == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((upBound - poiVal)/sigma,2);

   //std::cout << " sigma = " << sigma << std::endl;

   double out = Phi_m(poiVal,mu_prime,std::numeric_limits<double>::infinity(),sigma,compatCode) - 1;

   //if(fOutputLevel<=0) Info("asymptoticPValue","sigma=%f k_low=%f k_high=%f ", sigma, k_low, k_high);
     
   //std::cout << " out = " << out << std::endl;
   //std::cout << " k_low = " << k_low << std::endl;

   //go through the 4 'regions' ... only two of which will apply
   if( k <= k_low ) {
      out += ROOT::Math::gaussian_cdf(sqrt(k)-Lambda_y) + Phi_m(poiVal,mu_prime,Lambda_y - sqrt(k),sigma,compatCode);
   } else {
      double Lambda_low = (poiVal - lowBound)*(poiVal + lowBound - 2.*mu_prime)/(sigma*sigma);
      double sigma_low = 2.*(poiVal - lowBound)/sigma;
      out +=  ROOT::Math::gaussian_cdf((k-Lambda_low)/sigma_low) + Phi_m(poiVal,mu_prime,(Lambda_low - k)/sigma_low,sigma,compatCode);
   }
   //std::cout << " out = " << out << std::endl;
   //std::cout << " k_high = " << k_high << std::endl;

   if( k <= k_high ) {
      out += ROOT::Math::gaussian_cdf(sqrt(k)+Lambda_y) - Phi_m(poiVal,mu_prime,Lambda_y + sqrt(k),sigma,compatCode);
   } else {
      double Lambda_high = (poiVal - upBound)*(poiVal + upBound - 2.*mu_prime)/(sigma*sigma);
      double sigma_high = 2.*(upBound-poiVal)/sigma;
      out +=  ROOT::Math::gaussian_cdf((k-Lambda_high)/sigma_high) - Phi_m(poiVal,mu_prime,(k - Lambda_high)/sigma_high,sigma,compatCode);
   }


   return 1. - out;
}
