
#include "TRooFit/Utils.h"

ClassImp(TRooFit)

TRooFit::TRooFit() : TObject() { }

RooAbsPdf* TRooFit::BuildModel(TRooAbsH1& pdf, RooAbsData& data) {
  return pdf.buildConstraints( *data.get() , "", true );
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