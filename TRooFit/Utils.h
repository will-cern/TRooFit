//static methods that are useful for various activities
//they live in a class so that the root dictionary is generated for them
//so you can use them on the command line

#ifndef TROOFIT
#define TROOFIT

#include "RooAbsPdf.h"
#include "TRooH1.h"

#include "RooStats/ModelConfig.h"

#include "TFile.h"

class TRooFit : public TObject  {

  public:
    TRooFit();
    
    //takes a pdf and dataset and builds a model from it (i.e. the pdf + constraint terms)
    static RooAbsPdf* BuildModel(TRooAbsH1& pdf, RooAbsData& data);

    //creates a roostats ModelConfig from a workspace containing a model (pdf), data, and give name of parameter of interest
    static RooStats::ModelConfig* CreateModelConfig(RooWorkspace& w, const char* modelName, const char* dataName, const char* poiName );
   
   
    //minimize a function using a good retry strategy
    static RooFitResult* minimize(RooAbsReal* nll, bool save=true);
   
    //updates the asymmetric errors of the given fitResult with the minos errors of the specified pars
    static RooFitResult* minos(RooAbsReal* nll, const RooArgSet& pars,RooFitResult* unconditionalFitResult = 0);
   
    //do a series of minos fits, with progress parameters held constant as determined by the groups
    static std::vector<RooFitResult*> minos_series(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult=0, bool runInitialMinos=false);
   
    //do uncertainty breakdown
    //caller owns lists
    static const std::map<TString,RooArgList*> breakdown(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult=0, bool runInitialMinos=false, bool doSTATCORR=false);
    
    static double findSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val_guess, double val_best, double N_sigma, double precision, bool mode = 0);
   
   
    
   
    static void setRecommendedDefaultOptions();
    
    static void SetDebugFile(TFile* debugFile) { m_debugFile = debugFile; }
    static TFile* GetDebugFile() { return m_debugFile; }
   
    ClassDef(TRooFit,1);
    
   private:
    static TFile* m_debugFile;

};

#endif