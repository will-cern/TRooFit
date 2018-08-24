//static methods that are useful for various activities
//they live in a class so that the root dictionary is generated for them
//so you can use them on the command line

#ifndef TROOFIT
#define TROOFIT

#define protected public
#include "RooFitResult.h"
#undef protected

#include "RooAbsPdf.h"
#include "TRooH1.h"

#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"

#include "TFile.h"

class TRooFit : public TObject  {

  public:
    TRooFit();
    
    ///Static methods of TRooFit, can be used on their own ...
    //returns a significance for an observation given an expectation and a relative uncert (can be asymmetric)
    static double significance(double obs, double exp, double relUncert, double relUncertDown);
    
    
    //set the recommended options for the minimization methods 
    static void setRecommendedDefaultOptions();
    //create a NLL method with the recommended settings
    static RooAbsReal* createNLL(RooAbsPdf* pdf, RooAbsData* data, const RooArgSet* gobs);
    //minimize a function using a good retry strategy ... if save is true then return a RooFitResult
    static RooFitResult* minimize(RooAbsReal* nll, bool save=true, bool hesse=true);
    //updates the asymmetric errors of the given fitResult with the minos errors of the specified pars
    //This method improves on the built-in Minos() method of RooFit, which was found to be unstable
    static RooFitResult* minos(RooAbsReal* nll, const RooArgSet& pars,RooFitResult* unconditionalFitResult = 0, bool respectBoundaries=false);
    //do a series of minos fits, with progressive parameters held constant as determined by the groups
    //this is used when computing a breakdown of uncertainties
    //runInitialMinos should be true if minos has not already been run on the unconditionalFitResult, if it was provided
    static std::vector<RooFitResult*> minos_series(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult=0, bool runInitialMinos=false);
    //do uncertainty breakdown returning lists of the pars with the asymmetric errors set to the error corresponding to each group
    //there is always a 'TOTAL' group for the total asymmetric uncertainty
    static const std::map<TString,RooArgList*> breakdown(RooAbsReal* nll, const RooArgSet& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult=0, bool runInitialMinos=false, bool doSTATCORR=false);
    //same as above method but for a single par, so just return list of up and down values ...
    static const std::map<TString,std::pair<double,double>> breakdown(RooAbsReal* nll, const RooRealVar& pars, std::vector<TString> groups, RooFitResult* unconditionalFitResult=0, bool runInitialMinos=false);
    
    
    static std::pair<RooAbsData*,RooArgSet*> generateAsimovDataset(RooAbsPdf* thePdf, RooAbsData* data, const RooArgSet* gobs=0);
    static std::pair<RooAbsData*,RooArgSet*> generateAsimovDataset(RooAbsPdf* thePdf, const RooArgSet* obs, const RooArgSet* gobs=0);
    static std::pair<RooAbsData*,RooArgSet*> generateToy(RooAbsPdf* model, RooAbsData* data, const RooArgSet* gobs=0, bool doBinned=true);
    static std::pair<RooAbsData*,RooArgSet*> generateToy(RooAbsPdf* model, const RooArgSet* obs, const RooArgSet* gobs=0, bool doBinned=true);
    
    
    //attempts to create a roostats ModelConfig from a workspace containing a model (pdf), data, and give name of parameter of interest
    static RooStats::ModelConfig* CreateModelConfig(RooWorkspace& w, const char* modelName, const char* dataName, const char* poiName );
    
    
    
    //Following is only used when building models with TRooFit classes
    //takes a pdf and dataset and builds a model from it (i.e. the pdf + constraint terms)
    static RooAbsPdf* BuildModel(TRooAbsH1& pdf, RooAbsData& data);
   
    static TObject& msg() { if(m_msgObj==0) m_msgObj=new TObject; return *m_msgObj; }
   
    static void SetColorByName(const char* name, Int_t color) { m_colorsByName[name]=color;}
    static Int_t GetColorByName(const char* name, bool generateColor=false);
    static void PrintColorsByName() { for(auto c : m_colorsByName) { std::cout << c.first << " : " << c.second << std::endl; } }
    
    static void SetDebugFile(TFile* debugFile) { m_debugFile = debugFile; }
    static TFile* GetDebugFile() { return m_debugFile; }
   
   
    static double GetBakerCousins(RooAbsPdf* model, RooAbsData* data);
   
   
   
   
   
    //TRooFit( RooWorkspace& w, const char* poiNames=0 ); 
   
   //function that removes uncertainties on a given histogram
   static void RemoveErrors(TH1* hist, double relUncertThreshold=99999.);
   
   //this is just like RooAbsCollection::setAttribAll but it also dirties all the elements.
   static void setAttribAll(RooAbsCollection& coll, const char* name, Bool_t value);
   
   //some basic compatibility functions 
    static bool TwoSided(double mu, double mu_hat) { return mu==mu_hat; }
    static bool OneSidedPositive(double mu, double mu_hat) { return mu_hat>=mu; }
    static bool OneSidedNegative(double mu, double mu_hat) { return mu_hat<=mu; }
    static bool OneSidedAbsolute(double mu, double mu_hat) { return fabs(mu_hat) >= fabs(mu); }
   
    ClassDef(TRooFit,1);
    
   private:
    static double findSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val_guess, double val_best, double N_sigma, double precision, bool mode = 0);
   
    static TFile* m_debugFile;

    static TObject* m_msgObj;


    static std::map<TString,Int_t> m_colorsByName;
    
  
  
   public:
    class Asymptotics : public TObject {
     public:
      Asymptotics();
    
      
    
      
    
      static Double_t Phi_m(double mu, double mu_prime, double a, double sigma, bool (*_compatibilityFunction)(double mu, double mu_hat));
      static Double_t PValue(double k, double mu, double mu_prime, double sigma_mu, double mu_low, double mu_high, bool (*_compatibilityFunction)(double mu, double mu_hat));
      
      
      static Double_t nullPValue(double k, RooRealVar* mu_test, double sigma_mu, bool (*_compatibilityFunction)(double mu, double mu_hat) ) {
        return PValue(k,mu_test->getVal(),mu_test->getVal(),sigma_mu,mu_test->getMin(),mu_test->getMax(),_compatibilityFunction);
      }
      static Double_t altPValue(double k, RooRealVar* mu_test, double alt_mu, double sigma_mu, bool (*_compatibilityFunction)(double mu, double mu_hat)) {
        return PValue(k,mu_test->getVal(),alt_mu,sigma_mu,mu_test->getMin(),mu_test->getMax(),_compatibilityFunction);
      }
    
    
      static TObject& msg() { return TRooFit::msg(); }
      
      ClassDef(TRooFit::Asymptotics,1);
  
  };
  
};

#endif