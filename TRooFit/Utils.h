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
   
   
   static Double_t asymptoticPValue(double k, RooRealVar* mu, double mu_prime, double sigma, int compatCode);
   
   
    //TRooFit( RooWorkspace& w, const char* poiNames=0 ); 
   
   //function that removes uncertainties on a given histogram
   static void RemoveErrors(TH1* hist, double relUncertThreshold=99999.);
   
   
    ClassDef(TRooFit,1);
    
   private:
    static double findSigma(RooAbsReal* nll, double nll_min, RooRealVar* par, double val_guess, double val_best, double N_sigma, double precision, bool mode = 0);
   
    static TFile* m_debugFile;

    static TObject* m_msgObj;


    static Double_t Phi_m(double mu, double mu_prime, double a, double sigma, int compatCode);

    //RooAbsReal* m_nll = 0;
    
    static std::map<TString,Int_t> m_colorsByName;
    
  
};

#endif