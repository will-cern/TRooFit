
#include "TRooFit/TRooChi2Constraint.h"
#include "RooChi2Var.h"


using namespace std;

ClassImp(TRooChi2Constraint)
//ClassImp(TRooGaussian)

TRooChi2Constraint::TRooChi2Constraint(const char *name, const char *title, RooAbsPdf& pdf, RooDataHist& hdata,
	    Bool_t extended, const char* rangeName, const char* addCoefRangeName)
	    : RooAbsPdf(name,title), fChi2("chi2","chi2",this,true,false,false/*owns arg*/)
{
  RooChi2Var* chi2 = new RooChi2Var("chi2","chi2",pdf,hdata,extended,rangeName,addCoefRangeName);
  fChi2.setArg(*chi2); //will take ownership of it, so will manage its deletion
}


//_____________________________________________________________________________
TRooChi2Constraint::TRooChi2Constraint(const TRooChi2Constraint& other, const char* name) : 
  RooAbsPdf(other,name),fChi2("chi2",this,other.fChi2)
{
  // Copy constructor
}



//_____________________________________________________________________________
TRooChi2Constraint::~TRooChi2Constraint()
{
  // Destructor
}

