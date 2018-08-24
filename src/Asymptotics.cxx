

#include "TRooFit/Utils.h"
#include "Math/ProbFuncMathCore.h"

TRooFit::Asymptotics::Asymptotics() : TObject() { }

Double_t TRooFit::Asymptotics::Phi_m(double mu, double mu_prime, double a, double sigma, bool (*compatFunc)(double mu, double mu_hat) ) {
   //implemented only for special cases of compatibility function
   if(compatFunc==TwoSided) {
      return 0;
   } 

   if(compatFunc==OneSidedPositive) {
      //ignore region below x*sigma+mu_prime = mu ... x = (mu - mu_prime)/sigma 
      if(a < (mu-mu_prime)/sigma) return 0;
      return ROOT::Math::gaussian_cdf(a) - ROOT::Math::gaussian_cdf((mu-mu_prime)/sigma);
   } else if(compatFunc==OneSidedNegative) {
      //ignore region above x*sigma+mu_prime = mu ... 
      if(a > (mu-mu_prime)/sigma) return ROOT::Math::gaussian_cdf((mu-mu_prime)/sigma);
      return ROOT::Math::gaussian_cdf(a);
   } else if(compatFunc==OneSidedAbsolute) {
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
Double_t TRooFit::Asymptotics::PValue(double k, double poiVal, double poi_primeVal, double sigma, double lowBound, double upBound, bool (*compatFunc)(double mu, double mu_hat)) {
  if(k<0) return 1.;
  if(k==0) {
    if(compatFunc==OneSidedNegative) return 0.5; //when doing discovery (one-sided negative) use a 0.5 pValue
    return 1.; //case to catch the delta function that ends up at exactly 0 for the one-sided tests 
  }
  
  if(sigma==0 && (lowBound!=-std::numeric_limits<double>::infinity() || upBound != std::numeric_limits<double>::infinity())) {
      msg().Error("TRooFit::Asymptotic::PValue","You must specify a sigma_mu parameter when parameter of interest is bounded");
      return -1;
   }
  
      //get the poi value that defines the test statistic, and the poi_prime hypothesis we are testing 
   //when setting limits, these are often the same value 



   Double_t Lambda_y = 0;
   if(fabs(poiVal-poi_primeVal)>1e-9) Lambda_y = (poiVal-poi_primeVal)/sigma;

   

   Double_t k_low = (lowBound == -std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((poiVal - lowBound)/sigma,2);
   Double_t k_high = (upBound == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((upBound - poiVal)/sigma,2);

   //std::cout << " sigma = " << sigma << std::endl;

   double out = Phi_m(poiVal,poi_primeVal,std::numeric_limits<double>::infinity(),sigma,compatFunc) - 1;

   //if(fOutputLevel<=0) msg().Info("PValue","sigma=%f k_low=%f k_high=%f ", sigma, k_low, k_high);
     
   //std::cout << " out = " << out << std::endl;
   //std::cout << " k_low = " << k_low << std::endl;

   //go through the 4 'regions' ... only two of which will apply
   if( k <= k_low ) {
      out += ROOT::Math::gaussian_cdf(sqrt(k)-Lambda_y) + Phi_m(poiVal,poi_primeVal,Lambda_y - sqrt(k),sigma,compatFunc);
   } else {
      double Lambda_low = (poiVal - lowBound)*(poiVal + lowBound - 2.*poi_primeVal)/(sigma*sigma);
      double sigma_low = 2.*(poiVal - lowBound)/sigma;
      out +=  ROOT::Math::gaussian_cdf((k-Lambda_low)/sigma_low) + Phi_m(poiVal,poi_primeVal,(Lambda_low - k)/sigma_low,sigma,compatFunc);
   }
   //std::cout << " out = " << out << std::endl;
   //std::cout << " k_high = " << k_high << std::endl;

   if( k <= k_high ) {
      out += ROOT::Math::gaussian_cdf(sqrt(k)+Lambda_y) - Phi_m(poiVal,poi_primeVal,Lambda_y + sqrt(k),sigma,compatFunc);
   } else {
      double Lambda_high = (poiVal - upBound)*(poiVal + upBound - 2.*poi_primeVal)/(sigma*sigma);
      double sigma_high = 2.*(upBound-poiVal)/sigma;
      out +=  ROOT::Math::gaussian_cdf((k-Lambda_high)/sigma_high) - Phi_m(poiVal,poi_primeVal,(k - Lambda_high)/sigma_high,sigma,compatFunc);
   }


   return 1. - out;
}
