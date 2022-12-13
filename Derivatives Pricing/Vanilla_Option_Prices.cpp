#include <cmath>             // standard mathematical library
#include <algorithm>             // defining the max() operator
#include <vector>           // STL vector templates
#include <iostream>
#include <random>
using namespace std;

//2 old functions from Prof. G
double option_price_put_american_binomial( const double&, const double&, const double&, const double&, const double&, const int&);

double option_price_put_european_binomial( const double&, const double&, const double&, const double&, const double&, const int&);

//4 new functions:

// USING INLINE TO MAKE A COPY OF THIS FUNCTION AND MAKE CALLING IT 10,000 times FASTER
inline double two_antithetic_paths( const double&, const double&, const double&, const double&, const double&);

double option_price_put_european_Monte_Carlo( const double&, const double&, const double&, const double&, const double&, const int&);

double BS_put( const double&, const double&, const double&, const double&, const double&);

double cfd_normal_at_x(double);





int main()
{
    cout<<option_price_put_european_binomial( 80, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_european_binomial( 90, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_european_binomial( 100, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_european_binomial( 110, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_european_binomial( 120, 100, .05, 0.30, 1, 12)<<endl;

    cout<<option_price_put_european_Monte_Carlo( 80, 100, .05, 0.30, 1, 40000000)<<endl; //using inline function inside "option_price_put_european_Monte_Carlo"
    cout<<option_price_put_european_Monte_Carlo( 90, 100, .05, 0.30, 1, 40000000)<<endl; //using inline function inside "option_price_put_european_Monte_Carlo"
    cout<<option_price_put_european_Monte_Carlo( 100, 100, .05, 0.30, 1, 40000000)<<endl; //using inline function inside "option_price_put_european_Monte_Carlo"
    cout<<option_price_put_european_Monte_Carlo( 110, 100, .05, 0.30, 1, 40000000)<<endl; //using inline function inside "option_price_put_european_Monte_Carlo"
    cout<<option_price_put_european_Monte_Carlo( 120, 100, .05, 0.30, 1, 40000000)<<endl; //using inline function inside "option_price_put_european_Monte_Carlo"

    cout<<BS_put( 80, 100, .05, 0.30, 1)<<endl;
    cout<<BS_put( 90, 100, .05, 0.30, 1)<<endl;
    cout<<BS_put( 100, 100, .05, 0.30, 1)<<endl;
    cout<<BS_put( 110, 100, .05, 0.30, 1)<<endl;
    cout<<BS_put( 120, 100, .05, 0.30, 1)<<endl;

    
    cout<<option_price_put_american_binomial( 80, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_american_binomial( 90, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_american_binomial( 100, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_american_binomial( 110, 100, .05, 0.30, 1, 12)<<endl;
    cout<<option_price_put_american_binomial( 120, 100, .05, 0.30, 1, 12)<<endl;
}









double option_price_put_american_binomial( const double& S,
                        const double& K,
                        const double& r,
                        const double& sigma,
                        const double& t,
                        const int& steps) {
   double R = exp(r*(t/steps));
   double Rinv = 1.0/R; //... this is exp(-r*(t/steps));
   double u = exp(sigma*sqrt(t/steps));
   double d = 1.0/u;
   double p_up = (R-d)/(u-d);
   double p_down = 1.0-p_up;

   vector<double> prices(steps+1);       // price of underlying
   prices[0] = S*pow(d, steps);  // fill in the endnodes.
   double uu = u*u;
   for (int i=1; i<=steps; ++i) prices[i] = uu*prices[i-1]; // fill in the endnodes.

   vector<double> call_values(steps+1);       // value of corresponding call
   for (int i=0; i<=steps; ++i) call_values[i] = max(0.0, (K-prices[i])); // call payoffs at maturity   : CHANGE TO K-S

   for (int step=steps-1; step>=0; --step) {
      for (int i=0; i<=step; ++i) {
     call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
     prices[i] = d*prices[i+1]; //CAN ALSO DO "prices[i] = u*prices[i];" here !
     call_values[i] = max(call_values[i],K-prices[i]);       // check for exercise   : CHANGE TO K-S
      };
   };
   return call_values[0];
};





double option_price_put_european_binomial( const double& S,     // spot price
                        const double& K,     // exercice price
                        const double& r,     // interest rate
                        const double& sigma, // volatility
                        const double& t,     // time to maturity
                        const int& steps){  // no steps in binomial tree
   double R = exp(r*(t/steps));            // interest rate for each step
   double Rinv = 1.0/R;                    // inverse of interest rate
   double u = exp(sigma*sqrt(t/steps));    // up movement
   double uu = u*u;
   double d = 1.0/u;
   double p_up = (R-d)/(u-d);
   double p_down = 1.0-p_up;
   vector<double> prices(steps+1);       // price of underlying
   prices[0] = S*pow(d, steps);  // fill in the endnodes.
   for (int i=1; i<=steps; ++i) prices[i] = uu*prices[i-1]; // fill in the endnodes.
    
   vector<double> call_values(steps+1);       // value of corresponding call
   for (int i=0; i<=steps; ++i) call_values[i] = max(0.0, (K-prices[i])); // call payoffs at maturity  : CHANGE TO K-S
    
   for (int step=steps-1; step>=0; --step) {
      for (int i=0; i<=step; ++i) {
     call_values[i] = (p_up*call_values[i+1]+p_down*call_values[i])*Rinv;
      };
   };
   return call_values[0];
};








default_random_engine generator;
normal_distribution<double> distribution(0.0,1.0);

// USING INLINE TO MAKE A COPY OF THIS FUNCTION AND MAKE CALLING IT 10,000 times FASTER!!
inline double two_antithetic_paths( const double& S,     // spot price
                        const double& K,     // exercice price
                        const double& r,     // interest rate
                        const double& sigma, // volatility
                        const double& t){   // time to maturity
    //This function will output the average of two antithetic variates (a variance reduction technique) to simulate the payoff of // a European put
    
    double norm;
    double W;
    //double W2;
    double S_t;
    double S2_t;
    double payoff;
    double payoff2;
    double price;
    double price2;

    norm = distribution(generator);
    W = sqrt(t)*norm;
    //W2 =pow(t, 0.5)*-norm;
        
    S_t = S * exp(sigma*W + (r-(pow(sigma,2)/2))*t);
    //S2_t = S * exp(sigma*W2 + (r-(pow(sigma,2)/2))*t);
    S2_t = S * exp(sigma*-W + (r-(pow(sigma,2)/2))*t);

    if (K - S_t > 0)
        payoff = K - S_t;
    else
        payoff = 0;
        
    if (K - S2_t > 0)
        payoff2 = K - S2_t;
    else
        payoff2 = 0;

    price = exp(-r*t)*payoff;
    price2 = exp(-r*t)*payoff2;
        
    return (price+price2)/2;  //return the average of the two paths!
};


double option_price_put_european_Monte_Carlo( const double& S,     // spot price
                        const double& K,     // exercice price
                        const double& r,     // interest rate
                        const double& sigma, // volatility
                        const double& t,     // time to maturity
                        const int& sims){  // no sims in MC sim.... PUT AN EVEN INTEGER HERE PLEASE !!
    //This function will use antithetic variates (a variance reduction technique) to simulate the payoff of // a European put
    
    vector<double> avg_price(sims/2);
    double sum_of_avg_prices = 0; //counter
    
    for(int i=1; i<=sims/2;i++){
        avg_price[i] = two_antithetic_paths(S,K,r,sigma,t);
        sum_of_avg_prices += avg_price[i];
    }
    
    return sum_of_avg_prices/avg_price.size();
};








double cfd_normal_at_x(double x){
    double d1 = 0.0498673470;
    double d2 = 0.0211410061;
    double d3 = 0.0032776263;
    double d4 = 0.0000380036;
    double d5 = 0.0000488906;
    double d6 = 0.0000053830;
    
    double x2 = pow(x,2);
    double x3 = pow(x,3);
    double x4 = pow(x,4);
    double x5 = pow(x,5);
    double x6 = pow(x,6);
    if(x>=0)
        return( 1-0.5*pow(1+d1*x+d2*x2+d3*x3+d4*x4+d5*x5+d6*x6,-16) );
    else
        return(1-cfd_normal_at_x(-x));
};




double BS_put(const double& S,     // spot price
               const double& K,     // exercice price
               const double& r,     // interest rate
               const double& sigma, // volatility
               const double& t){     // time to maturity

    double d1 = ( log(S/K) + (r+pow(sigma,2)/2)*t )/(sigma*sqrt(t));
    double d2 = d1-sigma*sqrt(t);
    
    return( -S*cfd_normal_at_x(-d1)+K*exp(-r*t)*cfd_normal_at_x(-d2) );
};
