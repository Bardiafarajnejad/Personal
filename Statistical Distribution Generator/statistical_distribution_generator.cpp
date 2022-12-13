//
//  
//  Statistical_distribution_generator
//
//  Created by Bardia Farajnejad on 8/27/21.
//

#include <iostream>
#include <cmath> //for pow(x,y)
#include <iomanip> //for "setprecision(2) << fixed" and for "setw(10)"
#include <cstdlib>  //for srand and rand functions
#include <ctime>  // for time(0): # of seconds since Jan 1970
#include <vector>
#include <utility> // for outputing 2 doubles in "two_normal_BOX_MULLER_RV"
using namespace std;



double one_uniform_RV();
int one_bernouli_RV(double);
int one_binomial_RV(int, double);
int one_poisson_RV(double);

double one_exponential_RV(double);
pair<double, double> two_normal_BOX_MULLER_RV();
pair<double, double> two_normal_POLAR_MARSAGLIA_RV();

int main()
{
    srand(time(0)); //seed random number generator
    
    vector<double> uniforms(10000);
    double total = 0;
    for(int i=1;i<=uniforms.size();i++)
    {
        uniforms[i] = one_uniform_RV();
        total += uniforms[i];
    }
    
    cout << "Avg. of uniform RV's: "<< total/uniforms.size() << endl;
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> bernoulis(10000);
    total = 0;
    for(int i=1;i<=bernoulis.size();i++)
    {
        bernoulis[i] = one_bernouli_RV(0.3);
        total += bernoulis[i];
    }
    
    cout << "Avg. of bernouli RV's: "<< total/bernoulis.size() << endl;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> binomials(10000);
    total = 0;
    for(int i=1;i<=binomials.size();i++)
    {
        binomials[i] = one_binomial_RV(10,0.3);
        total += binomials[i];
    }
    
    cout << "Avg. of binomial RV's: "<< total/binomials.size() << endl;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> poissons(10000);
    total = 0;
    for(int i=1;i<=poissons.size();i++)
    {
        poissons[i] = one_poisson_RV(3);
        total += poissons[i];
    }
    
    cout << "Avg. of poisson RV's: "<< total/poissons.size() << endl;

    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> exponentials(10000);
    total = 0;
    for(int i=1;i<=exponentials.size();i++)
    {
        exponentials[i] = one_exponential_RV(3);
        total += exponentials[i];
    }
    
    cout << "Avg. of exponential RV's: "<< total/exponentials.size() << endl;

    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> normals(10000);
    pair<double,double> result;
    total = 0;
    for(int i=1;i<=normals.size()/2;i++)
    {
        result = two_normal_BOX_MULLER_RV();
        normals[2*i -1] = result.first;
        normals[2*i] = result.second;
        
        total += normals[2*i -1];
        total += normals[2*i];
    }
    
    cout << "Avg. of Normal RV's: "<< total/normals.size() << endl;
    
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<double> normals2(10000);
    pair<double,double> result2;
    total = 0;
    for(int i=1;i<=normals2.size()/2;i++)
    {
        result2 = two_normal_POLAR_MARSAGLIA_RV();
        normals2[2*i -1] = result2.first;
        normals2[2*i] = result2.second;
        
        total += normals2[2*i -1];
        total += normals2[2*i];
    }
    
    cout << "Avg. of Normal RV's: "<< total/normals2.size() << endl;

    return 0;
}



double one_uniform_RV()
{
    int m = pow(2, 15)-1;
    double x_n = rand() % m;
    double unif = (0.5 + x_n)/m;
    return unif;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int one_bernouli_RV(double p)
{
    double U = one_uniform_RV();
    if(U<p)
        return 1;
    else
        return 0;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int one_binomial_RV(int n, double p)
{
    vector<double> bernoulis(n);
    int binomial = 0;
    for(int j=1; j<=n; j++)
    {
        bernoulis[j] = one_bernouli_RV(p);
        binomial += bernoulis[j];
    }
    
    return binomial;
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int one_poisson_RV(double lambda)
{
    double z = exp(-lambda);
    double U = one_uniform_RV();
    int k = 0;
    double x = exp(-lambda);
    
    while(U>=x)
    {
        z *= lambda/(k+1);
        x += z;
        k++;
    }
    return k;
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double one_exponential_RV(double lambda)
{
    double U = one_uniform_RV();
    return (-1/lambda)*log(U);
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pair<double, double> two_normal_BOX_MULLER_RV()
{
    double U1 = one_uniform_RV();
    double U2 = one_uniform_RV();
    double pi = 2*acos(0.0);
    double z1 = sqrt(-2*log(U1))*cos(2*pi*U2);
    double z2 = sqrt(-2*log(U1))*sin(2*pi*U2);
    return make_pair(z1,z2);
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pair<double, double> two_normal_POLAR_MARSAGLIA_RV()
{
    double U1 = one_uniform_RV();
    double U2 = one_uniform_RV();
    double V1 = 2*U1 -1;
    double V2 = 2*U2 -1;
    double W = V1*V1 + V2*V2;
    if(W<=1)
    {
        double z1 = V1*sqrt(-2*log(W)/W);
        double z2 = V2*sqrt(-2*log(W)/W);
        return make_pair(z1,z2);
    }
    else
        return two_normal_POLAR_MARSAGLIA_RV();
};
