#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#define _USE_MATH_DEFINES
#include <vector>
#include <string>
# define M_PI           3.14159265358979323846 
#define hbarc  0.19733  //GeV*fm
#define str to_string
using namespace std;


//hypersurface storage
class Hypersurface{
    public:
    double  t, z, u_t, u_z, dsigma_t, dsigma_z,
            s, T, p , e;

    Hypersurface(double x0, double x3, double u_0, double u_3, 
                 double dsigma_0, double dsigma_3,
                 double s_, double T_, double p_ , double e_){

            t = x0*cosh(x3);
            z = x0*sinh(x3);
            u_t = cosh(x3)*u_0 - sinh(x3)*u_3/x0;
            u_z = -sinh(x3)*u_0 + cosh(x3)*u_3/x0;
            dsigma_t = cosh(x3)*dsigma_0 - sinh(x3)*dsigma_3/x0;
            dsigma_z = -sinh(x3)*dsigma_0 + cosh(x3)*dsigma_3/x0;
            s = s_;
            T = T_/hbarc;
            p = p_;
            e = e_;
    };
};

//Cooper-Frye Calculations 
class CooperFrye{
    public:
    vector<Hypersurface> hs;

    CooperFrye(string path){
        ifstream infile(path);
        string line;
        
        while (getline(infile, line)){
            double  x0, x3, u_0, u_3, dsigma_0, dsigma_3,
                    s, T, p , e;
            try{
                std::istringstream iss(line);
                if (iss >> x0 >> x3 >> u_0 >> u_3 >> dsigma_0
                        >> dsigma_3 >> s >> T >> p >> e) {

                    Hypersurface aux(x0, x3, u_0, u_3, 
                                    dsigma_0, dsigma_3,
                                    s, T, p, e);

                    hs.push_back(aux);
                }
            }catch (const std::exception& e) { continue; } 

        }
    }

    //Phase-Space Distribution
    double phaseSpaceDist(double x, double g, double a){
        return (g/pow(2*M_PI,3.))*(1./(exp(x) + a));
    }

    //Ed3N/d3p - Invariant One-Particle Distribution
    double onePartDist(double *p, double g, double a){
        double ed3nd3p = 0.;

        for(int i = 0; i < hs.size(); i++){
            double  pcdotu = p[0]*hs[i].u_t + p[3]*hs[i].u_z,
                    pcdotdsigma = p[0]*hs[i].dsigma_t + p[3]*hs[i].dsigma_z,
                    f = phaseSpaceDist(pcdotu/hs[i].T, g, a);
                    //cout<<"T = "<<hs[i].T<<endl;
            if(pcdotdsigma > 0) ed3nd3p += pcdotdsigma*f;
        }

        return ed3nd3p;
    };

    //dN/deta - Pseudorapidity Distribution
    double psRapidityDist(double eta, double g, double a, double m){
        double dndeta = 0.;
        double  mt0 = sqrt(m*m+pow(0.15/hbarc,2.)),
                mtf  = sqrt(m*m + pow(5.5/hbarc,2.)),
                n = 10000.,
                dmt = (mtf - mt0)/(n);

        double  mt = mt0,
                pt = sqrt(mt*mt-m*m),
                y = (1./2.)*log((sqrt(pow(pt*cosh(eta),2.) +m*m) + pt*sinh(eta))
                                /(sqrt(pow(pt*cosh(eta),2.) +m*m) - pt*sinh(eta))),
                p[4] = {mt*cosh(y), pt, pt, mt*sinh(y)},
                f1 = mt*onePartDist(p, g, a)*sqrt(pt*pt + p[3]*p[3])/p[0],
                f2 = 0.;
                //cout<<"pt = "<<pt<<endl;
                

        for(int i = 0; i<int(n); i++){
            if(i > 0)   f1 = f2;
            mt = mt0 + (i+1)*dmt; 
            pt = sqrt(mt*mt-m*m);
            y = (1./2.)*log((sqrt(pow(pt*cosh(eta),2.) +m*m) + pt*sinh(eta))
                            /(sqrt(pow(pt*cosh(eta),2.) +m*m) - pt*sinh(eta)));
            double p2[4] = {mt*cosh(y), pt, pt, mt*sinh(y)};
            f2 = mt*onePartDist(p2, g, a)*sqrt(pt*pt + p2[3]*p2[3])/p2[0];

            dndeta += (f1 + f2)*dmt/2;
        }

        return 2*M_PI*dndeta;
    }

};

// Smearing Kernel
double gaussian(double x, double sigma){
    return 1./sqrt(2*M_PI)*sigma*exp(-(1./2.)*pow(x/sigma,2.));
}

// Smearing Function
double smearing(double x_, double dx, double *x, double *f, int a){
    double f_ = 0., sigma = 1.,  part1 = 0., part2 = 0.;

    for(int i = 0; i<=a; i++){
        if(i > 0) part1 = part2;
        else part1 = gaussian(x[i] - x_, sigma)*f[i];
        part2 = gaussian(x[i+1] - x_, sigma)*f[i+1];
        f_ += (part1+part2)*dx/2.;
    }

    return f_;
}


int main(int argc, char** argv){

    CooperFrye cp(argv[1]);

    //proton+ 0.938 2 1
    //kaon+ 0.494 1 -1
    //pion+ 0.140 1 -1

    double m = 0.140, //0.494
           g = 1.,
           a = -1.;

    double  eta0 = -6., etaf = 6., n = 350., deta = (etaf - eta0)/n;

    double eta[int(n)+1];
    double dndeta[int(n) + 1];

    ofstream outfile("dndeta.dat");
    for(int i = 0; i <= int(n); i++){
        eta[i] = eta0 + i*deta;
        dndeta[i] =   cp.psRapidityDist(eta[i], 1., -1., 0.140/hbarc) 
                    + cp.psRapidityDist(eta[i], 1., -1., 0.494/hbarc) 
                    + cp.psRapidityDist(eta[i], 2., 1., 0.938/hbarc);
        outfile<<eta[i]<<" "<<dndeta[i]<<endl;
    }
    outfile.close();

    ofstream outfile2("sdndeta.dat");
    for(int i = 0; i <= int(n); i++){
        outfile2<<eta[i]<<" "<<smearing(eta[i], deta, eta, dndeta, int(n))<<endl;
    }
    outfile2.close();





    return 0;
}
