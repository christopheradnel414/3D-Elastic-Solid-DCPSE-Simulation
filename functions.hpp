#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include <string>

using namespace std;

bool is_number(const std::string& s);

double findmax(vector<double> vect);

double findmin(vector<double> vect);

void calcStressStrain  (int boundaryplotmethod,
                        vector<double>& epsxx,
                        vector<double>& epsyy,
                        vector<double>& epszz,
                        vector<double>& epsxy,
                        vector<double>& epsxz,
                        vector<double>& epsyz,
                        vector<double>& sigmax,
                        vector<double>& sigmay,
                        vector<double>& sigmaz,
                        vector<double>& sigmaxy,
                        vector<double>& sigmaxz,
                        vector<double>& sigmayz,
                        vector<double>& vonmises,
                        vector<double> dUxdX,
                        vector<double> dUxdY,
                        vector<double> dUxdZ,
                        vector<double> dUydX,
                        vector<double> dUydY,
                        vector<double> dUydZ,
                        vector<double> dUzdX,
                        vector<double> dUzdY,
                        vector<double> dUzdZ,
                        double lamda,
                        double Mu,
                        vector<int> dispXi,
                        vector<int> dispYi,
                        vector<int> dispZi,
                        vector<int> forceXi,
                        vector<int> forceYi,
                        vector<int> forceZi);

void writeresult   (int boundaryplotmethod,
                    vector<double> x,
                    vector<double> y,
                    vector<double> z,
                    vector<double> vonmises,
                    vector<double> sigmax,
                    vector<double> sigmay,
                    vector<double> sigmaz,
                    vector<double> sigmaxy,
                    vector<double> sigmaxz,
                    vector<double> sigmayz,
                    vector<double> Ux,
                    vector<double> Uy,
                    vector<double> Uz,
                    double h,
                    vector<int> dispXi,
                    vector<int> dispYi,
                    vector<int> dispZi,
                    vector<int> forceXi,
                    vector<int> forceYi,
                    vector<int> forceZi,
                    string resultfilename);

#endif