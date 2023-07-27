#include <sstream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iostream>
#include "functions.hpp"

using namespace std;

bool is_number(const std::string& s)
{
    long double ld;
    return((std::istringstream(s) >> ld >> std::ws).eof());
}

double findmax(vector<double> vect)
{
    int size = vect.size();
    double maxval = vect[0];
    for (int i = 0; i < size; i++)
    {
        if (vect[i] > maxval)
        {
            maxval = vect[i];
        }
    }
    return maxval;
}

double findmin(vector<double> vect)
{
    int size = vect.size();
    double minval = vect[0];
    for (int i = 0; i < size; i++)
    {
        if (vect[i] < minval)
        {
            minval = vect[i];
        }
    }
    return minval;
}

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
                        vector<int> forceZi)
{
    cout << "Calculating Stress and Strain ... ";

    auto start =chrono::steady_clock::now();
    
    int particlenumber = dUxdX.size();
    int j = 0;
    if (boundaryplotmethod == 1)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            epsxx.push_back(dUxdX[i]);
            epsyy.push_back(dUydY[i]);
            epszz.push_back(dUzdZ[i]);
            epsxy.push_back(0.5*(dUxdY[i] + dUydX[i]));
            epsxz.push_back(0.5*(dUxdZ[i] + dUzdX[i]));
            epsyz.push_back(0.5*(dUydZ[i] + dUzdY[i]));
            
            sigmax.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsxx[j]));
            sigmay.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsyy[j]));
            sigmaz.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epszz[j]));
            sigmaxy.push_back(2*Mu*epsxy[j]);
            sigmaxz.push_back(2*Mu*epsxz[j]);
            sigmayz.push_back(2*Mu*epsyz[j]);
            vonmises.push_back(sqrt(pow(sigmax[j],2) + pow(sigmay[j],2) + pow(sigmaz[j],2) - sigmax[j]*sigmay[j] - sigmax[j]*sigmaz[j] - sigmay[j]*sigmaz[j] + 3*(pow(sigmaxy[j],2) + pow(sigmaxz[j],2) + pow(sigmayz[j],2))));
            j = j + 1;
        }
    }
    else if (boundaryplotmethod == 2)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            if ((dispXi[i] == 0) && ((dispYi[i] == 0) && (dispZi[i] == 0)))
            {
                epsxx.push_back(dUxdX[i]);
                epsyy.push_back(dUydY[i]);
                epszz.push_back(dUzdZ[i]);
                epsxy.push_back(0.5*(dUxdY[i] + dUydX[i]));
                epsxz.push_back(0.5*(dUxdZ[i] + dUzdX[i]));
                epsyz.push_back(0.5*(dUydZ[i] + dUzdY[i]));
                
                sigmax.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsxx[j]));
                sigmay.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsyy[j]));
                sigmaz.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epszz[j]));
                sigmaxy.push_back(2*Mu*epsxy[j]);
                sigmaxz.push_back(2*Mu*epsxz[j]);
                sigmayz.push_back(2*Mu*epsyz[j]);
                vonmises.push_back(sqrt(pow(sigmax[j],2) + pow(sigmay[j],2) + pow(sigmaz[j],2) - sigmax[j]*sigmay[j] - sigmax[j]*sigmaz[j] - sigmay[j]*sigmaz[j] + 3*(pow(sigmaxy[j],2) + pow(sigmaxz[j],2) + pow(sigmayz[j],2))));
                j = j + 1;
            }
        }
    }
    else if (boundaryplotmethod == 3)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            if ((forceXi[i] == 0) && ((forceYi[i] == 0) && ((forceZi[i] == 0) && ((dispXi[i] == 0) && ((dispYi[i] == 0) && (dispZi[i] == 0))))))
            {
                epsxx.push_back(dUxdX[i]);
                epsyy.push_back(dUydY[i]);
                epszz.push_back(dUzdZ[i]);
                epsxy.push_back(0.5*(dUxdY[i] + dUydX[i]));
                epsxz.push_back(0.5*(dUxdZ[i] + dUzdX[i]));
                epsyz.push_back(0.5*(dUydZ[i] + dUzdY[i]));
                
                sigmax.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsxx[j]));
                sigmay.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsyy[j]));
                sigmaz.push_back(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epszz[j]));
                sigmaxy.push_back(2*Mu*epsxy[j]);
                sigmaxz.push_back(2*Mu*epsxz[j]);
                sigmayz.push_back(2*Mu*epsyz[j]);
                vonmises.push_back(sqrt(pow(sigmax[j],2) + pow(sigmay[j],2) + pow(sigmaz[j],2) - sigmax[j]*sigmay[j] - sigmax[j]*sigmaz[j] - sigmay[j]*sigmaz[j] + 3*(pow(sigmaxy[j],2) + pow(sigmaxz[j],2) + pow(sigmayz[j],2))));
                j = j + 1;
            }
        }
    }
    else
    {
        cout << "Boundary plot method selection INVALID" <<endl;
        exit(EXIT_FAILURE);
    }
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

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
                    string resultfilename)
{
    cout << "Writing Result File ... ";

    auto start =chrono::steady_clock::now();

    ofstream file;
    file.open(resultfilename);
    file<<h<<endl;

    int particlenumber = x.size();
    int j = 0;
    if (boundaryplotmethod == 1)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            file<<x[i]<<"\t";
            file<<y[i]<<"\t";
            file<<z[i]<<"\t";
            file<<vonmises[j]<<"\t";
            file<<sigmax[j]<<"\t";
            file<<sigmay[j]<<"\t";
            file<<sigmaz[j]<<"\t";
            file<<sigmaxy[j]<<"\t";
            file<<sigmaxz[j]<<"\t";
            file<<sigmayz[j]<<"\t";
            file<<Ux[i]<<"\t";
            file<<Uy[i]<<"\t";
            file<<Uz[i]<<endl;
            j++;
        }
    }
    else if (boundaryplotmethod == 2)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            if ((dispXi[i] == 0) && ((dispYi[i] == 0) && (dispZi[i] == 0)))
            {
                file<<x[i]<<"\t";
                file<<y[i]<<"\t";
                file<<z[i]<<"\t";
                file<<vonmises[j]<<"\t";
                file<<sigmax[j]<<"\t";
                file<<sigmay[j]<<"\t";
                file<<sigmaz[j]<<"\t";
                file<<sigmaxy[j]<<"\t";
                file<<sigmaxz[j]<<"\t";
                file<<sigmayz[j]<<"\t";
                file<<Ux[i]<<"\t";
                file<<Uy[i]<<"\t";
                file<<Uz[i]<<endl;
                j++;
            }
            else
            {
                file<<"nan";
            }
            file<<endl;
        }
    }
    else if (boundaryplotmethod == 3)
    {
        for (int i = 0; i < particlenumber; i++)
        {
            if ((forceXi[i] == 0) && ((forceYi[i] == 0) && ((forceZi[i] == 0) && ((dispXi[i] == 0) && ((dispYi[i] == 0) && (dispZi[i] == 0))))))
            {
                file<<x[i]<<"\t";
                file<<y[i]<<"\t";
                file<<z[i]<<"\t";
                file<<vonmises[j]<<"\t";
                file<<sigmax[j]<<"\t";
                file<<sigmay[j]<<"\t";
                file<<sigmaz[j]<<"\t";
                file<<sigmaxy[j]<<"\t";
                file<<sigmaxz[j]<<"\t";
                file<<sigmayz[j]<<"\t";
                file<<Ux[i]<<"\t";
                file<<Uy[i]<<"\t";
                file<<Uz[i]<<endl;
                j++;
            }
            else
            {
                file<<"nan";
            }
            file<<endl;
        }
    }
    file.close();
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}