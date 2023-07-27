#include <chrono>
#include <string>
#include <iostream>
#include "functions.hpp"
#include "readgeom.hpp"
#include "neighbourfind.hpp"
#include "DCPSEgeneral.hpp"
#include "MatrixOperation.hpp"
#include "readparameter.hpp"
#include "eigen-3.3.7/Eigen/Dense"

using namespace std;

int main()
{
    auto start = chrono::steady_clock::now();

    // Reading Solver Parameters Data
    int boundaryplotmethod;
    string geomfilename;
    string geomparamfilename;
    string resultfilename;
    int neighbourmethod;
    int neighbournumber;
    double rcScale;
    int method_solver;
    int preconditioner_iterative_solver;
    double tolerance_iterative_solver;
    int iteration_limit_iterative_solver;
    int n_cpu_iterative_solver;

    readsolverparameterfile(boundaryplotmethod,geomfilename,geomparamfilename,
                            resultfilename,neighbourmethod,neighbournumber,
                            rcScale,method_solver,preconditioner_iterative_solver,
                            tolerance_iterative_solver,iteration_limit_iterative_solver,
                            n_cpu_iterative_solver);

    // Initializing the geometry parameters
    double E;
    double v;
    double h;

    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> nx;
    vector<double> ny;
    vector<double> nz;
    vector<int> nxi;
    vector<int> nyi;
    vector<int> nzi;
    vector<double> dispX;
    vector<double> dispY;
    vector<double> dispZ;
    vector<int> dispXi;
    vector<int> dispYi;
    vector<int> dispZi;
    vector<double> forceX;
    vector<double> forceY;
    vector<double> forceZ;
    vector<int> forceXi;
    vector<int> forceYi;
    vector<int> forceZi;
    vector<double> bx;
    vector<double> by;
    vector<double> bz;

    // Reading the geometry file
    readgeom   (E,v,h,x,y,z,nx,ny,nz,nxi,nyi,nzi,dispX,dispY,dispZ,dispXi,dispYi,
                dispZi,forceX,forceY,forceZ,forceXi,forceYi,forceZi,bx,by,bz,
                geomfilename,geomparamfilename);

    int particlenumber = x.size();

    // Nearest Neighbour Search
    double rc = rcScale*h;
    vector<vector<int>> neighbour(particlenumber);
    if (neighbourmethod == 1)
    {
        neighbourfind(neighbour,x,y,z,rc);
    }
    else if (neighbourmethod == 2)
    {
        neighbourfindlimited(neighbour,x,y,z,rc,neighbournumber);
    }
    for (int i = 0; i < particlenumber; i++)
    {
        neighbour[i].push_back(i);
    }

    // DCPSE Operator
    vector<MatrixXd> Ainv = DCPSEcalcAinv(x,y,z,neighbour,h);
    
    vector<vector<double>> EtaDx = DCPSEcalcEta(x,y,z,neighbour,h,1,0,0,Ainv);
    vector<vector<double>> EtaDy = DCPSEcalcEta(x,y,z,neighbour,h,0,1,0,Ainv);
    vector<vector<double>> EtaDz = DCPSEcalcEta(x,y,z,neighbour,h,0,0,1,Ainv);
    vector<vector<double>> EtaDx2 = DCPSEcalcEta(x,y,z,neighbour,h,2,0,0,Ainv);
    vector<vector<double>> EtaDy2 = DCPSEcalcEta(x,y,z,neighbour,h,0,2,0,Ainv);
    vector<vector<double>> EtaDz2 = DCPSEcalcEta(x,y,z,neighbour,h,0,0,2,Ainv);
    vector<vector<double>> EtaDxDy = DCPSEcalcEta(x,y,z,neighbour,h,1,1,0,Ainv);
    vector<vector<double>> EtaDxDz = DCPSEcalcEta(x,y,z,neighbour,h,1,0,1,Ainv);
    vector<vector<double>> EtaDyDz = DCPSEcalcEta(x,y,z,neighbour,h,0,1,1,Ainv);

    // Calculating Lame Constants
    double lamda = v*E/((1+v)*(1-2*v));
    double Mu = E/(2*(1+v));

    // Equation Matrix Modelling
    VectorXd b;
    SparseMatrix<double> A;
    ArrangeEquationMatrix  (nx,ny,nz,nxi,nyi,nzi,dispX,dispY,dispZ,dispXi,dispYi,dispZi,forceX,forceY,
                            forceZ,forceXi,forceYi,forceZi,bx,by,bz,lamda,Mu,particlenumber,neighbour,
                            EtaDx,EtaDy,EtaDz,EtaDx2,EtaDy2,EtaDz2,EtaDxDy,EtaDxDz,EtaDyDz,h,b,A);
    
    // Solving Matrix Model
    vector<double> Ux(particlenumber);
    vector<double> Uy(particlenumber);
    vector<double> Uz(particlenumber);
    SolveMatrix(A,b,Ux,Uy,Uz,method_solver,preconditioner_iterative_solver,
                tolerance_iterative_solver,iteration_limit_iterative_solver,
                n_cpu_iterative_solver);

    cout <<"Max Positive X displacement: "<< findmax(Ux) << endl;
    cout <<"Max Negative X displacement: "<< findmin(Ux) << endl;
    cout <<"Max Positive Y displacement: "<< findmax(Uy) << endl;
    cout <<"Max Negative Y displacement: "<< findmin(Uy) << endl;
    cout <<"Max Positive Z displacement: "<< findmax(Uz) << endl;
    cout <<"Max Negative Z displacement: "<< findmin(Uz) << endl;
    
    // Solving for Displacement Gradient
    vector<double> dUxdX = DCPSEcalcGrad(Ux,neighbour,h,1,0,0,Ainv,EtaDx);
    vector<double> dUxdY = DCPSEcalcGrad(Ux,neighbour,h,0,1,0,Ainv,EtaDy);
    vector<double> dUxdZ = DCPSEcalcGrad(Ux,neighbour,h,0,0,1,Ainv,EtaDz);
    vector<double> dUydX = DCPSEcalcGrad(Uy,neighbour,h,1,0,0,Ainv,EtaDx);
    vector<double> dUydY = DCPSEcalcGrad(Uy,neighbour,h,0,1,0,Ainv,EtaDy);
    vector<double> dUydZ = DCPSEcalcGrad(Uy,neighbour,h,0,0,1,Ainv,EtaDz);
    vector<double> dUzdX = DCPSEcalcGrad(Uz,neighbour,h,1,0,0,Ainv,EtaDx);
    vector<double> dUzdY = DCPSEcalcGrad(Uz,neighbour,h,0,1,0,Ainv,EtaDy);
    vector<double> dUzdZ = DCPSEcalcGrad(Uz,neighbour,h,0,0,1,Ainv,EtaDz);

    // Solving for Stress and Strain
    vector<double> epsxx;
    vector<double> epsyy;
    vector<double> epszz;
    vector<double> epsxy;
    vector<double> epsxz;
    vector<double> epsyz;
    vector<double> sigmax;
    vector<double> sigmay;
    vector<double> sigmaz;
    vector<double> sigmaxy;
    vector<double> sigmaxz;
    vector<double> sigmayz;
    vector<double> vonmises;

    calcStressStrain(boundaryplotmethod,epsxx,epsyy,epszz,epsxy,epsxz,epsyz,sigmax,sigmay,sigmaz,sigmaxy,
                     sigmaxz,sigmayz,vonmises,dUxdX,dUxdY,dUxdZ,dUydX,dUydY,dUydZ,dUzdX,dUzdY,dUzdZ,lamda,
                     Mu,dispXi,dispYi,dispZi,forceXi,forceYi,forceZi);

    // Printing Result
    writeresult(boundaryplotmethod,x,y,z,vonmises,sigmax,sigmay,sigmaz,sigmaxy,sigmaxz,sigmayz,Ux,Uy,Uz,h,
                dispXi,dispYi,dispZi,forceXi,forceYi,forceZi,resultfilename);

    // Calculating Elapsed Time
    auto end = chrono::steady_clock::now();
    cout << "Total Elapsed Time: " <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms"<< endl;

}