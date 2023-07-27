#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>
#include "eigen-3.3.7/Eigen/Dense"
#include "eigen-3.3.7/Eigen/Sparse"
#include "eigen-3.3.7/Eigen/OrderingMethods"
#include "eigen-3.3.7/Eigen/IterativeLinearSolvers"
#include "MatrixOperation.hpp"

using namespace std;
using namespace Eigen;

void ArrangeEquationMatrix (vector<double> nx,
                            vector<double> ny,
                            vector<double> nz,
                            vector<int> nxi,
                            vector<int> nyi,
                            vector<int> nzi,
                            vector<double> dispX,
                            vector<double> dispY,
                            vector<double> dispZ,
                            vector<int> dispXi,
                            vector<int> dispYi,
                            vector<int> dispZi,
                            vector<double> forceX,
                            vector<double> forceY,
                            vector<double> forceZ,
                            vector<int> forceXi,
                            vector<int> forceYi,
                            vector<int> forceZi,
                            vector<double> bx,
                            vector<double> by,
                            vector<double> bz,
                            double lamda,
                            double Mu,
                            int particlenumber,
                            vector<vector<int>> neighbour,
                            vector<vector<double>> EtaDx,
                            vector<vector<double>> EtaDy,
                            vector<vector<double>> EtaDz,
                            vector<vector<double>> EtaDx2,
                            vector<vector<double>> EtaDy2,
                            vector<vector<double>> EtaDz2,
                            vector<vector<double>> EtaDxDy,
                            vector<vector<double>> EtaDxDz,
                            vector<vector<double>> EtaDyDz,
                            double e,
                            VectorXd& bresult,
                            SparseMatrix<double>& Aresult)
{
    cout << "Arranging Equation Matrix ... ";

    auto start =chrono::steady_clock::now();
    
    typedef Triplet<double> T;
    vector<T> tripletlist;

    VectorXd b = VectorXd::Zero(3*particlenumber);
    
    int row;
    int col;
    double data;

    double temp1;
    double temp2;
    double temp3;
    double temp4;
    double temp5;
    double temp6;
    double temp7;
    double temp8;
    double temp9;

    for (int i = 0; i < particlenumber; i++)
    {
        // Displacement Boundary X
        if (dispXi[i] == 1)
        {
            row = 3*i;
            col = 3*i;
            data = 1.0;

            tripletlist.push_back(T(row,col,data));

            b(3*i) = dispX[i];
        }

        // Displacement Boundary Y
        if (dispYi[i] == 1)
        {
            row = 3*i + 1;
            col = 3*i + 1;
            data = 1.0;

            tripletlist.push_back(T(row,col,data));

            b(3*i+1) = dispY[i];
        }

        // Displacement Boundary Z
        if (dispZi[i] == 1)
        {
            row = 3*i + 2;
            col = 3*i + 2;
            data = 1.0;

            tripletlist.push_back(T(row,col,data));

            b(3*i+2) = dispZ[i];
        }

        // Force Boundary X
        if (forceXi[i] == 1)
        {
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;
            
            for (int j = 0; j < neighbour[i].size(); j++)
            {
                temp1 = temp1 + (lamda + 2*Mu)*nx[i]*EtaDx[i][j]/e + Mu*ny[i]*EtaDy[i][j]/e + Mu*nz[i]*EtaDz[i][j]/e;
            
                temp2 = temp2 + lamda*nx[i]*EtaDy[i][j]/e + Mu*ny[i]*EtaDx[i][j]/e;
                
                temp3 = temp3 + lamda*nx[i]*EtaDz[i][j]/e + Mu*nz[i]*EtaDx[i][j]/e;
                
                row = 3*i;
                col = 3*neighbour[i][j];
                data = (lamda + 2*Mu)*nx[i]*EtaDx[i][j]/e + Mu*ny[i]*EtaDy[i][j]/e + Mu*nz[i]*EtaDz[i][j]/e;

                tripletlist.push_back(T(row,col,data));
            
                row = 3*i;
                col = 3*neighbour[i][j]+1;
                data = lamda*nx[i]*EtaDy[i][j]/e + Mu*ny[i]*EtaDx[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i;
                col = 3*neighbour[i][j]+2;
                data = lamda*nx[i]*EtaDz[i][j]/e + Mu*nz[i]*EtaDx[i][j]/e;

                tripletlist.push_back(T(row,col,data));
            
                if (j == neighbour[i].size()-1)
                    {
                        row = 3*i;
                        col = 3*i;
                        data = temp1;

                        tripletlist.push_back(T(row,col,data));
                        
                        row = 3*i;
                        col = 3*i+1;
                        data = temp2;

                        tripletlist.push_back(T(row,col,data));
                        
                        row = 3*i;
                        col = 3*i+2;
                        data = temp3;

                        tripletlist.push_back(T(row,col,data));
                    }
                    
            }
            b(3*i) = forceX[i];
        }

        // Force Boundary Y
        if (forceYi[i] == 1)
        {
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;
            
            for (int j = 0; j < neighbour[i].size(); j++)
            {               
                temp1 = temp1 + lamda*ny[i]*EtaDx[i][j]/e + Mu*nx[i]*EtaDy[i][j]/e;
            
                temp2 = temp2 + (lamda + 2*Mu)*ny[i]*EtaDy[i][j]/e + Mu*nx[i]*EtaDx[i][j]/e + Mu*nz[i]*EtaDz[i][j]/e;
                
                temp3 = temp3 + lamda*ny[i]*EtaDz[i][j]/e + Mu*nz[i]*EtaDy[i][j]/e;
                
                row = 3*i+1;
                col = 3*neighbour[i][j];
                data = lamda*ny[i]*EtaDx[i][j]/e + Mu*nx[i]*EtaDy[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+1;
                col = 3*neighbour[i][j]+1;
                data = (lamda + 2*Mu)*ny[i]*EtaDy[i][j]/e + Mu*nx[i]*EtaDx[i][j]/e + Mu*nz[i]*EtaDz[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+1;
                col = 3*neighbour[i][j]+2;
                data = lamda*ny[i]*EtaDz[i][j]/e + Mu*nz[i]*EtaDy[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                if (j == neighbour[i].size()-1)
                {
                    row = 3*i+1;
                    col = 3*i;
                    data = temp1;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+1;
                    col = 3*i+1;
                    data = temp2;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+1;
                    col = 3*i+2;
                    data = temp3;

                    tripletlist.push_back(T(row,col,data));
                }
            }
            b(3*i+1) = forceY[i];
        }

        // Force Boundary Z
        if (forceZi[i] == 1)
        {
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;
            
            for (int j = 0; j < neighbour[i].size(); j++)
            {               
                temp1 = temp1 + lamda*nz[i]*EtaDx[i][j]/e + Mu*nx[i]*EtaDz[i][j]/e;
            
                temp2 = temp2 + lamda*nz[i]*EtaDy[i][j]/e + Mu*ny[i]*EtaDz[i][j]/e;
                
                temp3 = temp3 + (lamda + 2*Mu)*nz[i]*EtaDz[i][j]/e + Mu*nx[i]*EtaDx[i][j]/e + Mu*ny[i]*EtaDy[i][j]/e;
                
                row = 3*i+2;
                col = 3*neighbour[i][j];
                data = lamda*nz[i]*EtaDx[i][j]/e + Mu*nx[i]*EtaDz[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+2;
                col = 3*neighbour[i][j]+1;
                data = lamda*nz[i]*EtaDy[i][j]/e + Mu*ny[i]*EtaDz[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+2;
                col = 3*neighbour[i][j]+2;
                data = (lamda + 2*Mu)*nz[i]*EtaDz[i][j]/e + Mu*nx[i]*EtaDx[i][j]/e + Mu*ny[i]*EtaDy[i][j]/e;

                tripletlist.push_back(T(row,col,data));
                
                if (j == neighbour[i].size()-1)
                {
                    row = 3*i+2;
                    col = 3*i;
                    data = temp1;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+2;
                    col = 3*i+1;
                    data = temp2;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+2;
                    col = 3*i+2;
                    data = temp3;

                    tripletlist.push_back(T(row,col,data));
                }
            }
            b(3*i+2) = forceZ[i];
        }
        
        // Inner Particles
        if ((forceXi[i] == 0) && ((forceYi[i] == 0) && ((forceZi[i] == 0) && ((dispXi[i] == 0) && ((dispYi[i] == 0) && (dispZi[i] == 0))))))
        {
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;
            temp4 = 0.0;
            temp5 = 0.0;
            temp6 = 0.0;
            temp7 = 0.0;
            temp8 = 0.0;
            temp9 = 0.0;
            
            for (int j = 0; j < neighbour[i].size(); j++)
            {
                // X direction
                temp1 = temp1 + (lamda + 2*Mu)*EtaDx2[i][j]/pow(e,2) + Mu*EtaDy2[i][j]/pow(e,2) + Mu*EtaDz2[i][j]/pow(e,2);
                
                temp2 = temp2 + (lamda + Mu)*EtaDxDy[i][j]/pow(e,2);
                
                temp3 = temp3 + (lamda + Mu)*EtaDxDz[i][j]/pow(e,2);

                row = 3*i;
                col = 3*neighbour[i][j];
                data = (lamda + 2*Mu)*EtaDx2[i][j]/pow(e,2) + Mu*EtaDy2[i][j]/pow(e,2) + Mu*EtaDz2[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i;
                col = 3*neighbour[i][j]+1;
                data = (lamda + Mu)*EtaDxDy[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i;
                col = 3*neighbour[i][j]+2;
                data = (lamda + Mu)*EtaDxDz[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                
                // Y direction
                temp4 = temp4 + (lamda + Mu)*EtaDxDy[i][j]/pow(e,2);
                
                temp5 = temp5 + (lamda + 2*Mu)*EtaDy2[i][j]/pow(e,2) + Mu*EtaDx2[i][j]/pow(e,2) + Mu*EtaDz2[i][j]/pow(e,2);
                
                temp6 = temp6 + (lamda + Mu)*EtaDyDz[i][j]/pow(e,2);
                
                row = 3*i+1;
                col = 3*neighbour[i][j];
                data = (lamda + Mu)*EtaDxDy[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+1;
                col = 3*neighbour[i][j]+1;
                data = (lamda + 2*Mu)*EtaDy2[i][j]/pow(e,2) + Mu*EtaDx2[i][j]/pow(e,2) + Mu*EtaDz2[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+1;
                col = 3*neighbour[i][j]+2;
                data = (lamda + Mu)*EtaDyDz[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                
                // Z direction
                temp7 = temp7 + (lamda + Mu)*EtaDxDz[i][j]/pow(e,2);
                
                temp8 = temp8 + (lamda + Mu)*EtaDyDz[i][j]/pow(e,2);
                
                temp9 = temp9 + (lamda + 2*Mu)*EtaDz2[i][j]/pow(e,2) + Mu*EtaDx2[i][j]/pow(e,2) + Mu*EtaDy2[i][j]/pow(e,2);
                
                row = 3*i+2;
                col = 3*neighbour[i][j];
                data = (lamda + Mu)*EtaDxDz[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+2;
                col = 3*neighbour[i][j]+1;
                data = (lamda + Mu)*EtaDyDz[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                row = 3*i+2;
                col = 3*neighbour[i][j]+2;
                data = (lamda + 2*Mu)*EtaDz2[i][j]/pow(e,2) + Mu*EtaDx2[i][j]/pow(e,2) + Mu*EtaDy2[i][j]/pow(e,2);

                tripletlist.push_back(T(row,col,data));
                
                
                if (j == neighbour[i].size()-1)
                {   
                    row = 3*i;
                    col = 3*i;
                    data = temp1;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i;
                    col = 3*i+1;
                    data = temp2;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i;
                    col = 3*i+2;
                    data = temp3;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+1;
                    col = 3*i;
                    data = temp4;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+1;
                    col = 3*i+1;
                    data = temp5;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+1;
                    col = 3*i+2;
                    data = temp6;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+2;
                    col = 3*i;
                    data = temp7;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+2;
                    col = 3*i+1;
                    data = temp8;

                    tripletlist.push_back(T(row,col,data));
                    
                    row = 3*i+2;
                    col = 3*i+2;
                    data = temp9;

                    tripletlist.push_back(T(row,col,data));

                }
            }
            b(3*i) = bx[i];
            b(3*i+1) = by[i];
            b(3*i+2) = bz[i];
        }
        
    }
    SparseMatrix<double> A(3*particlenumber,3*particlenumber);
    A.setFromTriplets(tripletlist.begin(), tripletlist.end());

    Aresult = A;
    bresult = b;

    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

void SolveMatrix   (SparseMatrix<double> A,
                    VectorXd b,
                    vector<double>& Ux,
                    vector<double>& Uy,
                    vector<double>& Uz,
                    int method,
                    int preconditioner_iterative_solver,
                    double tolerance_iterative_solver,
                    int iteration_limit_iterative_solver,
                    int n_cpu_iterative_solver)
{
    cout << "Solving Matrix Model ";

    auto start =chrono::steady_clock::now();

    A.makeCompressed();

    if (method == 1)
    {
        cout << "(Direct LU) ... " <<endl;
        
        SparseLU <SparseMatrix<double>, COLAMDOrdering<int>> solver;

        solver.compute(A);

        VectorXd U = solver.solve(b);

        int j = 0;
        int k = 0;
        int l = 0;

        for (int i = 0; i < b.rows(); i++)
        {
            if (i%3 == 0)
            {
                Ux[j] = U(i);
                j++;
            }
            if (i%3 == 1)
            {
                Uy[k] = U(i);
                k++;
            }
            if (i%3 == 2)
            {
                Uz[l] = U(i);
                l++;
            }
        }
    }
    else if (method == 2)
    {
        initParallel();
        omp_set_num_threads(n_cpu_iterative_solver);
        setNbThreads(n_cpu_iterative_solver);

        if (preconditioner_iterative_solver == 1)
        {
            cout << "(BiCGSTAB Jacobian) ... ";

            BiCGSTAB <SparseMatrix<double,RowMajor>> solver;

            solver.setTolerance(tolerance_iterative_solver);

            if (iteration_limit_iterative_solver == 0)
            {
                solver.setMaxIterations(b.rows()/5);
            }
            else
            {
                solver.setMaxIterations(iteration_limit_iterative_solver);
            }

            solver.compute(A);

            VectorXd U = solver.solve(b);
            cout <<endl<< "#iterations:     " << solver.iterations() << endl;
            cout << "estimated error: " << solver.error()      << endl;

            int j = 0;
            int k = 0;
            int l = 0;

            for (int i = 0; i < b.rows(); i++)
            {
                if (i%3 == 0)
                {
                    Ux[j] = U(i);
                    j++;
                }
                if (i%3 == 1)
                {
                    Uy[k] = U(i);
                    k++;
                }
                if (i%3 == 2)
                {
                    Uz[l] = U(i);
                    l++;
                }
            }
        }
        else if (preconditioner_iterative_solver == 2)
        {
            cout << "(BiCGSTAB IncompleteLUT) ... ";

            BiCGSTAB <SparseMatrix<double,RowMajor>,IncompleteLUT<double>> solver;
            solver.preconditioner().setDroptol(0.001);

            solver.setTolerance(tolerance_iterative_solver);

            if (iteration_limit_iterative_solver == 0)
            {
                solver.setMaxIterations(b.rows()/5);
            }
            else
            {
                solver.setMaxIterations(iteration_limit_iterative_solver);
            }

            solver.compute(A);

            VectorXd U = solver.solve(b);
            cout <<endl<< "#iterations:     " << solver.iterations() << endl;
            cout << "estimated error: " << solver.error()      << endl;

            int j = 0;
            int k = 0;
            int l = 0;

            for (int i = 0; i < b.rows(); i++)
            {
                if (i%3 == 0)
                {
                    Ux[j] = U(i);
                    j++;
                }
                if (i%3 == 1)
                {
                    Uy[k] = U(i);
                    k++;
                }
                if (i%3 == 2)
                {
                    Uz[l] = U(i);
                    l++;
                }
            }
        }
        else
        {
            cout << "Solving Error : Iterative Solver Preconditioner Undefined" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        cout << "Solving Error : Solver Method Undefined" << endl;
        exit(EXIT_FAILURE);
    }
    auto end = chrono::steady_clock::now();
    cout << "Matrix Solver Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}