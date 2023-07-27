#ifndef MATRIXOPERATION_HPP
#define MATRIXOPERATION_HPP
#define NDEBUG
#define EIGEN_NO_STATIC_ASSERT

#include <vector>
#include "eigen-3.3.7/Eigen/Dense"
#include "eigen-3.3.7/Eigen/Sparse"
#include "eigen-3.3.7/Eigen/OrderingMethods"
#include "eigen-3.3.7/Eigen/IterativeLinearSolvers"

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
                            SparseMatrix<double>& Aresult);

void SolveMatrix   (SparseMatrix<double> A,
                    VectorXd b,
                    vector<double>& Ux,
                    vector<double>& Uy,
                    vector<double>& Uz,
                    int method,
                    int preconditioner_iterative_solver,
                    double tolerance_iterative_solver,
                    int iteration_limit_iterative_solver,
                    int n_cpu_iterative_solver);

#endif