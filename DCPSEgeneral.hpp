#ifndef DCPSEGENERAL_HPP
#define DCPSEGENERAL_HPP

#include <vector>
#include "eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

vector<MatrixXd> DCPSEcalcAinv (vector<double> x,
                                vector<double> y,
                                vector<double> z,
                                vector<vector<int>> neighbour,
                                double e);

vector<vector<double>> DCPSEcalcEta(vector<double> x,
                                    vector<double> y,
                                    vector<double> z,
                                    vector<vector<int>> neighbour,
                                    double e,
                                    int m,
                                    int n,
                                    int o,
                                    vector<MatrixXd> Ainv);

vector<double> DCPSEcalcGrad   (vector<double> w,
                                vector<vector<int>> neighbour,
                                double e,
                                int m,
                                int n,
                                int o,
                                vector<MatrixXd> Ainv,
                                vector<vector<double>> Eta);

#endif