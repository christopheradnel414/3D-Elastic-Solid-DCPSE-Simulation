#ifndef READPARAMETER_HPP
#define READPARAMETER_HPP

#include <string>
#include <vector>

using namespace std;

void readsolverparameterfile   (int& boundaryplotmethod,
                                string& geomfilename,
                                string& geomparamfilename,
                                string& resultfilename,
                                int& neighbourmethod,
                                int& neighbournumber,
                                double& rcScale,
                                int& method,
                                int& preconditioner_iterative_solver,
                                double& tolerance_iterative_solver,
                                int& iteration_limit_iterative_solver,
                                int& n_cpu_iterative_solver);
#endif