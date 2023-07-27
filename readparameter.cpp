#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "readparameter.hpp"

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
                                int& n_cpu_iterative_solver)
{
    ifstream inFile;
    inFile.open("solverparameter.txt");

    if (!inFile) 
    {
    cerr << "Unable to open file solverparameter.txt";
    exit(1);   // call system to stop
    }

    string line;
    vector<string> data;

    for (int i = 0; i < 12; i++)
    {
        getline(inFile,line);
        string s;
        vector<string> v;
        stringstream ss(line);
        while (getline(ss, s, ':'))
            v.push_back(s); // store token string in the vector
        data.push_back(v[1].substr(1, v[1].size()-2)); // don't use the first char (reserved for space)
    }
    inFile.close();

    geomfilename = data[0];
    geomparamfilename = data[1];
    resultfilename = data[2];
    neighbourmethod = stoi(data[3]);
    rcScale = stod(data[4]);
    neighbournumber = stoi(data[5]);
    boundaryplotmethod = stoi(data[6]);
    method = stoi(data[7]);
    preconditioner_iterative_solver = stoi(data[8]);
    tolerance_iterative_solver = stod(data[9]);
    iteration_limit_iterative_solver = stoi(data[10]);
    n_cpu_iterative_solver = stoi(data[11]);

    cout << "Geometry Data File: " <<geomfilename<<endl;
    cout << "Geometry Parameter File: " <<geomparamfilename<<endl;
    cout << "Result File: " <<resultfilename<<endl;
    cout << "Neighbourfind Method: " <<neighbourmethod<<endl;
    cout << "Cutoff Radius Scale: " <<rcScale<<endl;
    cout << "Desired Neighbour Number: " <<neighbournumber<<endl;
    cout << "Boundary Plot Method: " <<boundaryplotmethod<<endl;
    cout << "Linear System Solver Method: "<<method<<endl;
    cout << "Iterative Solver Preconditioner: "<<preconditioner_iterative_solver<<endl;
    cout << "Iterative Solver Tolerance: "<<tolerance_iterative_solver<<endl;
    cout << "Iterative Solver Iteration Limit: "<<iteration_limit_iterative_solver<<endl;
    cout << "Iterative Solver Multi-Threading CPU: "<<n_cpu_iterative_solver<<endl<<endl;

}