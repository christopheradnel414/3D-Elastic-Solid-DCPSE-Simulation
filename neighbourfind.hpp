#ifndef NEIGHBOURFIND_HPP
#define NEIGHBOURFIND_HPP

#include <vector>

using namespace std;

void zip   (vector<int> sorted,
            vector<double> sorter,
            vector<pair<int,double>> &zipped);

void unzip (vector<int> &sorted,
            vector<double> &sorter,
            vector<pair<int,double>> zipped);


void neighbourfind (vector<vector<int>>& neighbour,
                    vector<double> x,
                    vector<double> y,
                    vector<double> z,
                    double rc);

void neighbourfindlimited  (vector<vector<int>>& neighbour,
                            vector<double> x,
                            vector<double> y,
                            vector<double> z,
                            double rc,
                            int targetnumber);

#endif