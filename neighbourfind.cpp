#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "neighbourfind.hpp"

using namespace std;

void zip   (vector<int> sorted,
            vector<double> sorter,
            vector<pair<int,double>> &zipped)
{
    for (int i = 0; i < sorted.size(); i++)
    {
        zipped.push_back(make_pair(sorted[i],sorter[i]));
    }
}

void unzip (vector<int> &sorted,
            vector<double> &sorter,
            vector<pair<int,double>> zipped)
{
    for (int i = 0; i < sorted.size(); i++)
        {
            sorted[i] = zipped[i].first;
            sorter[i] = zipped[i].second;
        }
}


void neighbourfind (vector<vector<int>>& neighbour,
                    vector<double> x,
                    vector<double> y,
                    vector<double> z,
                    double rc)
{
    cout <<"Finding Particle's Neighbour (unlimited) ...";

    auto start = chrono::steady_clock::now();

    double dX;
    double dY;
    double dZ;
    double dR;

    int particlenumber = x.size();
    
    for (int i = 0; i < particlenumber-1; i++)
    {
        for (int j = i+1; j < particlenumber; j++)
        {
            dX = x[i] - x[j];
            dY = y[i] - y[j];
            dZ = z[i] - z[j];
            dR = sqrt(pow(dX,2) + pow(dY,2) + pow(dZ,2));
            if (dR < rc)
            {
                neighbour[i].push_back(j);
                neighbour[j].push_back(i);
            }

        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

void neighbourfindlimited  (vector<vector<int>>& neighbour,
                            vector<double> x,
                            vector<double> y,
                            vector<double> z,
                            double rc,
                            int targetnumber)
{
    cout <<"Finding Particle's Neighbour (limited to "<<targetnumber<<" neighbours) ...";

    auto start = chrono::steady_clock::now();

    double dX;
    double dY;
    double dZ;
    double dR;

    int particlenumber = x.size();

    vector<vector<double>> dRdata(particlenumber);
    vector<pair<int,double>> zipped;

    vector<vector<int>> neighbourtemp(particlenumber);
    
    for (int i = 0; i < particlenumber-1; i++)
    {
        for (int j = i+1; j < particlenumber; j++)
        {
            dX = x[i] - x[j];
            dY = y[i] - y[j];
            dZ = z[i] - z[j];
            dR = sqrt(pow(dX,2) + pow(dY,2) + pow(dZ,2));
            if (dR < rc)
            {
                neighbourtemp[i].push_back(j);
                neighbourtemp[j].push_back(i);
                dRdata[i].push_back(dR);
                dRdata[j].push_back(dR);
            }
        }
    }

    for (int i = 0; i < particlenumber; i++)
    {
        if (neighbourtemp[i].size() < targetnumber)
        {
            cout<<"Neighbour Search Error: Not enough neighbour"<<endl;
            exit(EXIT_FAILURE);
        }

        zip(neighbourtemp[i],dRdata[i],zipped);
        
        sort(begin(zipped), end(zipped), 
        [&](const auto& a, const auto& b)
        {
            return a.second < b.second;
        });

        unzip(neighbourtemp[i],dRdata[i],zipped);
        zipped.clear();

        for (int j = 0; j < targetnumber; j++)
        {
            neighbour[i].push_back(neighbourtemp[i][j]);
        }
    }
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}