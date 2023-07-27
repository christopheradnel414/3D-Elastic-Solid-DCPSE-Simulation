#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include "readgeom.hpp"
#include "functions.hpp"

using namespace std;

void readgeom  (double& E,
                double& v,
                double& h,
                vector<double>& x,
                vector<double>& y,
                vector<double>& z,
                vector<double>& nx,
                vector<double>& ny,
                vector<double>& nz,
                vector<int>& nxi,
                vector<int>& nyi,
                vector<int>& nzi,
                vector<double>& dispX,
                vector<double>& dispY,
                vector<double>& dispZ,
                vector<int>& dispXi,
                vector<int>& dispYi,
                vector<int>& dispZi,
                vector<double>& forceX,
                vector<double>& forceY,
                vector<double>& forceZ,
                vector<int>& forceXi,
                vector<int>& forceYi,
                vector<int>& forceZi,
                vector<double>& bx,
                vector<double>& by,
                vector<double>& bz,
                string geomfilename,
                string geomparamfilename)
{

    // Reading the geometry data file
    cout <<"Reading geometry file ... ";

    auto start = chrono::steady_clock::now();

    ifstream file;

    vector<double> param_temp;
    double reader;
    file.open(geomparamfilename);

    while (file >> reader)
    {
        param_temp.push_back(reader);
    }
    file.close();

    E = param_temp[0];
    v = param_temp[1];
    h = param_temp[2];
    
    string xs,ys,zs,nxs,nys,nzs,dispXs,dispYs,dispZs,forceXs,forceYs,forceZs,bxs,bys,bzs;
    file.open(geomfilename);

    while (file>>xs>>ys>>zs>>nxs>>nys>>nzs>>dispXs>>dispYs>>dispZs>>forceXs>>forceYs>>forceZs>>bxs>>bys>>bzs)
    {
        // reading x and y and z
        x.push_back(stod(xs));
        y.push_back(stod(ys));
        z.push_back(stod(zs));

        // reading nx
        if (is_number(nxs))
        {
            nx.push_back(stod(nxs));
            nxi.push_back(1);
        }
        else
        {
            nx.push_back(0.0);
            nxi.push_back(0);
        }
        
        // reading ny
        if (is_number(nys))
        {
            ny.push_back(stod(nys));
            nyi.push_back(1);
        }
        else
        {
            ny.push_back(0.0);
            nyi.push_back(0);
        }

        // reading nz
        if (is_number(nzs))
        {
            nz.push_back(stod(nzs));
            nzi.push_back(1);
        }
        else
        {
            nz.push_back(0.0);
            nzi.push_back(0);
        }

        // reading dispX
        if (is_number(dispXs))
        {
            dispX.push_back(stod(dispXs));
            dispXi.push_back(1);
        }
        else
        {
            dispX.push_back(0.0);
            dispXi.push_back(0);
        }

        // reading dispY
        if (is_number(dispYs))
        {
            dispY.push_back(stod(dispYs));
            dispYi.push_back(1);
        }
        else
        {
            dispY.push_back(0.0);
            dispYi.push_back(0);
        }

        // reading dispZ
        if (is_number(dispZs))
        {
            dispZ.push_back(stod(dispZs));
            dispZi.push_back(1);
        }
        else
        {
            dispZ.push_back(0.0);
            dispZi.push_back(0);
        }

        // reading forceX
        if (is_number(forceXs))
        {
            forceX.push_back(stod(forceXs));
            forceXi.push_back(1);
        }
        else
        {
            forceX.push_back(0.0);
            forceXi.push_back(0);
        }

        // reading forceY
        if (is_number(forceYs))
        {
            forceY.push_back(stod(forceYs));
            forceYi.push_back(1);
        }
        else
        {
            forceY.push_back(0.0);
            forceYi.push_back(0);
        }

        // reading forceZ
        if (is_number(forceZs))
        {
            forceZ.push_back(stod(forceZs));
            forceZi.push_back(1);
        }
        else
        {
            forceZ.push_back(0.0);
            forceZi.push_back(0);
        }

        // reading bx and by and bz
        bx.push_back(stod(bxs));
        by.push_back(stod(bys));
        bz.push_back(stod(bzs));
        
    }
    file.close();
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
    cout << "Number of Particles: " <<x.size()<<endl;
}