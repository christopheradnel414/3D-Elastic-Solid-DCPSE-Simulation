# 3D-Elastic-Solid-DCPSE-Simulation-
This repository is an extension of my 2D Linear Elasticity Solver project (https://github.com/christopheradnel414/2D-Elastic-Solid-DCPSE-Simulation) using DCPSE (Discretization Corrected Particle Strength Exchange) to the 3D domain. The structure of this repository is largely similar to the 2D variant with the exception of the result visualization code, which is using Matlab here.

<img height="200" alt="Screenshot 2023-07-27 at 01 41 52" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/5c4fefec-4dc4-4e98-9f0a-f1adc8308847">
<img height="200" alt="Screenshot 2023-07-27 at 01 42 37" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/d1102c33-31ef-4102-a29d-c8b783a08384">
<img height="200" alt="Screenshot 2023-07-27 at 01 43 04" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/73cdbbb2-c255-4954-9f84-d14165c83a41">
<img height="200" alt="Screenshot 2023-07-27 at 01 43 33" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/9c1c06c7-b3b7-47ab-8d09-ed260dc2e9fc">

This repository mainly consists of 3 parts:
1. 3D Linear Elastostatics DCPSE Source Code (main directory)
2. Example Geometry Generation Code (geometries folder)
3. Result Plotter Code (result folder)

# Main Solver
The main Solver is written in C++ and utilizes the Eigen (https://eigen.tuxfamily.org) library as the sparse linear algebra engine. The main solver takes into input a .txt file containing point clouds, each anotated with the constraint and condition of the nodes (inner node, dirichlet boundary node, neumann boundary node, etc) and outputs a .txt file of the simulation results containing the node displacements, stresses, and strains. This solver works by solving the Static Navier-Cauchy equations (https://en.wikipedia.org/wiki/Linear_elasticity) using DCPSE discretization method. The 3D Static Navier-Cauchy equations can be written as follows:

1. Static Navier-Cauchy for Inner Nodes:
<img height="150" alt="Screenshot 2023-07-27 at 00 42 36" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/29a2d53f-3812-4522-a876-ef5df9c4b80e">

2. Static Navier-Cauchy for Neumann Nodes (Traction/Force Boundary):
<img height="150" alt="Screenshot 2023-07-27 at 00 42 46" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/81041bd6-2f87-450e-94a0-3a000ed93ba0">

3. Static Navier-Cauchy for Dirichlet Nodes (Displacement Boundary):
<img height="100" alt="Screenshot 2023-07-27 at 00 42 58" src="https://github.com/christopheradnel414/3D-Elastic-Solid-DCPSE-Simulation/assets/41734037/20cfaae0-bac4-4e53-b32c-ccf551d9a302">

There are many methods that can be used to solve the Static Navier-Cauchy Equations (e.g, finite difference, finite element, SPH, etc) which I discussed in more detail on my journal (https://journals.itb.ac.id/index.php/jets/article/view/17606). Here, we are only interested in using DCPSE as our discretization method to solve the Navier-Cauchy equation. DCPSE operator as introduced by Schrader,2010 (https://publications.mpi-cbg.de/Schrader_2010_4838.pdf) approximates the value and spatial gradients of a field function as a linear combination of the neighbouring nodes function values. This is written as follows:

<img height="200" alt="Screenshot 2023-07-27 at 01 03 56" src="https://github.com/christopheradnel414/2D-Elastic-Solid-DCPSE-Simulation/assets/41734037/d503a3fb-bb3c-4904-a4cd-be936ec17188">

The detail on the computation of the DCPSE operator is given on the journal. However, the main idea is that to approximate the spatial gradient at any particular point, it is simply the weighted sum of the neighbouring nodes' function value. The weights $\eta_{ij}$ can be easily computed as detailed on the journal and this leaves us only to solve the huge system of linear equations after discretizing the Navier-Cauchy equations using DCPSE. To solve this huge system system of equations, we are using Eigen library's sparse matrix solver capability.

To compile the solver, you can simply run the Makefile in the main directory using the "make" command:
```
make
```

Afterwards, you can run the compiled executable as follows:
```
./DCPSE_SOLID_3D_Solver
```

To modify the solver parameters (e.g, size of cutoff radius, input location, output location, linear solver method, etc), you can modify the solverparameter.txt file.

# Geometry Generation
Here, geometry is represented as point clouds with annotations to denote whether it is an inner node, dirichlet node, or neumann node. Examples of the geometry generation scripts (Python) are given in the geometries folder. The main output of these geometry generation scripts are "geom_data.txt" and "geom_param.txt". "geom_data.txt" stores the point cloud data representing the geometry along with the node annotations while "geom_param.txt" is a very short .txt file containing the material properties of the structure.

# Result Plotter
To visualize the results of the simulation, we are using Matlab to make use of the 3D scatter point viewer capabilities (especially the drag-and-rotate capabilities). This viewer script is written as Viewer_Result.m.
