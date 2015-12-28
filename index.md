---
layout: index
---
## 1. Overview
MkniX is a code designed for the simulation of nonlinear flexible multibody systems, including thermal effects. It is organised as a library but, as it already has all necesary I/O functions, building a binary is trivial. Several examples of applications using the MkniX library are provided in the `tests/` folder. The results can be viewed in the 3D interactive GUI, [MkniXPost](http://daniel-iglesias.github.io/mknixpost/).

## 2. Description
The library is divided in four modules: Systems, Core, Simulation, and Reader. Each of them has a specific functionality, as detailed in the following sections. 

## 2.1. Systems module
Hierarchical and recursive scheme of systems and subsystems with no limit in depth. Bodies are linked by joints and are independent of each other in terms of formulation and degrees of freedom, which allows mixing rigid and flexible solid bodies with different discretizations in the same model.

### 2.1.1. Rigid bodies
* 0D mass points.
* 1D bars and chains (in 2D or 3D spaces).
* 2D and 3D generic solid bodies:
	* CoG and main in\ertia vectors.
	* Meshed 3D body.
* CAD (STEP, STL) body [TBD].

### 2.1.2. Continuum bodies
* 2D and 3D Thermal rigid bodies.
* 2D and 3D Flexible nonlinear bodies using ANCF with FEM or Meshfree methods.
* Weakly coupled thermo-mechanical flexible bodies.

### 2.1.3. Joints and contacts
* Constraint methods:
	+ Penalty.
	+ Augmented Lagrangian.
* Types of joints (can be defined at any point in space):
	+ Spherical.
	+ Spherical with clearances.
	+ Revolute (pin joint) [in-work].
	+ Cylindrical slider [in-work].
	+ Prismatic slider [in-work].
	+ Planar [in-work].
* Contact detection:
	+ Auto generated contact elements by α-shapes.
	+ Generalized normal contact force (penalty potential) w/o friction.
	+ 1D wire boundary interfaces.
	+ 2D faceted boundary interfaces [TBD].
	
### 2.1.4 Loads
* Mechanical:
	+ Generalized acceleration (gravity).
	+ Point loads (non-following).
	+ Surface loads [in-work].
	+ Volumetric (body) loads.
* Thermal:
	+ Heat flux (surface) [in-work]:
		- Explicitly defined in any boundary.
		- Convection effects (HTC and Tbulk as inputs).
	+ Heat flow (volume):
		- Explicitly defined in local or global frames.
		- Calculated from radiation maps (see below).
	+ Beam particle models:
		- Bi-gaussian distributions defined by a vector in space.
* Radiation [in-work]
	+ 2D or 3D gamma maps.
	+ Fixed in space or associated to a moving frame of any body.
	+ Accumulated damage monitoring (non-shadowed body).
	
## 2.2. Core module
Implements the continuum mechanics residual contributions for solid bodies. It also provides the generic point and node classes used by rigid bodies and joints.

### 2.2.1. Formulations
* Rigid body dynamics.
* Linear elasticity (small displacements).
* Nonlinear elasticity in a total Lagrangian approach (large strains and large displacements).
* Thermal formulation for rigid bodies (small /large displacements).
* Linear thermoelasticity.
* Nonlinear thermoelasticity in a total Lagrangian form (static/dynamics weakly coupled).

### 2.2.2. Discretization methods for continuum bodies
* Finite element:
	+ 2D triangle, multiple order.
	+ 3D tetrahedron, multiple order.
* Meshfree Galerkin methods:
	+ Radial Basis interpolated functions.
	+ Moving least squares approximated shape functions.

### 2.2.3. Material models for all elastic/thermoelastic formulations
* Isotropic elastic (linear) Hooke model
* Isotropic hyperelastic (nonlinear) Saint-Venant Kirchhoff
* All materials available in 3D and 2D (plain stress [TBD] or strain) spaces
* Nonlinear properties (i.e. temperature dependent)

## 2.3. Simulation module
It is the main interface of the Mknix code. The module controls the time stepping, integration method, and the convergence criteria for the nonlinear solver. Most procedures are provided by the numerical library LMX. It also owns and controls the global matrices and the top level systems.

### 2.3.1 Analysis types and integrators
* Static and quasi-static.
* Implicit dynamics:
	+ Central differences.
	+ Newmark.
	+ HTT-α [in work].
	+ Linear multistep integrators:
		- BDF (Order 1-5).
		- Adams Moulton (O 1-5).
* Explicit dynamics:
	+ Adams Bashforth (O 1-5).
	+ Dual order integration for thermo-mechanical dynamic analysis.

### 2.3.2. Linear solvers
* Direct:
	- Gauss elimination for dense matrices
	- LU decomposition methods (SuperLU and Gmm::lu_solve)
* Iterative Krylov methods:
	- Conjugate gradient for symmetric matrices
	- Preconditioned conjugate gradient
	- Generalized minimal residual method (GMRES)

### 2.3.3. Nonlinear solvers
* Newton-Raphson
* Tangent method (quasi-Newton, BFGS) [TBD]

## 2.4. Reader module
Controls the input, output by bridging the Analysis’ solvers required functions with the Systems’ formulation. 

### 2.4.1 Input readers
One input file referring to other geometry and mesh files is usually used. The file is structured in sections, each one read by independent object classes and member functions:
* General class reader for Systems´ structure, materials, loads and Analysis.
* Rigid body reader, including mass properties and geometry.
* Flexbody reader, with formulation type, geometry (with or without mesh file), assigned material number, and discretization parameters.
* Constraints reader, for joints definition between bodies.

Accepted geometry files:
* Gmsh format.
* Salome mesh format.
* STL Cad models [in work].

### 2.4.2 Output writers
* Results summary including displacements and temperatures, nodal stress, energies, etc.
* Simulation times (real stopwatch output).
* Reader output, used mainly for checking and debugging.
* Shapefunctions data for checking consistency of the approximation.

## 3. Building and testing
See README.md for a description of requirements and procedure for building the MkniX library and using it for creating binary executables.

## 4. License
MkniX library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

Tests are distributed under GPL version 3 of the License, or (at your option) any later version.
