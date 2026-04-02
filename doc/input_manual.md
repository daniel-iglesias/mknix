# mknix Input File Manual

## Overview

An mknix input file is a plain-text file read sequentially by keyword. Keywords are case-sensitive. Comments begin with `//` and extend to the end of the line. The general structure is:

```
TITLE <name>
WORKINGDIR <path>
DIMENSION <dim>
MATERIALS
  ...
ENDMATERIALS
SYSTEM <name>
  ...
ENDSYSTEM
ANALYSIS
  ...
ENDANALYSIS
```

---

## Top-Level Commands

### `TITLE`
```
TITLE <name>
```
Sets the simulation title (single token, no spaces).

---

### `WORKINGDIR`
```
WORKINGDIR <path>
```
Changes the working directory. All subsequent file paths are relative to this directory.

---

### `DIMENSION`
```
DIMENSION <2|3>
```
Sets the spatial dimension of the problem.

---

### `GRAVITY`
```
GRAVITY <gx> <gy> <gz>
```
Sets the gravity vector components.

---

### `CONTACT`
```
CONTACT <GLOBAL|NONE>
```
Sets the contact detection strategy.

---

### `VISUALIZATION`
```
VISUALIZATION <ON|OFF>
```
Enables or disables visualization output.

---

### `OUTPUT`
```
OUTPUT MATRICES
```
Enables output of system matrices to file.

---

### `SMOOTHING`
```
SMOOTHING <OFF|LOCAL|CONSTANT|GLOBAL>
```
Sets the stress-smoothing strategy.

---

### `INITIALTEMPERATURE`
```
INITIALTEMPERATURE <value>
```
Sets the global initial temperature for all nodes.

---

## `MATERIALS` Section

```
MATERIALS
  // one or more material definitions
ENDMATERIALS
```

### `PLSTRAIN` (Plane-strain mechanical material)
```
PLSTRAIN <mat_id> <E> <nu> <density>
```
- `mat_id` – integer material identifier
- `E` – Young's modulus
- `nu` – Poisson's ratio
- `density` – mass density

### `THERMAL` (Thermal material, constant properties)
```
THERMAL <mat_id> <Cp> <kappa> <beta> <density>
```
- `Cp` – heat capacity
- `kappa` – thermal conductivity
- `beta` – thermal expansion coefficient
- `density` – mass density

### `POROSITY` (Porous medium extension; requires a `THERMAL` entry for the same `mat_id`)
```
POROSITY <mat_id> <Rsl> <T_bulk>
```
- `Rsl` – fluid-solid thermal resistance
- `T_bulk` – bulk fluid temperature

### `FILES` (Temperature-dependent thermal properties from files)

```
FILES
  CAPACITY     <mat_id> <filename>
  CONDUCTIVITY <mat_id> <filename>
  RESISTANCE   <mat_id> <filename>
ENDFILES
```

Each file is a two-column, whitespace-separated text file with no header:

| File type | Column 1 | Column 2 |
|---|---|---|
| `CAPACITY` | temperature | heat capacity |
| `CONDUCTIVITY` | temperature | thermal conductivity |
| `RESISTANCE` | distance | interface thermal resistance |

---

## `SYSTEM` Section

```
SYSTEM <name>
  RIGIDBODIES ... ENDRIGIDBODIES
  FLEXBODIES  ... ENDFLEXBODIES
  BODYPOINTS  ... ENDBODYPOINTS
  JOINTS      ... ENDJOINTS
  LOADS       ... ENDLOADS
  ENVIRONMENT ... ENDENVIRONMENT
  MOTION      <node_id> ... ENDMOTION
  SCALE       <sx> <sy> <sz>
  MIRROR      <x|y|z>
  SHIFT       <dx> <dy> <dz>
  SIGNALS     ... ENDSIGNALS
ENDSYSTEM
```

`SCALE`, `MIRROR`, and `SHIFT` apply immediately to all nodes defined so far.

---

## `RIGIDBODIES` Section

```
RIGIDBODIES
  PENALTY | AUGMENTED
  ALPHA <value>
  MASSPOINT <name> ... ENDMASSPOINT
  BAR       <name> ... ENDBAR
  CHAIN     <name> ... ENDCHAIN
  GENERIC2D <name> ... ENDGENERIC2D
  GENERIC3D <name> ... ENDGENERIC3D
ENDRIGIDBODIES
```

### `PENALTY` / `AUGMENTED`
Sets the constraint enforcement method (penalty or augmented Lagrange). Default is penalty.

### `ALPHA`
```
ALPHA <value>
```
Sets the penalty / augmented-Lagrange stiffness factor.

---

### `MASSPOINT`
```
MASSPOINT <name>
  NODEA <x> <y> <z>
  MASS  <value>
ENDMASSPOINT
```

---

### `BAR`
```
BAR <name>
  NODEA  <x> <y> <z>
  NODEB  <x> <y> <z>
  MASS   <value>
  OUTPUT ENERGY
ENDBAR
```

---

### `CHAIN`
```
CHAIN <name>
  NODEA      <x> <y> <z>
  NODEB      <x> <y> <z>
  MASS       <value>
  SEGMENTS   <n>
  LENGTH     <value>
  TIMELENGTH <time> <length>  // repeatable
  OUTPUT     ENERGY
ENDCHAIN
```
`TIMELENGTH` may appear multiple times to define a time-varying chain length.

---

### `GENERIC2D`
```
GENERIC2D <name>
  MASS     <value>
  IXX      <value>
  IYY      <value>
  IXY      <value>
  DENSITY  <factor>
  POSITION <xCoG> <yCoG> <theta>
  TRIANGLES FILE <meshfile>
  OUTPUT   ENERGY
ENDGENERIC2D
```
`meshfile` is a mesh file with format: first line `<nNodes> <nCells>`, then `<id> <x> <y> <z>` per node, then connectivity.

---

### `GENERIC3D`
```
GENERIC3D <name>
  MASS        <value>
  IXX <value>  IYY <value>  IZZ <value>
  IXY <value>  IYZ <value>  IXZ <value>
  DENSITYFACTOR <factor>
  POSITION    <xCoG> <yCoG> <zCoG>
  OUTPUT      ENERGY
ENDGENERIC3D
```

---

## `FLEXBODIES` Section

```
FLEXBODIES
  MESHFREE <name> ... ENDMESHFREE
  FEMESH   <name> ... ENDFEMESH
  SHARENODES <nameA> <nameB>
ENDFLEXBODIES
```

`SHARENODES` makes two bodies share their node lists (for coupled domains).

### `MESHFREE` / `FEMESH`

Both use the same sub-keywords:

```
MESHFREE <name>
  FORMULATION    <THERMAL|MECHANICAL|THERMOMECHANICAL>
  METHOD         <EFG|RPIM>
  OUTPUT         <STRESS|ENERGY>
  INITIALTEMPERATURE <value>
  BOUNDARYGROUP  <bgName> ... ENDBOUNDARYGROUP
  NODES          ...
  CELLS          ...
  MESH           ...
ENDMESHFREE
```

#### `BOUNDARYGROUP`
```
BOUNDARYGROUP <name>
  METHOD <formulation> <nGPs> <alpha>
  FILE   <meshfile>
ENDBOUNDARYGROUP
```

#### `NODES` – Inline definition

**Rectangular patch:**
```
NODES
RECTANGULAR <nx> <ny> <x1> <y1> <x2> <y2>
```
Creates `(nx+1)*(ny+1)` nodes on a regular grid from `(x1,y1)` to `(x2,y2)`.

**Grid from file (triangles or quads):**
```
NODES
GRID TRIANGLES <meshfile>
NODES
GRID QUADS     <meshfile>
```

#### `CELLS` – Integration cell definition

```
CELLS
<mat_id> <nGPs> <alpha> <dc>
RECTANGULAR <mat_id> <dcx> <dcy> <nx> <ny> <x1> <y1> <x2> <y2>
```

#### `MESH` – Combined node + cell from file

```
MESH
<mat_id> <nGPs> <alpha>
TRIANGLES <meshfile>
QUADS     <meshfile>
```

---

## `BODYPOINTS` Section

```
BODYPOINTS
  <bodyName> <x> <y> <z>                         // rigid body point
  <bodyName> <x> <y> <z> <alpha> <dc>             // flex body point
ENDBODYPOINTS
```

Adds extra nodes to already-defined bodies.

---

## `JOINTS` Section

```
JOINTS
  PENALTY | AUGMENTED
  ALPHA <value>
  SPHERICAL       <name> ... ENDSPHERICAL
  DISTANCE        <name> ... ENDDISTANCE
  AXIS            <name> ... ENDAXIS
  CLEARANCE       <name> ... ENDCLEARANCE
  THERMALSPHERICAL <name> ... ENDTHERMALSPHERICAL
ENDJOINTS
```

All joint types share `NODEA` and `NODEB` sub-commands:
```
NODEA <bodyName>.<nodeId>
NODEB <bodyName>.<nodeId>
```
Use `GROUND` as the body name to fix a node to the ground frame.

### `AXIS` – additional sub-command
```
DIRECTION <x|y|z>
```

### `CLEARANCE` – additional sub-command
```
TOLERANCE <value>
```

---

## `LOADS` Section

```
LOADS
  FORCE           <body>.<node> <fx> <fy> <fz>
  THERMALFLUENCE  <body>.<node> <value>
  THERMALOUTPUT   <body>.<node>
  THERMALOUTPUT   MAX_INTERFACE_TEMP
  THERMALBODY     <bodyName> VALUE <value>
  THERMALBODY     <bodyName> FILE  <filename>
  THERMALFLUX1D   <body>.<boundaryGroup> ... ENDTHERMALFLUX1D
  RADIATION       ... ENDRADIATION
ENDLOADS
```

---

### `THERMALBODY`

Applies a volumetric heat source to a thermal body.

**Constant value:**
```
THERMALBODY <bodyName> VALUE <value>
```

**1-D spatial distribution from file:**
```
THERMALBODY <bodyName> FILE <filename>
```
The file contains two whitespace-separated columns (X coordinate, load value), one row per point:
```
0.0   1000.0
0.05  1500.0
0.10  1200.0
```

**2-D spatial distribution from file (regular grid):**
```
THERMALBODY <bodyName> FILE <filename>
```
The file uses a matrix layout. The first row lists the `key2` coordinate values (first cell is ignored). Each subsequent row starts with the `key1` coordinate, followed by the load values for each `key2`:
```
0.0   0.0    0.05   0.10
0.0   1000   1100   1200
0.05  1500   1600   1700
0.10  1200   1300   1400
```
The parser automatically distinguishes 1-D from 2-D based on the number of columns.

---

### `THERMALFLUX1D`

Applies a 1-D heat flux to a boundary group.

```
THERMALFLUX1D <bodyName>.<boundaryGroupName>
  FILE     <fluxFilename>
  TIMEFILE <timeFilename>
  SCALE    <factor>
ENDTHERMALFLUX1D
```

- `FILE` – two-column file: coordinate vs flux value
- `TIMEFILE` – two-column file: time vs scale factor
- `SCALE` – multiplies all flux values by a constant factor

---

### `RADIATION`

```
RADIATION
  STATIC3D | STATIC2D
  SKIPLINES   <n>
  LENGTHFACTOR <factor>
  SCALEAXIS   <sx> <sy> <sz>
  DOSEFACTOR  <factor>
  MAPFILE     <filename>
ENDRADIATION
```

`MAPFILE` contains rows of `<x> <y> [<z>] <dose>`. `STATIC3D` expects 4 columns; `STATIC2D` expects 3 columns (z is assumed 0). `SKIPLINES` discards header lines.

---

## `MOTION` Section

```
MOTION <nodeId>
  TIMECONF <time> <ux> <uy> <uz>   // repeatable
ENDMOTION
```

Prescribes time-varying displacement on a ground node. Multiple `TIMECONF` entries form a piecewise-linear trajectory.

---

## `SIGNALS` Section

```
SIGNALS
  MECHANICALOUTPUT <signalName> <bodyName>.<nodeId>
  MECHANICALINPUT  <signalName> <constraintName1> [<constraintName2> ...]
ENDSIGNALS
```

---

## `ANALYSIS` Section

```
ANALYSIS
  STATIC                  ... ENDSTATIC
  DYNAMIC                 ... ENDDYNAMIC
  THERMALSTATIC           ... ENDTHERMALSTATIC
  THERMALDYNAMIC          ... ENDTHERMALDYNAMIC
  THERMOMECHANICALDYNAMIC ... ENDTHERMOMECHANICALDYNAMIC
ENDANALYSIS
```

All analysis types share `EPSILON` and `TIME`. Transient types also require `INTEGRATOR`.

### `STATIC`
```
STATIC
  EPSILON <tol>
  TIME    <t_end>
ENDSTATIC
```

### `THERMALSTATIC`
```
THERMALSTATIC
  EPSILON <tol>
  TIME    <t_end>
ENDTHERMALSTATIC
```

### `THERMALDYNAMIC`
```
THERMALDYNAMIC
  EPSILON    <tol>
  INTEGRATOR <BDF-1|BDF-2|...>
  TIME       <t0> <tf> <dt>
ENDTHERMALDYNAMIC
```

### `THERMOMECHANICALDYNAMIC`
```
THERMOMECHANICALDYNAMIC
  EPSILON    <tol>
  INTEGRATOR <type>
  TIME       <t0> <tf> <dt>
ENDTHERMOMECHANICALDYNAMIC
```

### `DYNAMIC`
```
DYNAMIC
  EPSILON    <tol>
  INTEGRATOR NEWMARK         <beta> <gamma>
  INTEGRATOR NEWMARK-ALPHA   <beta> <gamma>
  INTEGRATOR HHT-SIMPLE      <alpha>
  INTEGRATOR HHT-GENERALIZED <alpha>
  TIME       <t0> <tf> <dt>
ENDDYNAMIC
```

---

## Complete Minimal Example

```
// Thermal transient analysis of a two-triangle FEM mesh
TITLE 2triangfem

DIMENSION 2

MATERIALS
  THERMAL 1 1348 224 0 1750
ENDMATERIALS

SYSTEM 2triangles
  FLEXBODIES
    FEMESH triangle
      FORMULATION THERMAL
      MESH
      1 1 3.5
      TRIANGLES 2triangles.dat
    ENDFEMESH
  ENDFLEXBODIES

  JOINTS
  ENDJOINTS

  LOADS
    THERMALFLUENCE triangle.2 5E6
  ENDLOADS
ENDSYSTEM

ANALYSIS
  THERMALDYNAMIC
    EPSILON    1E-6
    INTEGRATOR BDF-1
    TIME       0.0 10E-3 1E-3
  ENDTHERMALDYNAMIC
ENDANALYSIS
```
