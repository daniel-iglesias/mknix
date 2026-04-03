# MkniX

MkniX is a simulation library for nonlinear thermo-mechanical multibody systems using FEM and mesh-free methods. It includes as well a standalone executable for forward simulations, and is used in other projects for such as [WHAM](https://github.com/daniel-iglesias/wham) and [ALICIA](https://github.com/daniel-iglesias/alicia) for real-time and inverse thermal operational monitoring, respectively. 

For details on the formulation, applications and benchmarking, you can consult the author's [PhD Thesis [1]](https://oa.upm.es/39121/1/Daniel_Iglesias_Ibanez.pdf). Please use any of the references bellow for citing this work.

## Installation

MkniX is built using cmake (version 2.8 or greater).

Build using (from the MkniX root directory):

```
cmake -B<build directory> -H. -DCMAKE_BUILD_TYPE=<Debug|Release>
```

i.e.

```
cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Debug
```

## Tests

MkniX tests can be run (after building) using the following:

```
cd <build directory>/tests
ctest
```

i.e.

```
cd build/tests
ctest
```

## Usage

Run MkniX with an input file as the sole argument:

```
mknix <input_file.mknix>
```

For a full description of the input file format and all available commands (materials, rigid/flexible bodies, joints, loads, motion, signals, and analysis types), see [doc/input_manual.md](doc/input_manual.md).

## References
[1] Iglesias, Daniel. On the application of meshfree methods to the nonlinear dynamics of multibody systems, 2017. [doi: 10.20868/upm.thesis.39121](https://doi.org/10.20868/upm.thesis.39121)

[2] Iglesias, D., García Orden, J.C. Galerkin meshfree methods applied to the nonlinear dynamics of flexible multibody systems. Multibody Syst Dyn 25, 203–224 (2011). [doi: 10.1007/s11044-010-9224-9](https://doi.org/10.1007/s11044-010-9224-9).

[3] D. Iglesias, J. C. García Orden, B. Brañas, J.M. Carmona, J. Molla. Application of Galerkin meshfree methods to nonlinear thermo-mechanical simulation
of solids under extremely high pulsed loading, Fusion Engineering and Design, Vol.
88(9-10), 2744–2747, 2013. [doi: 10.1016/j.fusengdes.2013.02.158](https://doi.org/10.1016/j.fusengdes.2013.02.158)

## License

Copyright (C) 2015 by Daniel Iglesias

MkniX is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

MkniX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.
