# MkniX

TODO: Write a project description

## Installation

MkniX is built using cmake (version 2.8 or greater).

Build using (from the MkniX root directory):

```
cmake -B<build directory> -H. -DCMAKE_BUILD_TYPE=<Debug|Release>
```

i.e.

```
cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
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

TODO: Write some usage instructions

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
