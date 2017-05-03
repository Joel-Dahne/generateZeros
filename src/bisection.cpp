/*
  Copyright (C) 2017  Joel Dahne
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//C-XSC libraries
#include "interval.hpp"
#include "complex.hpp"
#include "cinterval.hpp"
#include "civector.hpp"
#include "cimatrix.hpp"

//Private libraries
#include "citaylor.hpp"
#include "list.h"
#include "function.hpp"

//Standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <utility>

//OPENMP
#include "omp.h"

using namespace cxsc;
using namespace taylor;
using namespace std;


const interval PI(Pi());       // An enclosure of pi.
const interval PI2(sqr(Pi()));       // An enclosure of pi.
const cinterval I(interval(0), interval(1)); //An enclosure of i

int main(int argc, char * argv[]) {
  
  civector domain(2);
  domain[1] = cinterval(interval(atof(argv[1]), atof(argv[2])),
                        interval(atof(argv[3]), atof(argv[4])));
  domain[2] = cinterval(interval(atof(argv[5]), atof(argv[6])),
                        interval(atof(argv[7]), atof(argv[8])));
  
  
  List<civector> * currentWorkList;
  currentWorkList = new List<civector>;

  *currentWorkList+=domain;
  
  return 0;
}
