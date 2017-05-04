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

//****************************************
//Handling the jacobians
//****************************************

//Calculate the jacobian
cimatrix jacobian(const civector &domain, bool &ok, const interval &parameter) {
  cimatrix J(2, 2);

  citaylor z2(1, 0);
  
  //Compute df1/dz1 and df2/dz1
  citaylor z1 = citaylor(1, domain[1]);
  z2 = domain[2];
  J[1][1] = get_j_derive(function1(z1, z2, ok, parameter), 1);
  J[2][1] = get_j_derive(function2(z1, z2, ok, parameter), 1);

  //Compute df1/dz2 and df2/dz2
  z1 = domain[1];
  z2 = citaylor(1, domain[2]);
  J[1][2] = get_j_derive(function1(z1, z2, ok, parameter), 1);
  J[2][2] = get_j_derive(function2(z1, z2, ok, parameter), 1);

  return J;
}

//Computes the determinant of a 2x2 matrix
cinterval det(const cimatrix &M) {
  return M[1][1]*M[2][2] - M[1][2]*M[2][1];
}

//Calculate the inverse
cimatrix inverse(const cimatrix &M, bool &ok) {
  cinterval determinant = det(M);
  
  if (0 <= determinant) {
    ok = false;
    return M;
  }

  cimatrix inverse(2, 2);

  inverse[1][1] = M[2][2];
  inverse[1][2] = -M[1][2];
  inverse[2][1] = -M[2][1];
  inverse[2][2] = M[1][1];

  inverse /= determinant;  

  return inverse;
}

//****************************************
//End of Handling the Jacobians
//****************************************

//****************************************
//Handling civectors in different ways
//****************************************

void splitDomain(const civector &domain, civector &D1, civector &D2) {
  //Find the widest interval
  cvector widths = diam(domain);
  real width(0);
  int index;
  int part;

  for (int i = 1; i <= 2; ++i) {
    if (max(Re(widths[i]), Im(widths[i])) > width) {
      index = i;
      if (Re(widths[i]) > Im(widths[i])) {
        width = Re(widths[i]);
        part = 0;
      }
      else {
        width = Im(widths[i]);
        part = 1;
      }
    }
  }

  //Split the domain along the widest interval
  D1 = domain;
  D2 = domain;

  if (part == 0) {
    interval re = Re(domain[index]);
    D1[index] = cinterval(interval(Inf(re), Mid(re)), Im(domain[index]));
    D2[index] = cinterval(interval(Mid(re), Sup(re)), Im(domain[index]));
  } else {
    interval im = Im(domain[index]);
    D1[index] = cinterval(Re(domain[index]), interval(Inf(im), Mid(im)));
    D2[index] = cinterval(Re(domain[index]), interval(Mid(im), Sup(im)));
  }
}

//Determine if two intervals are disjoint
bool disjoint(interval &x1, interval &x2) {
  return Sup(x1) < Inf(x2) || Inf(x1) > Sup(x2); 
}

//Determines if to complex interval vector are disjoint
bool disjoint(civector &z1, civector &z2) {
  return disjoint(Re(z1[1]), Re(z2[1])) || disjoint(Im(z1[1]), Im(z2[1])) ||
    disjoint(Re(z1[2]), Re(z2[2])) || disjoint(Im(z1[2]), Im(z2[2]));
}

//Compute the mid point of a domain
civector Mid(civector &domain) {
  civector midPoint(2);

  for (int i = 1; i <= 2; i++) {
    midPoint[i] = cinterval(mid(domain[i]));
  }
  
  return midPoint;
}

//****************************************
//End of Handling civectors in different ways
//****************************************

//Check if the enclosure of the function on the domain contains a zero
bool zeroInDomain(civector &domain, bool &ok, interval parameter) {
  citaylor z1, z2;
  cinterval f1, f2;

  z1 = citaylor(0, domain[1]);
  z2 = citaylor(0, domain[2]);
  f1 = get_j_coef(function1(z1, z2, ok, parameter), 0);
  f2 = get_j_coef(function2(z1, z2, ok, parameter), 0);

  return (0<=f1 && 0<=f2 && ok);
}


civector N(civector &domain, bool &ok, interval parameter) {
  civector mid = Mid(domain);

  return mid - inverse(jacobian(domain, ok, parameter), ok)*function(mid, ok, parameter);  
}

civector validate(civector &domain, bool &isZero, bool &ok, interval &parameter) {

  civector newDomain;
  bool succCalc = true;
  isZero = false;
  ok = false;

  for (int i = 0; i < 10; i++) {
    newDomain = N(domain, succCalc, parameter);
    
    if (!succCalc) {
      //Failed calculations
      return domain;
    } else if (newDomain[1] <= domain[1] && newDomain[2] <= domain[2]) {
      //Is zero
      ok = true;
      isZero = true;
    } else if (disjoint(domain, newDomain)) {
      //Is not zero
      ok = true;
      return domain;
    }
    
    domain = newDomain & domain;
  }

  return domain;  
}

int main(int argc, char * argv[]) {

  string usage = "Usage: ./bisection [OPTION]... -- [DOMAIN] \n\
Find all zeros in a domain by using a combination of bisection\n\
and the Newton interval method\n\n\
Options are:\n\
  -p <value> set the parameter for the function to this\n\
  -v verbose - print more information\n                               \
Domain should be given after the options as:\n\
inf(real1) sup(real1) inf(im1) sup(im1) inf(real2) sup(real2) inf(im2) sup(im2)\n\n\
Example: ./bisection -v -- -1 1 -2 2 -3 3 -4 4\n\
This will use verbose output and find all zeros in the domain:\n\
[-1, 1] + i[-2, 2] x [-3, 3] + i[-4, 4]";

  //Read arguments
  //Read options
  //Set parameter for the function
  int parameterSet = 0;
  interval parameter = interval(0);

  //Verbose output
  int verbose = 0;

  //Do not print errors
  opterr = 0;
  int c;

  while ((c = getopt (argc, argv, "p:v")) != -1)
    switch (c)
      {
      case 'p':
        parameterSet = 1;
        parameter = interval(atof(optarg));
        break;
      case 'v':
        verbose = 1;
        break;
      case '?':
        if (optopt == 'p')
          cerr <<"Option -" << char(optopt) << " requires an argument.\n\n" << endl;
        else if (isprint (optopt))
          cerr << "Unknown option `-" << char(optopt) << "'.\n\n" << endl;
        else
          cerr << "Unknown option character `" << char(optopt) << "'.\n\n" << endl;
        cerr << usage << endl;
        exit(0);
      default:
        abort ();
      }

  //Read domain
  if (argc - optind < 8) {
    fprintf (stderr, "To few arguments, cannot read domain.\n\n");
    cerr << usage << endl;;
    exit(0);
  }
  
  civector domain(2);
  domain[1] = cinterval(interval(atof(argv[1]), atof(argv[2])),
                        interval(atof(argv[3]), atof(argv[4])));
  domain[2] = cinterval(interval(atof(argv[5]), atof(argv[6])),
                        interval(atof(argv[7]), atof(argv[8])));
  
  interval parameter = interval(0);
  
  List<civector> * currentWorkList;
  currentWorkList = new List<civector>;
  *currentWorkList+=domain;

  List<civector> * nextWorkList;
  nextWorkList = new List<civector>;

  bool done = false;

  // Permanent counters
  int step = 0;
  int zerosFound = 0;

  // Counters for each cycle
  int partsLeft = 1;
  int partsDone = 0;
  int partsFailed = 0;
  int partsDiscardedEnclosure = 0;
  int partsDiscardedNewton = 0;
  
  civector currentDomain, newDomain1, newDomain2, zero;
  cimatrix jac;
  
  bool ok, isZero, partDone, stepDone(false);
  bool partFailed, partDiscardedEnclosure, partDiscardedNewton;

  
  while(!done) {

    if (!IsEmpty(*currentWorkList)) {
      currentDomain = Pop(*currentWorkList);
      stepDone = false;
    } else {
      stepDone = true;
    }

    while (!stepDone) {
      ok = true;
      isZero = false;
      partDone = false;
      partFailed = false;
      partDiscardedEnclosure = false;
      partDiscardedNewton = false;
      
      if (zeroInDomain(currentDomain, ok, parameter)) {
        jac = jacobian(currentDomain, ok, parameter);

        if ((!(0 <= det(jac))) && ok) {
          zero = validate(currentDomain, isZero, ok, parameter);
          
          if (isZero && ok) {
            partDone = true;
            cout << "Zero: " << zero << endl;
          } else if (!isZero && ok) {
            partDiscardedNewton = true;
            partDone = true;
          }
        }
      } else {
        partDone = true;
        partDiscardedEnclosure = true;
      }

      if (!partDone) {
        splitDomain(currentDomain, newDomain1, newDomain2);
      }

#pragma omp critical
      {
        // Update counters
        if (partDiscardedEnclosure) {
          partsDiscardedEnclosure+=1;
        }
        if (partDiscardedNewton) {
          partsDiscardedNewton+=1;
        }
        if (partFailed) {
          partsFailed+=1;
        }
        if (isZero) {
          zerosFound+=1;
        }
        
        
        if (!partDone) {
          *nextWorkList+=newDomain1;
          *nextWorkList+=newDomain2;
          partsLeft+=2;
        } else {
          partsDone+=1;
        }
        
        if (!IsEmpty(*currentWorkList)) {
          currentDomain = Pop(*currentWorkList);
        } else {
          stepDone = true;
        }
      }
    }

    //Handle going to the next step

    cout << "Step " << step << endl;
    cout << "Parts done: " << partsDone << ", parts left: " << partsLeft << endl;
    cout << "Parts discarded from enclosure: " << partsDiscardedEnclosure << endl;
    cout << "Parts discarded from Newton " << partsDiscardedNewton << endl;
    cout << "Parts failed " << partsFailed << endl;
    
    //If not done set the next work list to the current one
    if (!IsEmpty(*nextWorkList) && !done) {
      delete currentWorkList;
      currentWorkList = nextWorkList;
      nextWorkList = new List<civector>;
    } else {
      done = true;
    }

    // Reset and update counter
    step+=1;
    partsDone = 0;
    partsLeft = 0;
    partsFailed = 0;
    partsDiscardedEnclosure = 0;
    partsDiscardedNewton = 0;  
  }

  cout << "Zeros found: " << zerosFound << endl;
  return 0;
}
