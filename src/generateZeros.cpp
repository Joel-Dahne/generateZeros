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
#include "function.hpp"

//Standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <unistd.h>

using namespace cxsc;
using namespace taylor;
using namespace std;

//****************************************
//Handling the jacobians
//****************************************

//Calculate the jacobian
cmatrix jacobian(cvector &domain, bool &ok, interval &parameter) {
  cmatrix J(2, 2);
  citaylor f1, f2;

  citaylor z2(1, 0);

  //Compute df1/dz1 and df2/dz1
  citaylor z1 = citaylor(1, domain[1]);
  z2 = domain[2];
  function (f1, f2, z1, z2, ok, parameter);
  J[1][1] = mid(get_j_derive(f1, 1));
  J[2][1] = mid(get_j_derive(f2, 1));

  //Compute df1/dz2 and df2/dz2
  z1 = domain[1];
  z2 = citaylor(1, domain[2]);
  function (f1, f2, z1, z2, ok, parameter);
  J[1][2] = mid(get_j_derive(f1, 1));
  J[2][2] = mid(get_j_derive(f2, 1));

  return J;
}

//Computes the determinant of a 2x2 matrix
complex det(cmatrix &M) {
  return M[1][1]*M[2][2] - M[1][2]*M[2][1];
}

//Calculate the inverse
cmatrix inverse(cmatrix &M, bool &ok) {
  complex determinant = det(M);

  if (determinant == 0) {
    ok = false;
    return M;
  }

  cmatrix inverse(2, 2);

  inverse[1][1] = M[2][2];
  inverse[1][2] = -M[1][2];
  inverse[2][1] = -M[2][1];
  inverse[2][2] = M[1][1];

  inverse /= determinant;

  return inverse;
}

//Create a random point in the domain
cvector createRandomPoint(civector &domain) {

  real re1 = InfRe(domain[1]) + ((float) rand() / (float) RAND_MAX)*diam(Re(domain[1]));
  real im1 = InfIm(domain[1]) + ((float) rand() / (float) RAND_MAX)*diam(Im(domain[1]));
  real re2 = InfRe(domain[2]) + ((float) rand() / (float) RAND_MAX)*diam(Re(domain[2]));
  real im2 = InfIm(domain[2]) + ((float) rand() / (float) RAND_MAX)*diam(Im(domain[2]));

  cvector z(2);
  z[1] = complex(re1, im1);
  z[2] = complex(re2, im2);

  return z;
}

void newtonStep(cvector &z, bool &ok, interval parameter) {
  cvector mid = midFunction(z, ok, parameter);
  cmatrix jac = jacobian(z, ok, parameter);
  z = z - inverse(jac, ok)*mid;
}

void newtonSteps(cvector &z, civector &domain, int iterations, bool &ok,
                 interval parameter = interval(0)) {

  for (int i = 0; i < iterations; ++i) {

    newtonStep(z, ok, parameter);

    if (!ok || !(z[1] <= domain[1]) || !(z[2] <= domain[2])) {
      //Outside of domain
      ok = false;
      break;
    }

  }
}

//****************************************
//Handling the jacobians
//****************************************

//Calculate the jacobian
cimatrix jacobian(const civector &domain, bool &ok, const interval &parameter) {
  cimatrix J(2, 2);
  citaylor f1, f2;

  citaylor z2(1, 0);

  //Compute df1/dz1 and df2/dz1
  citaylor z1 = citaylor(1, domain[1]);
  z2 = domain[2];
  function(f1, f2, z1, z2, ok, parameter);
  J[1][1] = get_j_derive(f1, 1);
  J[2][1] = get_j_derive(f2, 1);

  //Compute df1/dz2 and df2/dz2
  z1 = domain[1];
  z2 = citaylor(1, domain[2]);
  function(f1, f2, z1, z2, ok, parameter);
  J[1][2] = get_j_derive(f1, 1);
  J[2][2] = get_j_derive(f2, 1);

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

civector Mid(civector &domain) {
  civector midPoint(2);

  for (int i = 1; i <= 2; i++) {
    midPoint[i] = cinterval(mid(domain[i]));
  }

  return midPoint;
}

//Enclose a point in a domain of width tol
civector enclose(cvector &z, real &tol) {
  civector domain(2);

  domain[1] = cinterval(interval(Re(z[1]) - tol/2, Re(z[1]) + tol/2),
                        interval(Im(z[1]) - tol/2, Im(z[1]) + tol/2));
  domain[2] = cinterval(interval(Re(z[2]) - tol/2, Re(z[2]) + tol/2),
                        interval(Im(z[2]) - tol/2, Im(z[2]) + tol/2));

  return domain;
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

civector N(civector &domain, bool &ok, interval parameter) {
  civector mid = Mid(domain);

  return mid - inverse(jacobian(domain, ok, parameter), ok)*intervalFunction(mid, ok, parameter);
}

civector validate(cvector &z, bool &valid, real &tol,
                  interval &parameter) {

  bool ok(true);

  civector domain = enclose(z, tol);

  civector newDomain = N(domain, ok, parameter);

  for (int i = 0; i < 10; i++) {
    if (!ok) {
      //Failed calculations
      valid = false;
      return domain;
    } else if (newDomain[1] <= domain[1] && newDomain[2] <= domain[2]) {
      //Is zero
      valid = true;
      return newDomain;
    } else if (disjoint(domain, newDomain)) {
      //Is not zero
      valid = false;
      return domain;
    } else {
      //Undecided
      domain = newDomain & domain;
      newDomain = N(domain, ok, parameter);
    }
  }

  valid = false;
  return domain;
}

bool isNewZero(civector &zero, vector<civector> &zeros) {

  for (int i = 0; i < zeros.size(); i++) {
    if (!disjoint(zero, zeros[i])) {
      return false;
    }
  }

  return true;
}

int main(int argc, char * argv[]) {

  string  usage = "Usage: ./generateZeros [OPTION]... -- [DOMAIN] \n\
Generates and validates zeros for a function and a specified domain\n\n\
Options are:\n\
  -n <value> stop when this many zeros are found\n\
  -t <value> stop after this amount of time\n\
  -s <value> stop after this many steps\n\
  -l <value> stop when the number of steps since last found zero equals this\n\
  -p <value> set the parameter for the function to this\n\
  -v verbose - print more information\n\
  -i print the intervals containing zeros instead of the mid point\n\
Domain should be given after the options as:\n\
inf(real1) sup(real1) inf(im1) sup(im1) inf(real2) sup(real2) inf(im2) sup(im2)\n\n\
Example: ./generateZeros -v -n 100 -- -1 1 -2 2 -3 3 -4 4\n\
This will use verbose output, stop when 100 zeros are found and use the domain:\n\
[-1, 1] + i[-2, 2] x [-3, 3] + i[-4, 4]";

  //Read arguments

  //Read options

  //Set stop criterias
  //Stop when the number of zeros reached this
  int stopOnZeros = 0;
  int maxZeros = 0;
  //Stop when the time ran reaches this
  int stopOnTime = 0;
  int maxTime = 0;
  //Stop when the number of steps reaches this
  int stopOnSteps = 0;
  int maxSteps = 0;
  //Stop when the number of steps since last found reaches this
  int stopOnStepsLast = 0;
  int maxStepsLast = 0;

  //Set parameter for the function
  int parameterSet = 0;
  interval parameter = interval(0);

  //Verbose output
  int verbose = 0;

  //Print intervals instead of mid points
  int printIntervals = 0;

  //Do not print errors
  opterr = 0;
  int c;

  while ((c = getopt (argc, argv, "n:t:s:l:p:vi")) != -1)
    switch (c)
      {
      case 'n':
        stopOnZeros = 1;
        maxZeros = atoi(optarg);
        break;
      case 't':
        stopOnTime = 1;
        maxTime = atoi(optarg);
        break;
      case 's':
        stopOnSteps = 1;
        maxSteps = atoi(optarg);
        break;
      case 'l':
        stopOnStepsLast = 1;
        maxStepsLast = atoi(optarg);
        break;
      case 'p':
        parameterSet = 1;
        parameter = interval(atof(optarg));
        break;
      case 'v':
        verbose = 1;
        break;
      case 'i':
        printIntervals = 1;
        break;
      case '?':
        if (optopt == 't' || optopt == 'n' || optopt == 'l' || optopt == 's' || optopt == 'p')
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
  domain[1] = cinterval(interval(atof(argv[optind]), atof(argv[optind + 1])),
                        interval(atof(argv[optind + 2]), atof(argv[optind + 3])));
  domain[2] = cinterval(interval(atof(argv[optind + 4]), atof(argv[optind + 5])),
                        interval(atof(argv[optind + 6]), atof(argv[optind + 7])));


  //Start the rest of the program

  //Set parameter for the newton methods
  int iterations = 15;
  real tol = 1e-5;

  //Initiate random seed
  srand(time(NULL));

  //Set up counters
  int found = 0;
  int steps = 0;
  int stepsLast = 0;
  int startTime = time(NULL);

  vector<civector> zeros;
  zeros.reserve(maxZeros);


  while (true) {

    //Check for stop criterias
    if (steps%100 == 0) {
      if (stopOnZeros && found >= maxZeros)
        break;
      if (stopOnTime && (time(NULL) - startTime) >= maxTime)
        break;
      if (stopOnSteps && steps >= maxSteps)
        break;
      if (stopOnStepsLast && stepsLast > maxStepsLast)
        break;
    }

    //Create a random point
    cvector z = createRandomPoint(domain);

    //Peform a number of newton steps
    bool ok(true);
    newtonSteps(z, domain, iterations, ok, parameter);

    if (ok) {

      //Try to validate the zero with newtons validated method
      bool valid(false);
      civector zi = validate(z, valid, tol, parameter);

      if (valid && isNewZero(zi, zeros)) {
        //If the validation succeds add it to the list of zeros
        found++;
        stepsLast = -1;
        zeros.push_back(zi);

        //Print the found zero
        if (printIntervals) {
          if (verbose) {
            cout << found << ": " << zi;
            cout << "Current number of Steps: " << steps << endl;;
          } else {
            cout << zi << endl;
          }
        } else {
          real x1 = mid(Re(zi[1]));
          real y1 = mid(Im(zi[1]));
          real x2 = mid(Re(zi[2]));
          real y2 = mid(Im(zi[2]));
          if (verbose) {
            cout << found << ": " <<  x1 << " " << y1 << " " << x2 << " " << y2;
            cout << " (" << steps << ")" << endl;
          } else {
            cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
          }
        }
      }
    }

    ++steps;
    ++stepsLast;
  }

  return 0;
}
