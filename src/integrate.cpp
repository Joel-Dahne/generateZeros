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

//Return the jacobian
cimatrix jacobian(const citaylor &D1f1, const citaylor &D2f1,
                  const citaylor &D1f2, const citaylor &D2f2) {
  cimatrix J(2, 2);
  
  J[1][1] = get_j_derive(D1f1, 1);
  J[1][2] = get_j_derive(D2f1, 1);
  J[2][1] = get_j_derive(D1f2, 1);
  J[2][2] = get_j_derive(D2f2, 1);
  
  return J;
}

//The jacobian with the kth column replaced with
//[f_1, ..., f_n]^T
cimatrix jacobian_k(const citaylor &D1f1, const citaylor &D2f1,
                    const citaylor &D1f2, const citaylor &D2f2,
                    const int &side) {
  cimatrix J(2, 2);

  if (side/4 == 0) {
    J[1][1] = get_j_derive(D1f1, 0);
    J[1][2] = get_j_derive(D2f1, 0);
    J[2][1] = get_j_derive(D1f2, 1);
    J[2][2] = get_j_derive(D2f2, 1);
  } else {
    J[1][1] = get_j_derive(D1f1, 1);
    J[1][2] = get_j_derive(D2f1, 1);
    J[2][1] = get_j_derive(D1f2, 0);
    J[2][2] = get_j_derive(D2f2, 0);
  }
    
  return J;
}

//The jacobian with the kth column replaced with
//[f_1, ..., f_n]^T
cimatrix jacobian_k(const cimatrix &jacobian, const civector &function,
                    const int &side) {
  cimatrix J(2, 2);

  if (side/4 == 0) {
    J[1][1] = function[1];
    J[1][2] = function[2];
    J[2][1] = jacobian[2][1];
    J[2][2] = jacobian[2][2];
  } else {
    J[1][1] = jacobian[1][1];
    J[1][2] = jacobian[1][2];
    J[2][1] = function[1];
    J[2][2] = function[2];
  }
    
  return J;
}

//Computes the determinant of a 2x2 matrix
inline cinterval det(cimatrix &M) {
  return M[1][1]*M[2][2] - M[1][2]*M[2][1];
}

//****************************************
//End of Handling the Jacobians
//****************************************


//****************************************
//Handling civectors in different ways
//****************************************

//Splits the domain into two parts along the widest interval
//The new intervals are given in D1 and D2
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

//Take a part of the integral, split the domain and adjust the
//tolerance
void splitIntegrationPart(const pair<civector, pair<real, int> > &part,
                          pair<civector, pair<real, int> > &p1,
                          pair<civector, pair<real, int> > &p2) {
  civector D1, D2;
  real tol = part.second.first/2;
  int side = part.second.second;

  splitDomain(part.first, D1, D2);
  
  p1 = make_pair(D1, make_pair(tol, side));
  p2 = make_pair(D2, make_pair(tol, side));
}

//Computes the square of the length of the vector
//|z_1|^2 + ... + |z_n|^2
inline interval length2(civector &V) {
  return sqr(abs(V[1])) + sqr(abs(V[2]));
}

//Calculates the volume of the domain
//Assumes that the side corresponding to the int side
//has diameter 0
interval volume(const civector &V, const int &side) {
  
  cvector diameters = diam(V);
  interval volume(1);

  if (side < 2)
    volume *= Im(diameters[1])*Re(diameters[2])*Im(diameters[2]);
  else if (side < 4)
    volume *= Re(diameters[1])*Re(diameters[2])*Im(diameters[2]);
  else if (side < 6)
    volume *= Im(diameters[2])*Re(diameters[1])*Im(diameters[1]);
  else
    volume *= Re(diameters[2])*Re(diameters[1])*Im(diameters[1]);

  return volume;
}

//****************************************
//End of Handling civectors in different ways
//****************************************

//The integrand used in integration
cinterval integrand(civector &domain, const int &side,
                    bool &ok, interval &parameter) {

  citaylor z1, z2, constz1(1, 0), constz2(1, 0);
  
  z1 = citaylor(1, domain[1]);
  z2 = citaylor(1, domain[2]);
  constz1 = domain[1];
  constz2 = domain[2];

  citaylor D1f1, D2f1, D1f2, D2f2;

  D1f1 = function1(z1, constz2, ok, parameter);
  if (!ok) {
    ok = false;
    return cinterval(0);
  }
  D2f1 = function1(constz1, z2, ok, parameter);
  if (!ok) {
    ok = false;
    return cinterval(0);
  }
  D1f2 = function2(z1, constz2, ok, parameter);
  if (!ok) {
    ok = false;
    return cinterval(0);
  }
  D2f2 = function2(constz1, z2, ok, parameter);
  if (!ok) {
    ok = false;
    return cinterval(0);
  }
  
  civector enclosure(2);
  enclosure[1] = get_j_derive(D1f1, 0);
  enclosure[2] = get_j_derive(D1f2, 0);
    
  interval fLength2 = length2(enclosure);
   
  if (0 <= fLength2) {
    ok = false;
    return cinterval(0);
  }

  //cimatrix J = jacobian(domain, ok, parameter);
  cimatrix J = jacobian(D1f1, D2f1, D1f2, D2f2);
  
  //cimatrix J_k = jacobian_k(J, enclosure, side);
  cimatrix J_k = jacobian_k(D1f1, D2f1, D1f2, D2f2, side);

  return det(J)*conj(det(J_k))/sqr(fLength2);
}

int main(int argc, char * argv[]) {

  string usage = "Usage: ./integrate [OPTION]... -- [DOMAIN] \n\
Calculate the number of zeros for a function in a domain by\n\
calculating the logarithmic integral\n\n\
Options are:\n\
  -p <value> set the parameter for the function to this\n\
  -t <value> set the tolerance to use, default 1\n\
  -s 0-7 integrate only a specific side, default is all sides\n\
  -v verbose - print more information\n\
Domain should be given after the options as:\n\
inf(real1) sup(real1) inf(im1) sup(im1) inf(real2) sup(real2) inf(im2) sup(im2)\n\n\
Example: ./integrate -v -- -1 1 -2 2 -3 3 -4 4\n\
This will use verbose output and calculate the integral for the domain:\n\
[-1, 1] + i[-2, 2] x [-3, 3] + i[-4, 4]";

  //Read arguments

  //Read options
  //Set parameter for the function
  int parameterSet = 0;
  interval parameter = interval(0);

  //Tolerance to use
  real tol(1);

  //Sides to integrate
  vector<int> sides;
  
  //Verbose output
  int verbose = 0;

  //Do not print errors
  opterr = 0;
  int c;

  while ((c = getopt (argc, argv, "p:t:s:v")) != -1)
    switch (c)
      {
      case 'p':
        parameterSet = 1;
        parameter = interval(atof(optarg));
        break;
      case 't':
        tol = atof(optarg);
        break;
      case 's':
        if (atoi(optarg) >= 0 && atoi(optarg) <=7) {
          sides.push_back(atoi(optarg));
        } else {
          cerr << "Option -s requres an integer between 0 and 7 as argument" << endl;
          cerr << usage << endl;
          exit(0);
        }
        break;
      case 'v':
        verbose = 1;
        break;
      case '?':
        if (optopt == 'p' || optopt == 't' || optopt == 's')
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


  //Set up the program

  //Create a list of the sides to integrate. Each element is a pair,
  //the first part is the side and the second part corresponds to the
  //location of the side
  if (sides.size() == 0) {
    for (int i = 0; i < 8; ++i) {
      sides.push_back(i);
    }
  }
  
  List<pair<civector, pair<real, int> > > currentWorkList;

  for (int j = 0; j < sides.size(); ++j) {
    int i = sides[j];
    civector side = domain;

    if (i == 0) {
      Re(side[1]) = interval(InfRe(side[1]));
    } else if (i == 1) {
      Re(side[1]) = interval(SupRe(side[1]));
    } else if (i == 2) {
      Im(side[1]) = interval(InfIm(side[1]));
    } else if (i == 3) {
      Im(side[1]) = interval(SupIm(side[1]));
    } else if (i == 4) {
      Re(side[2]) = interval(InfRe(side[2]));
    } else if (i == 5) {
      Re(side[2]) = interval(SupRe(side[2]));
    } else if (i == 6) {
      Im(side[2]) = interval(InfIm(side[2]));
    } else if (i == 7) {
      Im(side[2]) = interval(SupIm(side[2]));
    }
    
    currentWorkList += make_pair(side, make_pair(tol/(int)sides.size(), i));
  }
  
  //For storing the result of the integral
  cinterval integral(0);
  cinterval tmpIntegral(0);
  
  cinterval sideIntegral[8];
  cinterval tmpSideIntegral[8];
  for (int i = 0; i < 8; ++i) {
    sideIntegral[i] = cinterval(interval(0), interval(0));
    tmpSideIntegral[i] = cinterval(interval(0), interval(0));

  }

  //For storing number of calls
  int steps = 0;

  //Constant to multiply the integral width
  cinterval coefs[8];
  for (int i = 0; i < 8; ++i) {
    coefs[i] = cinterval(interval(0), 1/(-2*PI2));
    if (i == 0)
      coefs[i] *= -I;
    else if (i == 1)
      coefs[i] *= I;
    else if (i == 2)
      coefs[i] *= 1;
    else if (i == 3)
      coefs[i] *= -1;
    else if (i == 4)
      coefs[i] *= -I;
    else if (i == 5)
      coefs[i] *= I;
    else if (i == 6)
      coefs[i] *= 1;
    else
      coefs[i] *= -1;
  }

  //Initiate the citaylor class
  citaylor tmp;

  //Split the current integration parts a few times to avoid a problem
  //with parallization
  pair<civector, pair<real, int> > p1, p2;
  for (int i = 0; i < 24; ++i) {
    splitIntegrationPart(Pop(currentWorkList), p1, p2);
    currentWorkList += p1;
    currentWorkList += p2;
  }
  
#pragma omp parallel
  {
    //initiate variables
    //For storing integral part
    pair<civector, pair<real, int> > currentPart, newPart1, newPart2;
    //For storing info from currentPart
    civector currentDomain;
    real tol;
    int side;
    //For storint integral evaluations
    cinterval integrandEnclosure;
    interval domainVolume;
    cinterval currentIntegral;
    //For checking the result
    bool ok(true);
    //For storint the width of the integral
    real width;
    //For cheking when done
    bool done(false);
    //For checking if the integration succeded and satisfies the
    //tolerance
    bool partDone(false);
    
#pragma omp critical
    {
      if (omp_get_thread_num() == 0 && verbose)
        cout << "Number of threads: " << omp_get_num_threads() << endl;
      
      if (!IsEmpty(currentWorkList)) {
        currentPart = Pop(currentWorkList);
      } else {
        done = true;
      }
    }

    while(!done) {

      ok = true;
      
      currentDomain = currentPart.first;
      tol = currentPart.second.first;
      side = currentPart.second.second;

      integrandEnclosure = integrand(currentDomain, side, ok, parameter);
      domainVolume = volume(currentDomain, side);

      currentIntegral = integrandEnclosure*domainVolume*coefs[side];

      width = diam(Re(currentIntegral));

      partDone = width <= tol && ok;

      if (!partDone) {
        splitIntegrationPart(currentPart, newPart1, newPart2);
      }
      
#pragma omp critical
      {
        ++steps;
        if (!partDone) {
          currentWorkList += newPart1;
          currentWorkList += newPart2;
        } else {
          integral += currentIntegral;
          sideIntegral[side] += currentIntegral;
        }

        if (!IsEmpty(currentWorkList)) {
          currentPart = Pop(currentWorkList);
        } else {
          done = true;
        }
      }
    }

    if (verbose) {
#pragma omp critical
      {
        cout << "Thread number " << omp_get_thread_num() << " done" << endl;
      }
    }
  }

  for (int i = 0; i < sides.size(); i++) {
    int j = sides[i];
    cout << "Integral side " << j << ": " << sideIntegral[j] << endl;
  }
  
  cout << "Total integral: " << integral << endl;

  if (verbose) {
    cout << "Number of steps used: " << steps << endl;
  }
  
  return 0;
}