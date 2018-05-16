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

#include <iostream>

using namespace cxsc;
using namespace taylor;
using namespace std;

const interval PI(Pi());       // An enclosure of pi.
const interval PI2(sqr(Pi()));       // An enclosure of pi.
const cinterval I(interval(0), interval(1)); //An enclosure of i

//**********************************
//Implementation of the function analyzed
//*********************************

//Used in example 3
const interval F0 = interval(0.0534);
const interval w = interval(0.057);
const interval np = interval(3);
const interval phi = PI/2;

const interval wp = w/np;
const interval w0 = w;
const interval w1 = w + wp;
const interval w2 = w - wp;

const interval E0 = F0/2;
const interval E1 = -E0/2;
const interval E2 = E1;

const interval A0 = E0/w0;
const interval A1 = E1/w1;
const interval A2 = E2/w2;

const interval alpha0 = A0/w0;
const interval alpha1 = A1/w1;
const interval alpha2 = A2/w2;

const interval Ip = interval(0.5145);
const interval omega = interval(1);

// For example 6
citaylor A_simple(citaylor z) {
  return A0*cos(z);
}

citaylor ks_simple(citaylor z1, citaylor z2) {
  return A0*(sin(z2 + z1) - sin(z1))/z2;
}

citaylor A(citaylor z) {
  return A0*cos(np*z) + A1*cos((np + 1)*z) + A2*cos((np - 1)*z);
}

citaylor alpha(citaylor z) {
  return alpha0*sin(np*z) + alpha1*sin((np + 1)*z) + alpha2*sin((np - 1)*z);
}

citaylor ks(citaylor z1, citaylor z2) {
  return (alpha(z1) - alpha(z1 - z2))/z2;
}
//End of stuff used in example 3

// Function
void function(citaylor &f1, citaylor &f2, citaylor &z1, citaylor &z2, bool &ok,
              interval p) {
  /*
  // Example 1 - The identity function on C^2
  f1 = z1;
  f2 = z2;
  return;
  */
  /*
  // Example 2 - A polynomial function
  citaylor sqrz1 = sqr(z1);
  citaylor sqrz2 = sqr(z2);
  citaylor quadz2 = sqr(sqrz2);
  f1 = 4e-5*sqr(sqrz1)*z1*sqr(z2) + 2e-3*z1*quadz2 +
  2*sqrz1*z2 - z2 + real(0.75);
  f2 = 3e-4*z1*quadz2 - 7e-6*z1*sqrz1 + 2*z1*sqrz2 - z1 + real(0.75);
  return;
  */
  /*
  // Example 3 - The function f(z1, z2) = (sin(z1) + (z1)^2 +
  // e^{z2} - cos(2(z2)), cos(z1) + (z2)^3 + e^{2(z2)} - 2)
  f1 = sin(z1) + sqr(z1) + exp(z2) - cos(2*z2);
  f2 = cos(z1) + sqr(z2)*z2 + exp(2*z2) - 2;
  return;
  */
  /*
  // Example 3.5 - 1d real saddle point problem
    f1 = sqr(p + (A0*cos(np*z1) + A1*cos((np + 1)*z1) + A2*cos((np - 1)*z1))) + 2*Ip;
  //f1 = sqr(p + A0*cos(z1)) + 2*Ip;
  f2 = z2;
  return;
  */
  /*
  // Example 4 - Simplified 2d real saddle point problem
  citaylor sqrz2 = sqr(z2);
  citaylor sinz2 = sin(z2);
  citaylor sinz1z2 = sin(z1 - z2);
  citaylor z2cosz1z2 = z2*cos(z1 - z2);
  citaylor part = sinz1z2 - sinz2;
  f1 = real(2.25)*sqr(sinz1z2 - sinz2 + z2cosz1z2*z2) + sqrz2;
  f2 = real(2.25)*(sqr(cos(z1)*z2 + part) - sqr(part + z2cosz1z2*z2))
  - 5*sqrz2;
  return;
  */
  /*
  // Example 6 - 2d HHG - simple saddle point problem
  citaylor ksz1z2 = ks_simple(z1, z2);
  f1 = sqr(ksz1z2 + A_simple(z1)) + 2*Ip;
  f2 = sqr(ksz1z2 + A_simple(z2 + z1)) + 2*(Ip - omega);
  return;
  */

  // Example 6.5 - 2d HHG - saddle point problem
  citaylor ksz1z2 = ks(z1, z2);
  f1 = sqr(ksz1z2 + A(z1 - z2)) + 2*Ip;
  f2 = sqr(ksz1z2 + A(z1)) + 2*(Ip - omega);
  return;

  /*
  // Example 7 - 2d ATI - simple saddle point problem
  citaylor ksz1z2 = ks_simple(z1, z2);
  citaylor Az2z1 = A_simple(z2 + z1);
  f1 = sqr(ksz1z2 + A_simple(z1)) + 2*Ip;
  f2 = sqr(p + Az2z1) - sqr(ksz1z2 + Az2z1);
  return;
  */
  /*
  // Example 7 - 2d ATI - saddle point problem
  citaylor ksz1z2 = ks(z1, z2);
  citaylor Az2z1 = A(z2 + z1);
  f1 = sqr(ksz1z2 + A(z1)) + 2*Ip;
  f2 = sqr(p + Az2z1) - sqr(ksz1z2 + Az2z1);
  return;
  */
  /*
  // Example 8 - 2d real saddle point problem
  citaylor ksz1z2 = ks(z1, z2, ok);
  citaylor Az2 = A(z2);
  f1 = sqr(ksz1z2 + A(z1)) + 2*Ip;
  f2 = sqr(p + Az2) - sqr(ksz1z2 + Az2);
  return;
  */
  /*
  // Example 9 - The function f(z1, z2) = (z1^50 + z1^12 + -
  // sin(20z1)cos(12z1) - 1, z2)
  citaylor z1_2 = sqr(z1);
  citaylor z1_4 = sqr(z1_2);
  citaylor z1_8 = sqr(z1_4);
  citaylor z1_16 = sqr(z1_8);
  citaylor z1_32 = sqr(z1_16);
  f1 = z1_32*z1_16*z1_2 + z1_8*z1_4 - 5*sin(20*z1)*cos(12*z1) - 1;
  f2 = z2;
  */

  real Ip = 0.5;
  real omega = 1;
  real A0 = 1;
  citaylor pst = A0*(cos(z1) - cos(z1 - z2))/z2;
  f1 = pst + A0*sin(z1 - z2) + I*sqrt(2*Ip);
  f2 = pst + A0*sin(z1) + I*sqrt(omega - 2*Ip);

}

//*********************************
//These functions are used to make calling the above functions in
//different ways easier, they normally do not need to be changed.
//*********************************

//Get an interval enclosure of the function on the domain
civector intervalFunction(civector &domain, bool &ok,
                          interval p) {
  civector f(2);

  citaylor z1(0, 0), z2(0, 0);
  citaylor f1, f2;
  z1 = domain[1];
  z2 = domain[2];
  function(f1, f2, z1, z2, ok, p);
  f[1] = get_j_derive(f1, 0);
  f[2] = get_j_derive(f2, 0);

  return f;
}

//Evaluate the function in a point
cvector midFunction(cvector &z, bool &ok,
                    interval p) {
  cvector f(2);

  citaylor z1(0, 0), z2(0, 0);
  citaylor f1, f2;
  z1 = z[1];
  z2 = z[2];
  function(f1, f2, z1, z2, ok, p);
  f[1] = mid(get_j_derive(f1, 0));
  f[2] = mid(get_j_derive(f2, 0));

  return f;
}
