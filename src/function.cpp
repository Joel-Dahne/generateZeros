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

using namespace cxsc;
using namespace taylor;

const interval PI(Pi());       // An enclosure of pi.
const interval PI2(sqr(Pi()));       // An enclosure of pi.
const cinterval I(interval(0), interval(1)); //An enclosure of i

//**********************************
//Implementation of the function analyzed
//*********************************

//Used in example 3
const interval F0 = interval(0.0534);
const interval w = interval(0.057);
const interval np = interval(4);
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

// For example 6
citaylor A_simple(citaylor &z) {
  return A0*cos(z);
}

citaylor ks_simple_func(citaylor &z1, citaylor &z2) {
  return A0*(sin(z2) - sin(z1))/(z2 - z1);
}

civector ks_simple_i(cinterval z1i, cinterval z2i, int order) {
  citaylor z1 = citaylor(order, z1i);
  citaylor z2 = citaylor(order, z2i);
  return get_all_coef(A0*(sin(z2) - sin(z1))/(z2 - z1));
}

civector map_union(civector v1, civector v2) {
  for (int i = Lb(v1); i <= Ub(v1); i++) {
    v1[i] = v1[i] | v2[i];
  }
  return v1;
}

citaylor ks_simple(citaylor &z1, citaylor &z2) {
  citaylor z_diff = z2 - z1;
  if (0 <= get_j_coef(z_diff, 0)) {
    cinterval z1i = get_j_coef(z1, 0);
    cinterval z2i = get_j_coef(z2, 0);
    int order = get_order(z1);
    real epsilon = min(0.25, max(max(diam(Re(z1i)), diam(Im(z1i))),
                                 max(diam(Re(z2i)), diam(Im(z2i)))));

    cinterval square = z1i & z2i;
    civector result;

    // Handle the part inside of the square

    // Handle the line
    complex p = complex(InfRe(square), InfIm(square));
    cinterval eps_square = cinterval(interval(-epsilon, epsilon),
                                     interval(-epsilon, epsilon));

    // Initiate the result to something - CXSC cannot handle empty
    // intervals.
    result = get_all_coef(cos(citaylor(order, p + eps_square)));

    // To get an enclosure of the function on a small square centered
    // over the plane z1=z2 we use the following method. If the square
    // is centered at p and P is the hypercube obtained by moving
    // epsilon in every direction: NOT ACCURATE (BUT SOME
    // APPROXIMATION) A0 * cos(P) and derivative A0*sin(P)
    while ((Re(p) <= SupRe(square))) {
      while ((Im(p) <= SupIm(square))) {
        result = map_union(result,
                           get_all_coef(cos(citaylor(order, p + eps_square))));

        p = p + complex(0, epsilon);
      }

      p = SetIm(p, InfIm(square));
      p = p + complex(epsilon, 0);
    }

    // Handle the part in the square not on the line
    p = complex(InfRe(square), InfIm(square));
    cinterval z1_tmp, z2_tmp;
    while ((Re(p) < SupRe(square)) && (Im(p) < SupIm(square))) {
      // Up for real
      z1_tmp = cinterval(interval(Re(p) - epsilon, Re(p)),
                         Im(square));
      z2_tmp = cinterval(interval(min(Re(p) + epsilon, SupRe(square)),
                                  SupRe(square)),
                         Im(square));
      result = map_union(result, ks_simple_i(z1_tmp, z2_tmp, order));

      // Left for real
      z1_tmp = cinterval(interval(min(Re(p) + epsilon, SupRe(square)),
                                  SupRe(square)),
                         Im(square));
      z2_tmp = cinterval(interval(Re(p) - epsilon, Re(p)),
                         Im(square));
      result = map_union(result, ks_simple_i(z1_tmp, z2_tmp, order));

      // Up for imaginary
      z1_tmp = cinterval(Re(square),
                         interval(Im(p) - epsilon, Im(p)));
      z2_tmp = cinterval(Re(square),
                         interval(min(Im(p) + epsilon, SupIm(square)),
                                  SupIm(square)));
      result = map_union(result, ks_simple_i(z1_tmp, z2_tmp, order));

      // Left for imaginary
      z1_tmp = cinterval(Re(square),
                         interval(min(Im(p) + epsilon, SupIm(square)),
                                  SupIm(square)));
      z2_tmp = cinterval(Re(square),
                         interval(Im(p) - epsilon, Im(p)));
      result = map_union(result, ks_simple_i(z1_tmp, z2_tmp, order));

      p = p + complex(epsilon, epsilon);
    }

    //cout << "Result after rest in square: " << result << endl;

    // Handle the part outside of the square
    interval real1bl, real1tr, real2bl, real2tr;
    if (InfRe(z1i) < InfRe(z2i) - epsilon) {
      real1bl = interval(InfRe(z1i), InfRe(z2i) - epsilon);
      real2bl = Re(z2i);

      result = ks_simple_i(cinterval(real1bl, Im(z1i)),
                           cinterval(real2bl, Im(z2i)),
                           order);
    } else if (InfRe(z2i) < InfRe(z1i) - epsilon){
      real1bl = Re(z1i);
      real2bl = interval(InfRe(z2i), InfRe(z1i) - epsilon);

      result = ks_simple_i(cinterval(real1bl, Im(z1i)),
                           cinterval(real2bl, Im(z2i)),
                           order);
    }
    if (SupRe(z1i) > SupRe(z2i) + epsilon) {
      real1tr = interval(SupRe(z2i) + epsilon, SupRe(z1i));
      real2tr = Re(z2i);

      result = map_union(result,
                         ks_simple_i(cinterval(real1tr, Im(z1i)),
                                     cinterval(real2tr, Im(z2i)),
                                     order));
    } else if (SupRe(z2i) > SupRe(z1i) + epsilon){
      real1tr = Re(z1i);
      real2tr = interval(SupRe(z1i) + epsilon, SupRe(z2i));

      result = map_union(result,
                         ks_simple_i(cinterval(real1tr, Im(z1i)),
                                     cinterval(real2tr, Im(z2i)),
                                     order));
    }
    interval im1bl, im1tr, im2bl, im2tr;
    if (InfIm(z1i) < InfIm(z2i) - epsilon) {
      im1bl = interval(InfIm(z1i), InfIm(z2i) - epsilon);
      im2bl = Im(z2i);

      result = map_union(result,
                         ks_simple_i(cinterval(Re(z1i), im1bl),
                                     cinterval(Re(z2i), im2bl),
                                     order));
    } else if (InfIm(z2i) < InfIm(z1i) - epsilon){
      im1bl = Im(z1i);
      im2bl = interval(InfIm(z2i), InfIm(z1i) - epsilon);

      result = map_union(result,
                         ks_simple_i(cinterval(Re(z1i), im1bl),
                                     cinterval(Re(z2i), im2bl),
                                     order));
    }
    if (SupIm(z1i) > SupIm(z2i) + epsilon) {
      im1tr = interval(SupIm(z2i) + epsilon, SupIm(z1i));
      im2tr = Im(z2i);

      result = map_union(result,
                         ks_simple_i(cinterval(Re(z1i), im1tr),
                                     cinterval(Re(z2i), im2tr),
                                     order));
    } else if (SupIm(z2i) > SupIm(z1i) + epsilon){
      im1tr = Im(z1i);
      im2tr = interval(SupIm(z1i) + epsilon, SupIm(z2i));

      result = map_union(result,
                         ks_simple_i(cinterval(Re(z1i), im1tr),
                                     cinterval(Re(z2i), im2tr),
                                     order));
    }

    return citaylor(result);
  } else {
    return ks_simple_func(z1, z2);
  }
}

citaylor A(citaylor &z) {
  return A0*cos(w0*z + phi) + A1*cos(w1*z + phi) + A2*cos(w2*z + phi);
}

citaylor alpha(citaylor &z) {
  return alpha0*sin(w0*z + phi) + alpha1*sin(w1*z + phi) + alpha2*sin(w2*z + phi);
}

citaylor ks(citaylor &z1, citaylor &z2, bool &ok) {
  if (0 <= get_j_derive(z1, 0) - get_j_derive(z2, 0)) {
    ok = false;
    return z1;
  }
  citaylor zm = z1-z2;
  citaylor zp = z1+z2;
  citaylor part1 = 2*alpha0*cos(w0/2*zp + phi);
  citaylor part2 = 2*alpha1*cos(w1/2*zp + phi);
  citaylor part3 = 2*alpha2*cos(w2/2*zp + phi);

  if (0 <= get_j_coef(zm, 0)) {
    if (get_order(z1) == 0) {
      cinterval zm0 = get_j_coef(zm, 0);
      cinterval z;
      cinterval tmp1e, tmp2e, tmp3e;
      z = InfRe(zm0) + I*Im(zm0);
      tmp1e = sin(w0/2*z)/z;
      tmp2e = sin(w1/2*z)/z;
      tmp3e = sin(w2/2*z)/z;
      z = SupRe(zm0) + I*Im(zm0);
      tmp1e = tmp1e | sin(w0/2*z)/z;
      tmp2e = tmp2e | sin(w1/2*z)/z;
      tmp3e = tmp3e | sin(w2/2*z)/z;
      z = Re(zm0) + I*InfIm(zm0);
      tmp1e = tmp1e | sin(w0/2*z)/z;
      tmp2e = tmp2e | sin(w1/2*z)/z;
      tmp3e = tmp3e | sin(w2/2*z)/z;
      z = Re(zm0) + I*SupIm(zm0);
      tmp1e = tmp1e | sin(w0/2*z)/z;
      tmp2e = tmp2e | sin(w1/2*z)/z;
      tmp3e = tmp3e | sin(w2/2*z)/z;
      citaylor tmp1 = citaylor(0, tmp1e);
      citaylor tmp2 = citaylor(0, tmp2e);
      citaylor tmp3 = citaylor(0, tmp3e);
      part1 = part1*tmp1;
      part2 = part2*tmp2;
      part3 = part3*tmp3;
    } else {
      cinterval zm0 = get_j_coef(zm, 0);
      cinterval z;
      civector tmp1e(0, 1), tmp2e(0, 1), tmp3e(0, 1);
      z = InfRe(zm0) + I*Im(zm0);
      tmp1e[0] = sin(w0/2*z)/z;
      tmp2e[0] = sin(w1/2*z)/z;
      tmp3e[0] = sin(w2/2*z)/z;
      tmp1e[1] = w0/2*cos(w0/2*z)/z - sin(w0/2*z)/sqr(z);
      tmp2e[1] = w1/2*cos(w1/2*z)/z - sin(w1/2*z)/sqr(z);
      tmp3e[1] = w2/2*cos(w2/2*z)/z - sin(w2/2*z)/sqr(z);
      z = SupRe(zm0) + I*Im(zm0);
      tmp1e[0] = tmp1e[0] | sin(w0/2*z)/z;
      tmp2e[0] = tmp2e[0] | sin(w1/2*z)/z;
      tmp3e[0] = tmp3e[0] | sin(w2/2*z)/z;
      tmp1e[1] = tmp1e[1] | w0/2*cos(w0/2*z)/z - sin(w0/2*z)/sqr(z);
      tmp2e[1] = tmp2e[1] | w1/2*cos(w1/2*z)/z - sin(w1/2*z)/sqr(z);
      tmp3e[1] = tmp3e[1] | w2/2*cos(w2/2*z)/z - sin(w2/2*z)/sqr(z);
      z = Re(zm0) + I*InfIm(zm0);
      tmp1e[0] = tmp1e[0] | sin(w0/2*z)/z;
      tmp2e[0] = tmp2e[0] | sin(w1/2*z)/z;
      tmp3e[0] = tmp3e[0] | sin(w2/2*z)/z;
      tmp1e[1] = tmp1e[1] | w0/2*cos(w0/2*z)/z - sin(w0/2*z)/sqr(z);
      tmp2e[1] = tmp2e[1] | w1/2*cos(w1/2*z)/z - sin(w1/2*z)/sqr(z);
      tmp3e[1] = tmp3e[1] | w2/2*cos(w2/2*z)/z - sin(w2/2*z)/sqr(z);
      z = Re(zm0) + I*SupIm(zm0);
      tmp1e[0] = tmp1e[0] | sin(w0/2*z)/z;
      tmp2e[0] = tmp2e[0] | sin(w1/2*z)/z;
      tmp3e[0] = tmp3e[0] | sin(w2/2*z)/z;
      tmp1e[1] = tmp1e[1] | w0/2*cos(w0/2*z)/z - sin(w0/2*z)/sqr(z);
      tmp2e[1] = tmp2e[1] | w1/2*cos(w1/2*z)/z - sin(w1/2*z)/sqr(z);
      tmp3e[1] = tmp3e[1] | w2/2*cos(w2/2*z)/z - sin(w2/2*z)/sqr(z);

      tmp1e[1] = tmp1e[1]*get_j_coef(zm, 1);
      tmp2e[1] = tmp2e[1]*get_j_coef(zm, 1);
      tmp3e[1] = tmp3e[1]*get_j_coef(zm, 1);

      citaylor tmp1 = citaylor(tmp1e);
      citaylor tmp2 = citaylor(tmp2e);
      citaylor tmp3 = citaylor(tmp3e);
      part1 = part1*tmp1;
      part2 = part2*tmp2;
      part3 = part3*tmp3;
    }
  } else {
    part1 = part1*sin(w0/2*zm)/zm;
    part2 = part2*sin(w1/2*zm)/zm;
    part3 = part3*sin(w2/2*zm)/zm;
  }
  return -(part1 + part2 + part3);
  return -(alpha(z2) - alpha(z1))/(z2 - z1);
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

  // Example 3.5 - 1d real saddle point problem
  const real A0 = 1;
  const real A1 = 1;
  const real A2 = 1;
  const real Ip = 0.5;
  const int np = 3;
    f1 = sqr(p + (A0*cos(np*z1) + A1*cos((np + 1)*z1) + A2*cos((np - 1)*z1))) + 2*Ip;
  //f1 = sqr(p + A0*cos(z1)) + 2*Ip;
  f2 = z2;
  return;

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
  */

  // Example 6 - 2d HHG - simple saddle point problem
  // What should Omega be?
  real Omega = 1;
  citaylor ksz2z1 = ks_simple(z2, z1);
  f1 = sqr(ksz2z1 + A_simple(z1)) + 2*Ip;
  f2 = sqr(ksz2z1 + A_simple(z2)) + 2*(Ip + Omega);
  return;

  /*
  // Example 7 - 2d ATI - simple saddle point problem
  // What should Omega and p be?
  real Omega = 1;
  real p = 1;
  citaylor ksz2z1 = ks_simple(z2, z1, ok);
  f1 = sqr(ksz2z1 + A_simple(z1)) + 2*Ip;
  f2 = (p + A_simple(z2)) - sqr(ksz2z1 + A_simple(z2));
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
