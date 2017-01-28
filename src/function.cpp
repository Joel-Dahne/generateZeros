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

citaylor A(citaylor &z) {
  return A0*cos(w0*z + phi) + A1*cos(w1*z + phi) + A2*cos(w2*z + phi);
}

citaylor alpha(citaylor &z) {
  return alpha0*sin(w0*z + phi) + alpha1*sin(w1*z + phi) + alpha2*sin(w2*z + phi);
}

citaylor ks(citaylor &z1, citaylor &z2) {
  return -(alpha(z2) - alpha(z1))/(z2 - z1);
}
//End of stuff used in example 3

//Function 1
citaylor function1(citaylor &z1, citaylor &z2, bool &ok,
                   interval p) {
  
  //Example 1 - the identity function on C^2
  return z1;
  
  /*
  //Example 2 - The function f(z1, z2) = (sin(z1) + (z1)^2 +
  //e^{z2} - cos(2(z2)), cos(z1) + (z2)^3 + e^{2(z2)} - 2)
  return sin(z1) + sqr(z1) + exp(z2) - cos(2*z2);
  */
  /*
  //Example 2 - Hard coded derivative
  cinterval z1i = get_j_coef(z1, 0);
  cinterval z2i = get_j_coef(z2, 0);
  if (get_order(z1) == 0) {
    return citaylor(0, sin(z1i) + sqr(z1i) + exp(z2i) - cos(2*z2i));
  } else {
    if (get_j_coef(z1, 1) == cinterval(1)) {
      civector z(0, 1);
      z[0] = sin(z1i) + sqr(z1i) + exp(z2i) - cos(2*z2i);
      z[1] = cos(z1i) + 2*z1i;
      return citaylor(z);
    } else if (get_j_coef(z2, 1) == cinterval(1)) {
      civector z(0, 1);
      z[0] = sin(z1i) + sqr(z1i) + exp(z2i) - cos(2*z2i);
      z[1] = exp(z2i) + 2*sin(2*z2i);
      return citaylor(z);
    }
  }
  */
  /*
  //Example 3 - 2d real saddle point problem
  //First check that the function is well defined on the domain
  if (0 <= get_j_derive(z1, 0) - get_j_derive(z2, 0)) {
    ok = false;
    return z1;
  }
  //Evaluate the function
  return sqr(ks(z1, z2) + A(z1)) + 2*Ip;
  */
  /*
  //Example 4 - Simplified version of example 3
  return real(2.25)*sqr(sin(z1-z2) - sin(z2) + cos(z1-z2)*z2) + sqr(z2);
  */
  /*
  //Example 5 - A polynomial function
  return 4e-5*sqr(sqr(z1))*z1*sqr(z2)+2e-3*z1*sqr(sqr(z2))+2*sqr(z1)*z2-z2+real(0.75);
  */
}

//Function 2
citaylor function2(citaylor &z1, citaylor &z2, bool &ok,
                   interval p) {
  
  //Example 1 - the identity function on C^2
  return z2;
  
  /*
  //Example 2 - The function f(z1, z2) = (sin(z1) + (z1)^2 +
  //e^{z2} - cos(2(z2)), cos(z1) + (z2)^3 + e^{2(z2)} - 2)
  return cos(z1) + sqr(z2)*z2 + exp(2*z2) - 2;
  */
  /*
  //Example 2 - Hard coded derivative
  cinterval z1i = get_j_coef(z1, 0);
  cinterval z2i = get_j_coef(z2, 0);
  if (get_order(z2) == 0) {
    return citaylor(0, cos(z1i) + sqr(z2i)*z2i + exp(2*z2i) - 2);
  } else {
    if (get_j_coef(z1, 1) == cinterval(1)) {
      civector z(0, 1);
      z[0] = cos(z1i) + sqr(z2i)*z2i + exp(2*z2i) - 2;
      z[1] = -sin(z1i);
      return citaylor(z);
    } else if (get_j_coef(z2, 1) == cinterval(1)) {
      civector z(0, 1);
      z[0] = cos(z1i) + sqr(z2i)*z2i + exp(2*z2i) - 2;
      z[1] = real(3)*sqr(z2i) + real(2)*exp(2*z2i);
      return citaylor(z);
    }
  }
  */
  /*
  //Example 3 - 2d real saddle point problem
  //First check that the function is well defined on the domain
  if (0 <= get_j_derive(z1, 0) - get_j_derive(z2, 0)) {
    ok = false;
    return z2;
  }
  //Evaluate the function
  citaylor Az2 = A(z2);
  return sqr(p + Az2) - sqr(ks(z1, z2) + Az2);
  */
  /*
  //Example 4 - Simplified version of example 3
  citaylor part = sin(z1 - z2) - sin(z2);
  return real(2.25)*(sqr(cos(z1)*z2 + part) -
                     sqr(part + cos(z1-z2)*z2)) - 5*sqr(z2);
  */
  /*
  //Example 5 - A polynomial function
  return 3e-4*z1*sqr(sqr(z2))-7e-6*z1*sqr(z1)+2*z1*sqr(z2)-z1+real(0.75);
  */
}


//*********************************
//These functions are used to make calling the above functions in
//different ways easier, they normally do not need to be changed.
//*********************************

//Get an interval enclosure of the function on the domain
civector function(civector &domain, bool &ok,
                  interval p) {
  civector f(2);

  citaylor z1(0, 0), z2(0, 0);
  z1 = domain[1];
  z2 = domain[2];
  f[1] = get_j_derive(function1(z1, z2, ok, p), 0);
  f[2] = get_j_derive(function2(z1, z2, ok, p), 0);
  
  return f;
}

//Evaluate the function in a point
cvector midFunction(cvector &z, bool &ok,
                    interval p) {
  cvector f(2);

  citaylor z1(0, 0), z2(0, 0);
  z1 = z[1];
  z2 = z[2];
  f[1] = mid(get_j_derive(function1(z1, z2, ok, p), 0));
  f[2] = mid(get_j_derive(function2(z1, z2, ok, p), 0));

  return f;
}
