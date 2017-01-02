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
  cinterval A2 = cinterval(interval(1.5*1.5), interval(0));
  cinterval twoB = cinterval(interval(2*0.5), interval(0));
  cinterval twoC = cinterval(interval(2*2.5), interval(0));
  citaylor part = sin(z1 - z2) - sin(z2);
  
  return A2*(sqr(part) + 2*part*cos(z1-z2)*z2 +
             sqr(cos(z1-z2))*sqr(z2)) + twoB*sqr(z2);
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
  cinterval A2 = cinterval(interval(1.5*1.5), interval(0));
  cinterval twoB = cinterval(interval(2*0.5), interval(0));
  cinterval twoC = cinterval(interval(2*2.5), interval(0));
  citaylor part = sin(z1 - z2) - sin(z2);
  
  return z2*(A2*(2*part*cos(z1) + cos(z1) - 2*part*cos(z1-z2) -
                 sqr(cos(z1-z2))*z2) - twoC*z2);
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
