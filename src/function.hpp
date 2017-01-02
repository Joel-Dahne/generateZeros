#ifndef _FUNCTION_H
#define _FUNCTION_H

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

citaylor function1(citaylor &z1, citaylor &z2, bool &ok,
                   interval p);

citaylor function2(citaylor &z1, citaylor &z2, bool &ok,
                   interval p);

civector function(civector &domain, bool &ok,
                  interval p);

cvector midFunction(cvector &z, bool &ok,
                    interval p);

#endif
