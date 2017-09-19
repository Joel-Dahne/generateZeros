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
void function(citaylor &f1, citaylor &f2, citaylor &z1, citaylor &z2,
                  bool &ok, interval p);

citaylor function1(citaylor &z1, citaylor &z2, bool &ok,
                   interval p);

citaylor function2(citaylor &z1, citaylor &z2, bool &ok,
                   interval p);

civector intervalFunction(civector &domain, bool &ok,
                  interval p);

cvector midFunction(cvector &z, bool &ok,
                    interval p);

#endif
