# generateZeros

A program for finding and validating all zeros to a system of two
analytic functions in a rectangular domain.

## Getting Started

### C-XSC

To get the program running you will need an installed copy of
C-XSC. This is a C++ library for handling interval arithmetic's.
Information about the library and instructions for installing can be
found on their [website](http://www2.math.uni-wuppertal.de/~xsc/). The
version used for this program is 2.5.4.

Unfortunately the complex interval arithmetic in C-XSC is not thread
safe and to be able to use it in a multi-threaded environment a small
change to the code is needed. We will here describe the steps to do
so.

In the directory were you installed C-XSC, usually ~/cxsc, you will
need to change in the file `path-to-cxsc/include/cinterval.inl`. On
line 644 change from

```
return mult_operator(a,b);
```

to

```
return cinterval(Re(a)*Re(b)-Im(a)*Im(b), Re(a)*Im(b)+Im(a)*Re(b));
```

and on line 657 change from

```
return div_operator(a,b);
```

to

```
return (a * cinterval(Re(b), -Im(b))) / (sqr(Re(b)) + sqr(Im(b)));
```

This all you need to change and you will not need to recompile the
library.

### Compiling the program

To be able to compile the program you will first need to add the path
to your C-XSC installation in the `makefile`. Set the variable
`CXSCDIR` to point to the path of your installation, usually this is

```
CXSCDIR=/home/USERNAME/cxsc
```

To compile the program call `make all` while in the root of the
project. This will create the two programs `integrate` and
`generateZeros` in the directory `build`. You can also compile
optimized versions of the program, this gives longer compiling times
but improves the performance of the program. To do this call `make
build/integrateOpt` or `make build/generateZerosOpt` which will create
the two programs `integrateOpt` and `generateZerosOpt` in the `build`
directory.

## Examples

TODO

## Authors

TODO

## License

TODO
