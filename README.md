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

### A minimal example

We will here give a minimal example to check if the program compiled
correctly works.

When you download the code the preset function is f(z_1, z_2) = (z_1,
z_2), i.e the identity in the two dimensional complex space. This of
course have only one zero, at the origin. We will here use the program
to prove this.

After having ran `make all` go into the `build` directory and run the
program with

```
./integrate -- -1 1 -1 1 -1 1 -1 1
```

This tells the program compute the logarithmic integral for the domain
z_1 = [-1, 1] x i[-1, 1], z_2 = [-1, 1] x i[-1, 1]. The answer is the
number of zeros, counting multiplicity, in the domain. Since we know
that the number of zeros in the domain is exactly one it should output
an interval containing the only integer 1. The output is

```
([  0.603480,  1.524247],[ -0.409304,  0.409304])
```

indicating that the integral has a value with the real part inside the
first interval and the imaginary part inside the second interval. We
see that the value 1 is indeed inside this interval.

## Examples

TODO

## Authors

TODO

## License

TODO
