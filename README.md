# generateZeros

A program for finding and validating all zeros to a system of two
analytic functions in a rectangular domain.

## Getting Started

### C-XSC

To get the program running you will need an installed copy of
C-XSC. This is a C++ library for handling interval arithmetic's.
Information about the library and instructions for installing can be
found on their [http://www2.math.uni-wuppertal.de/~xsc/](website). The
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

## Examples

TODO

## Authors

TODO

## License

TODO
