# SYMmetric POWer elliptic curve L-functions

__SYMPOW__ is a mathematical program to compute special values of symmetric
power elliptic curve L-functions; it can compute up to about 64 digits
of precision.

## Installation

### Sage

Note that __SYMPOW__ is a dependency of [SageMath](https://www.sagemath.org/):
if [SageMath](https://www.sagemath.org/) is installed on your computer,
__SYMPOW__ must be already present.

### Distributions

Otherwise, when applicable, it is highly recommended to install __SYMPOW__ through
the package manager system of your Operating System.

```sh
# Debian derivatives
sudo apt install sympow
```

### By hand

As a (very) last resort, you can always build and install __SYMPOW__
as you are used to do with any other classical Open Source Software
that comes with its own building-installing machinery.

Build-essential GNU tools along [PARI/GP](https://pari.math.u-bordeaux.fr/)
must be present.

```sh
# overview
sh Configure
make all
make install
make distclean
```

## Usage examples

### basic usages

```sh
sympow -sp 2p16 -curve "[1,2,3,4,5]"
```
will compute L(Sym^2 E,edge) for E=[1,2,3,4,5] to 16 digits of precision
The result
```sh
 2n0: 8.370510845377639E-01
 2w0: 8.370510845377586E-01
```
consists of two calculations using different parameters to test the
functional equation (to see if these are sufficiently close).

```sh
# heavy computation
sympow -sp 3p12d2,4p8 -curve "[1,2,3,4,5]"
```
will compute the 0th-2nd derivatives of L(Sym^3 E,center) to 12 digits
and L(Sym^4 E,edge) to 8 digits.

```sh
# light alternatives
sympow -sp 3p12d2,4p8 -curve "[0,1,2,3,4]"
sympow -sp 3p12d2,4p8 -curve "[1,1,0,1,0]"
```

Special cases: When a curve has CM, Hecke power can be used instead
```sh
sympow -sp 7p12d1 -hecke -curve "[0,0,1,-16758,835805]"
```
will compute the 0th-1st derivatives of L(Sym^7 psi,special) to 12 digits

Bloch-Kato numbers can be obtained for powers not congruent to 0 mod 4:
```sh
sympow -sp 2bp16 -curve "[1,2,3,4,5]"
```
should return
```sh
 2n0: 4.640000000000006E+02
 2w0: 4.639999999999976E+02
```
which can be seen to be very close to the integer 464.

Modular degrees can be computed with the -moddeg option.
```sh
sympow -curve "[1,2,3,4,5]" -moddeg
```
should return
```sh
Modular Degree is 464
```

Analytic ranks can be computed with the -analrank option.
```sh
sympow -curve "[1,2,3,4,5]" -analrank
```
should return
```sh
Analytic Rank is 1 : L'-value 3.51873e+00
```

and (if the mesh file for the fifth derivative is present)
```sh
sympow -curve "[0,0,1,-79,342]" -analrank
```
should return
```sh
Analytic Rank is 5 : leading L-term 3.02857e+01
```

### add new symmetric powers

__SYMPOW__ keeps data for symmetric powers in the datafiles directories
If a pre-computed mesh of inverse Mellin transform values is not
available for a given symmetric power, __SYMPOW__ will fail.
If __GP__ is available, the command
```sh
sympow -new_data 2
```
will add the data for the 2nd symmetric power, while
```sh
sympow -new_data 3d2
```
will add the data for the 2nd derivative and 3rd symmetric power,
```sh
sympow -new_data 6d0h
```
will add the data for the 0th derivative of the 6th Hecke translate, and
```sh
sympow -new_data 4c
```
will add data for the 4th symmetric power for curves with CM
(these need to be done separately for powers divisible by 4).

The mesh files are stored in binary form, and thus endian-ness is a
worry when moving from one platform to another.

To enable modular degree computations, the 2nd symmetric powers must be
extant, and analytic rank requires derivatives of the 1st power. There are
lower-precision datafiles in the sympow.tar distribution to expedite the
computation of analytic rank and modular degree. These latter can also be
re-generated in about 5-10 minutes via running armd.sh as a shell script.

For higher precision, here is a test case; it gives 6.0e+2 to almost 64 digits:
```sh
sympow -sp 2bp64 -curve '[0,0,0,0,5]'
```

For Hecke powers, the special value for the (2N-1)st power is computed
using the datafile for the Nth translate, and the (2N)th power requires
both the Nth and (N+1)st.

## Notes

* Calculation of periods should not be trusted to the full 64 digits.
* Analytic rank does not save the ap when more than 10^9 terms are needed:
  in such a case, it is probably easier to run sympow as
```sh
    sympow -sp 1w0p6D0,1w0p6D2,1w0p6D4,1w0p6D6,1w0p6D8
```
  [for even signed curves, and however many derivatives are deemed
   relevant, and to a desired precision] and parse the output.

## Bug reports, feature requests and ideas

Please submit _bug reports_, _feature requests_ or _ideas_ through the GitLab website
[https://gitlab.com/rezozer/forks/sympow](https://gitlab.com/rezozer/forks/sympow) .

## History

The **quad-double** package was modified from David Bailey's package:
 crd.lbl.gov/~dhbailey/mpdist/

The **squfof** implementation was modified from Allan Steel's version
of Arjen Lenstra's original **LIP**-based code.

The **ec_ap** code was originally written for the kernel of MAGMA,
but was modified to use small integers when possible.

__SYMPOW__ was originally developed using [PARI/GP](https://pari.math.u-bordeaux.fr/),
but due to licensing difficulties, this was eliminated. __SYMPOW__ also does not
use the standard math libraries, thus eliminating possible license conflicts.
However, the user should ensure that there are no C library license issues.

__SYMPOW__ still can use __GP__ to compute the meshes of inverse Mellin
transforms (this is done when a new symmetric power is added to datafiles),
but the datafiles can alternatively be downloaded from the SYMPOW website.
The current license for __GP__ appears to allow this usage. Please ensure
that the __GP__ you have allows this.

## Versioning

The __forked source versioning__ is based on the original source versioning as follows:
 - the _major version number_ is incremented by one
 - the _minor version number_ remains unchanged
 - the _micro version number_ (or _patch version number_) numbers the fork version

Notice that the two first forked versions are special:
 - sympow forked version _2.023.0_ is the original source itself (modulo this README.md file)
 - sympow forked version _2.023.1_ is the source as patched by Debian on June 2018 (Debian version _1.023-9_)

## Legal mumbo-jumbo

### original code

Copyright: 2005-2018 Mark Watkins <watkins@maths.usyd.edu.au>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
  * Redistribution of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
  * Redistribution in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
  * If redistribution is done as a part of a compilation that has a more
 restrictive license (such as the GPL), then the fact that SYMPOW has
 a less restrictive license must be made clear to the recipient.
 For example, a line like (include bracketed text if SYMPOW is modified):
  "This compilation includes [a modification of] SYMPOW whose [original]
   code has a less-restrictive license than the entire compilation."
 should appear in a suitable place in the COPYING and/or LICENSE file.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

### forked code

Copyright: 2018 Jerome G. M. Benoit <jgmbenoit@rezozer.net>

