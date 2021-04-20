## Patch to v0.3.1
This is a patch to version 0.3.1. In this patch I have:

* fixed an 'ambiguous name' generating compilation errors on SOLARIS
* tidied C++ code for better readability

## Test environments
* Ubuntu 18.04 local install, R 3.4.4 and R-devel (20.12.2018)
* win-builder, R 3.4.4, R 3.5.1 and R-devel (19.12.2018)

## R CMD check results
There were no ERRORs nor WARNINGs.

There was 1 NOTE (all environments):
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Thomas Lugrin <thomas.lugrin@alumni.epfl.ch>’

New maintainer:
  Thomas Lugrin <thomas.lugrin@alumni.epfl.ch>
Old maintainer(s):
  Thomas Lugrin <thomas.lugrin@epfl.ch>

_I replaced a temporary e-mail address with a life-long alumni e-mail address._

There was 1 additional NOTE (win-builder, R 3.4.4 only):

* checking CRAN incoming feasibility ... NOTE

Possibly mis-spelled words in DESCRIPTION:
  Extremal (3:30)
  Heffernan (9:211)
  Tawn (9:225)
  extremal (9:38, 9:311, 9:463)
  pre (9:93)

_The word 'extremal' is commonly used in extreme value theory, here with the meaning 'sub-asymptotic'._

_'Heffernan' and 'Tawn' are proper nouns._

_'pre-' is a common prefix in English meaning prior, before, earlier._

## Downstream dependencies
There are currently no downstream dependencies for this package.

---
