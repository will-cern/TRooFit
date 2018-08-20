Introduction
============

TRooFit is a collection of pdf classes for RooFit, that are designed to behave like objects from ROOT. For example, TRooFit introduces:

    - TRooH1D : A pdf (inherits from RooAbsPdf) that you can fill and draw like a TH1D. The value of this object is a DENSITY.
    - TRooHF1D : A function (inherits from RooAbsReal). The value of this object is a VALUE.
    - TRooHStack : A sum of TRooH1D (or any other RooAbsReal) ... think of it as a stack of histograms
    - TRooHPdfStack : Like a TRooHStack but only accepts (extended) pdfs, and ensures each is normalized individually in the summation, i.e. pdfValue = Sum(i) coef_i*pdfValue_i, where coef_i is the fraction of total expected events of that pdf. The coef_i are automatically recomputed 

Setup (CMake)
=====

Setup any release that has ROOT, e.g. AthAnalysis, and then use acm to clone this package and compile it (this will compile the master version of TRooFit)

```
mkdir build source
cd build
acmSetup AthAnalysis,21.2,latest
acm clone_project TRooFit will/TRooFit
acm compile
````

To checkout the master, just omit the last argument in the clone_project command.



Documentation
========
Example notebooks are available to help you get started ....

1a. [Model building with TRooWorkspaces](https://cernbox.cern.ch/index.php/s/dxJoaZC0W131OgZ)

2a. [Model inspection general demo](https://cernbox.cern.ch/index.php/s/tzyUPsHqBSgKK1A)
2b. [Obtaining yields and uncertainties](https://cernbox.cern.ch/index.php/s/I7dwstfWVMD2CIk)
2c. [Parameter estimation and systematic breakdown](https://cernbox.cern.ch/index.php/s/Mi2MtN8TvAdHocf)


Binder
=======
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/will-cern/TRooFit-binder/master)

