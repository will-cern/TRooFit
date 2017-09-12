Introduction
============

TRooFit is a collection of pdf classes for RooFit, that are designed to behave like objects from ROOT. For example, TRooFit introduces:

    - TRooH1D : A pdf that you can fill and draw like a TH1D
    - TRooHStack : A sum of TRooH1D (or any other RooAbsReal) ... think of it as a stack of histograms
    - TRooHPdfStack : Like a TRooHStack but only accepts (extended) pdfs, and ensures each is normalized individually in the summation, i.e. pdfValue = Sum(i) coef_i*pdfValue_i, where coef_i is the fraction of total expected events of that pdf. The coef_i are automatically recomputed 


Setup
=====

TRooFit currently is designed to work in AthAnalysisBase 2.4 series. Set that release up, clone the project, then compile:

```
asetup AthAnalysisBase,2.4.31,here
git clone https://:@gitlab.cern.ch:8443/will/TRooFit.git
cd TRooFit/cmt
cmt config
make
```

Examples
========
So far, try to look at the macros provided in the demo dir

The documentation for the latest version of TRooFit is also available at: http://will.web.cern.ch/will/TRooFit/