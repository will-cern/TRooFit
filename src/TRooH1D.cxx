
#include "TRooFit/TRooH1D.h"


//___________________________________
/* BEGIN_HTML
<p>A TRooH1D is the TRooFit version of a TH1D</p>
<p>You must first create a RooFit variable (e.g. RooRealVar or RooCategory) before you can
 construct a TRooH1D</p>
 
<p>Examples:
<pre>
RooRealVar x("x","x",0,10); //create a continuous variable with range 0->10
TRooH1D h("h","h",x,5); //create a 5-bin TRooH1D, that will be a pdf for "x"
<br>
RooCategory c("c","c");
c.defineType("A");c.defineType("B"); //create a discrete variable, values are A or B
TRooH1D h2("h2","h2",c); //creates a 2-bin TRooH1D, that will be a pdf for "c"
</pre>
</p>


<p>The histogram can then be filled like any normal histogram:
<pre>
h.Fill(4,3);
h.SetBinContent(4,3);
</pre>
</p>

<p>If the observable was a discrete variable (RooCategory) you can fill using the labels:
<pre>
h.Fill("A",3);
h.SetBinContent("A",3);
</pre>
</p>

<p>When you use the TRooH1::Fill method, the histogram will automatically create a 
nuisance parameter to represent the statistical uncertainty for that bin. 
The nuisance parameter appears as a <i>shapeFactor</i> for that bin.</p>

<p>There are three types of factors:
<ul>
 <li><b><i>normFactors</i></b> : These multiply the content of every bin. Add one with TRooAbsH1::addNormFactor </li>
 <li><b><i>shapeFactors</i></b> : These multiply the content of a specific bin. Add one with TRooAbsH1::addShapeFactor</li>
 <li><b><i>transFactors</i></b> : These connect one TRooFit histogram to another, 
                  and are advanced use case (documentation to come).
                  There can only be one transFactor for a histogram.
                  If a histogram has a transFactor, the value of the histogram (before normFactor and shapeFactors)
                  is given by transFactor*denominatorHistogram</li>
</ul>

</p>
END_HTML */
//____________________________________


ClassImp(TRooH1D) 