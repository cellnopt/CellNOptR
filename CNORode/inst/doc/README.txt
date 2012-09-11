Notes for developers. 

In order to generate the doc, type:

    R CMD Sweave CNORfuzzy-vignette.Rnw
	pdflatex CNORfuzzy-vignette.tex
	pdflatex CNORfuzzy-vignette.tex

You will need the figures that are provided in this package. The figures related to the the full analysis are saved in the SVN. The figures of the toy model are not stored in the SVN since they can be generated in a couple of minutes. 

To speed up the doc generation (for debugging), edit the Rnw file and set eval=false in the relevant code snipset (search for TODO pattern).

The generation of all the plots takes about 5-6 hours with the current version
(March 2012). Can be reduced significantly by playing with params$Maxtime and
params$optimisation$maxtime and also N=1 instead of N=10 but results will not be
meaningful in the plotMeanFuzzy case.

TC., march 2012



