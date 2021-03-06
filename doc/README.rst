Installation
==================

Quick installation
----------------------

If you already have all R dependencies and work from the source, you could try::

    git clone https://github.com/cellnopt/CellNOptR.git
    cd CellNOptR
    # To install all packages, use the pyhon script
    python setup.py install
    # or to install CellNOptR only:
    cd packages
    cd CellNOptR
    R CMD INSTALL .

Standard installation
-------------------------

Before starting, make sure you have installed the latest version of R. For more information and download
of R, please refer to `R project page <http://www.r-project.org/>`_ . For more information about how to 
install R packages, please refer to `Installing package <http://cran.r-project.org/doc/manuals/R-admin.html#Installing-packages>`_
These packages rely on several Bioconductor package (e.g., RBGL, graph, methods, etc.). As an example, you can
install RGBL package by typing:
::

  source("http://bioconductor.org/biocLite.R")
  biocLite("RBGL")
  
To install CellNOptR, type::

  biocLite("CellNOptR")
  
or from a tar ball as follows:

    install.packages("path_to_CellNOptR/CellNOptR_1.0.0.tar.gz", + repos=NULL, type="source")

or, using the R GUI by clicking on "Packages & Data" then "Package installer", then choosing "local source"
from the dropdown menu, clicking "install", choosing CellNOptR.1.0.0.tar.gz
and finally clicking "open".

Then, you can also install other CellNOptR related packages::

   biocLite("CNORdt")
   biocLite("CNORfeeder")
   biocLite("CNORfuzzy")
   biocLite("CNORode")
