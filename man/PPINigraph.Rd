\name{PPINigraph}
\docType{data}
\alias{PPINigraph}
\title{Protein-protein interaction netwrok}
\description{The human protein-protein interaction network was built using a unified PPI dataset obtained as APID (Prieto,C. and De Las Rivas,J. 2006), by the combination of interactions coming from six source databases. The starting whole dataset was composed by 68488 human physical protein-protein interactions validated at least by one experimental method and reported in one article published in PubMed. From this dataset we obtained two PPI subsets with increasing confidence: a set of 28971 interactions validated by at least one binary experimental method (binary as defined in  (De Las Rivas,J. and Fontanillo,C. 2010)); a set 6033 interactions validated by at least two experimental methods, one of them binary.}
\usage{PPINigraph}
\format{PPINigraph is an igraph with proteins as nodes and undirected links as physical protein interactions.}
\source{This network was bult for the analysis performed in (Eduati,F. et al. 2012)}
\references{
\enumerate{
\item{F. Eduati, J. De Las Rivas, B. Di Camillo, G. Toffolo, J. Saez-Rodriguez. Integrating literature-constrained and data-driven inference of signalling networks. Bioinformatics, 28(18):2311-2317, 2012.}
\item{C. Prieto, J. De Las Rivas. APID: Agile Protein Interaction DataAnalyzer. Nucleic Acids Res., 34, W298-302, 2006.}
\item{J. De Las Rivas, C. Fontanillo. Protein-protein interactions essentials: key concepts to building and analyzing interactome networks. PLoS Comput.Biol., 6, e1000807, 2010.}
}
}
\keyword{datasets}
