\name{gseq} 
\alias{gseq} 
\title{Geometric Sequence Creator} 
\description{ 
   This simple function creates a geometric sequence
}
\usage{ 
gseq(from, to, by, addrem = FALSE, addzero = FALSE)
} 
\arguments{ 
   \item{from}{the initial value}
   \item{to}{the final value}
   \item{by}{the multiplication factor}
   \item{addrem}{if the final value is surpassed by multiplication, add.rem = TRUE will add it as the last element of the vector}
   \item{addzero}{add a zero at the beginning}
}
\examples{
### For thse sequence 0 .4, 4, 40, ..., 4e10
gseq(.4,4e10,10, addzero = TRUE)

### A truncated sequence without add.rem
gseq(1,100,2)

### A truncated sequence with add.rem
gseq(1,100,2,addrem = TRUE)

}
\keyword{file} 
