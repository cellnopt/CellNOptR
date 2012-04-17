#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double hill_function(double x,double n,double k)
{
	return(pow(x,n)/(pow(x,n)+pow(k,n)));
}
