/*
 * decimal2binary.c
 *
 *  Created on: 10 Oct 2011
 *      Author: davidh
 */

  #include <stdio.h>

  int* decimal2binary(int decimal_value,int nBits)
  {
	  int i=0;
	  int j=nBits-1;
	  int* binary_value=malloc(nBits*sizeof(int));

	  for(i=0;i<nBits;i++)
	  {
		  binary_value[j--]=decimal_value%2;
		  decimal_value=decimal_value/2;
	  }

	  while((decimal_value%2)!=0);
	  return(binary_value);

  }

