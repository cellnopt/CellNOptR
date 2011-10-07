How to call a C function from R. Case 1 : no arguments
==========================================================

Compile the C function with the following command:

    R CMD SHLIB  func1.c -o func1.so

Then, open R and load the function:

    dyn.load("func1.so")

and call it:

    .C("func1")


Case 2: arguments
====================


    void func1(double *array, int *n)
    {
      int i = 0;
      for (i=0; i<n; i++)
         printf("%f\n", array[i]);

      printf("%d \n", *n);
    }

compile it with the same method as above.
Then in R, type:

    .C("func1", as.double(c(1,2)), as.integer(2))

