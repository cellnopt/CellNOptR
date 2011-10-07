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


