//============================================================================
// Name        : simulatorT1.c
// Author      : Thomas Cokelaer
// Version     : 0.1
// Copyright   :
// Description :
//============================================================================

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>

// keep this NA > 1 and integer
#define NA 100

SEXP simulatorT1 (

    SEXP nStimuli_in,
    SEXP nInhibitors_in,
    SEXP nCond_in,
    SEXP nReacs_in,
    SEXP nSpecies_in,
    SEXP nSignals_in,
    SEXP nMaxInputs_in,

    SEXP finalCube_in,
    SEXP ixNeg_in,
    SEXP ignoreCube_in,
    SEXP maxIx_in,

    SEXP valueGCube_in,
    SEXP valueKCube_in,
    SEXP valueNCube_in,

    SEXP indexSignals_in,
    SEXP indexStimuli_in,
    SEXP indexInhibitors_in,

    SEXP valueInhibitors_in,
    SEXP valueStimuli_in
) {

    SEXP simResults;

    int counter = 0;
    int i = 0;
    int j = 0;
    float curr_max = 0;
    float or_max = 0;
    int selection[40];
    int selCounter = 0;
    double *rans;

    int nStimuli = INTEGER(nStimuli_in)[0];
    int nInhibitors = INTEGER(nInhibitors_in)[0];
    int nCond = INTEGER(nCond_in)[0];
    int nReacs = INTEGER(nReacs_in)[0];
    int nSpecies = INTEGER(nSpecies_in)[0];
    int nSignals = INTEGER(nSignals_in)[0];
    int nMaxInputs = INTEGER(nMaxInputs_in)[0];
    int nCube2 = INTEGER(nMaxInputs_in)[0];
    int nCube1;// defined later.


    int *maxIx;
    maxIx = (int*) malloc(nReacs * sizeof(int));
    for (i = 0; i < nReacs; i++) {
        maxIx[i] = INTEGER(maxIx_in)[counter++];
    }

    counter = 0;
    int *indexStimuli;
    indexStimuli = (int*) malloc(nStimuli * sizeof(int));
    for (i = 0; i < nStimuli; i++) {
        indexStimuli[i] = INTEGER(indexStimuli_in)[counter++];
    }

    counter = 0;
    int *indexInhibitors;
    indexInhibitors = (int*) malloc(nInhibitors * sizeof(int));
    for (i = 0; i < nInhibitors; i++) {
        indexInhibitors[i] = INTEGER(indexInhibitors_in)[counter++];
    }

    counter = 0;
    int *indexSignals;
    indexSignals = (int*) malloc(nSignals * sizeof(int));
    for (i = 0; i < nSignals; i++) {
        indexSignals[i] = INTEGER(indexSignals_in)[counter++] ; 
    }

    counter=0;
    int **finalCube;
    finalCube = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        finalCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            finalCube[i][j] = INTEGER(finalCube_in)[j*nReacs+i];
        }
    }

    counter=0;
    int **ixNeg;
    ixNeg = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        ixNeg[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            ixNeg[i][j] = INTEGER(ixNeg_in)[j*nReacs+i];
        }
    }

    counter=0;
    int **ignoreCube;
    ignoreCube = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        ignoreCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            ignoreCube[i][j] = INTEGER(ignoreCube_in)[j*nReacs+i];
        }
    }

    counter=0;
    float **valueInhibitors;
    valueInhibitors = (float**) malloc(nCond * sizeof(float*));
    for (i = 0; i < nCond; i++) {
        valueInhibitors[i] = (float*) malloc(nInhibitors * sizeof(float));
        for (j = 0; j < nInhibitors; j++) {
            valueInhibitors[i][j] = REAL(valueInhibitors_in)[nCond*j+i];
        }
    }

    counter=0;
    float **valueStimuli;
    valueStimuli = (float**) malloc(nCond * sizeof(int*));
    for (i = 0; i < nCond; i++) {
        valueStimuli[i] = (float*) malloc(nStimuli * sizeof(float));
        for (j = 0; j < nStimuli; j++) {
            valueStimuli[i][j] = REAL(valueStimuli_in)[nCond*j+i];
        }
    }



    counter = 0;
    float **valueGCube;
    nCube1 = nCond*nReacs;
    valueGCube = (float**) malloc(nCube1 * sizeof(float*));
    for(i = 0; i<nReacs; i++){
        for(j = 0; j<nCond; j++){
            valueGCube[counter] = (float*) malloc(nCube2 * sizeof(float));
            for (int k=0; k<nCube2; k++){
                valueGCube[counter][k] = REAL(valueGCube_in)[i+k*nReacs];
            }
            counter++;
        }
    }

    counter = 0;
    float **valueKCube;
    nCube1 = nCond*nReacs;
    valueKCube = (float**) malloc(nCube1 * sizeof(float*));
    for(i = 0; i<nReacs; i++){
        for(j = 0; j<nCond; j++){
            valueKCube[counter] = (float*) malloc(nCube2 * sizeof(float));
            for (int k=0; k<nCube2; k++){
                valueKCube[counter][k] = REAL(valueKCube_in)[i+k*nReacs];
            }
            counter++;
        }
    }


    counter = 0;
    float **valueNCube;
    nCube1 = nCond*nReacs;
    valueNCube = (float**) malloc(nCube1 * sizeof(float*));
    for(i = 0; i<nReacs; i++){
        for(j = 0; j<nCond; j++){
            valueNCube[counter] = (float*) malloc(nCube2 * sizeof(float));
            for (int k=0; k<nCube2; k++){
                valueNCube[counter][k] = REAL(valueNCube_in)[i+k*nReacs];
            }
            counter++;
        }
    }

    // build this pow matrix once for all. slightly faster that way.
    float **powkn;
    powkn = (float**) malloc(nCube1 * sizeof(float*));
    for(i = 0; i<nReacs*nCond; i++){
        powkn[i] = (float*) malloc(nCube2 * sizeof(float));
        for (j=0; j<nCube2; j++){
            powkn[i][j] = pow(valueKCube[i][j], valueNCube[i][j]);
        }
    }

    //============================================================================

    // fill end_ix - how many reactions have each species as output
    int end_ix[nSpecies];
    int count_species=0;
    for(i = 0; i < nSpecies; i++) {
        for(j = 0; j < nReacs; j++) {
            if(i == maxIx[j]) {
                count_species++;
            }
        }
        end_ix[i] = count_species;
        count_species = 0;
    }

    // see stop conditions
    float test_val = 1e-3;

    // create an initial values matrix
    float init_values[nCond][nSpecies];
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            init_values[i][j] = NA;
        }
    }

    // set the initial values of the stimuli
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nStimuli; j++) {
            init_values[i][indexStimuli[j]] = valueStimuli[i][j];
        }
    }

    // flip and redefine inhibitors
    if(nInhibitors) {
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nInhibitors; j++) {
                valueInhibitors[i][j] = 1 - valueInhibitors[i][j];
                if(valueInhibitors[i][j] == 1) {
                    valueInhibitors[i][j] = NA;
                }
            }
        }
    }

    // set the initial values of the inhibitors
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nInhibitors; j++) {
            init_values[i][indexInhibitors[j]] = valueInhibitors[i][j];
        }
    }

    // initialize main loop
    float output_prev[nCond][nSpecies];
    float new_input[nCond][nSpecies];
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            new_input[i][j] = (float)init_values[i][j];
        }
    }

    int term_check_1 = 1;
    float term_check_2 = 1;
    int count = 1;
    float diff;

    // define the temp data
    float temp_store[nCond * nReacs][nMaxInputs];
    float transval[nCond * nReacs][nMaxInputs];


    //============================================================================
    float x;
    float transval_temp;



    while(term_check_1 && term_check_2) {

        // copy to outputPrev
        memcpy(output_prev, new_input, sizeof(output_prev));

        // fill temp store
        // this is different to R version, through a single loop
        // with conditions
        int track_cond = 0; // track condition
        int track_reac = 0; // track reaction
        for(i = 0; i < nCond * nReacs; i++) {
            for(j = 0; j < nMaxInputs; j++) {
                // initial values of each input
                temp_store[i][j] = output_prev[track_cond][finalCube[track_reac][j]];

                if(ignoreCube[track_reac][j]) {
                    temp_store[i][j] = NA;
                }

                x = temp_store[i][j];
                if (x == NA){
                    transval_temp = NA;
                }
                else {
                    // slightly faster
                    float powxn = pow(x, valueNCube[i][j]);
                    if ( (powxn + powkn[i][j]) == 0.){
                        transval_temp = NA;
                    }
                    else{
                        transval_temp = valueGCube[i][j] * (1. + powkn[i][j]) * powxn / ( powxn + powkn[i][j] );
                    }
                }


                if(ixNeg[track_reac][j]) {
                    // flip the values of the neg inputs
                    if (transval_temp!=NA){
                        transval_temp = 1 - transval_temp;
                    }
                }
                transval[i][j] =transval_temp;
            }

            track_cond++;
            if((track_cond == nCond)) {
                track_cond = 0;
                track_reac++;
            }
        }


        // compute the AND gates (find the min 0/1 of each row)
        float output_cube[nCond][nReacs]; // declare output_cube

        float current_min;
        int dial_reac = 0;
        int dial_cond = 0;
        for(i = 0; i < nCond * nReacs; i++) {
            current_min = transval[i][0];
            for(j = 1; j < nMaxInputs; j++) {

            //TODO TO CLEAN. SEE boolean case.
            // this code is currently same as Melody's but we may want to incorporate/change similarly 
            //to C code in the boolean case
                // if statement below is for AND gates with any NA  input
                // in this case output should always be NA
            //    if(transval[i][j] == NA && ignoreCube[dial_reac][j] == 0) {
            //        current_min = NA;
            //        break;
                //}
                //else
                //
                //rifgr now, original code is buggy and ignoreNA, so let us do
                //the same.
                if (transval[i][j]!=NA){
                     if(transval[i][j] < current_min) {current_min = transval[i][j];}
                }
            }

            output_cube[dial_cond][dial_reac] = current_min;
            dial_cond++;
            if(dial_cond==nCond) {dial_cond = 0; dial_reac++;}
        }

        // compute the OR gates and reinitialize new_input
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                new_input[i][j] = NA;
            }
        }

        // declare vector to store 'selection' (R)
        selCounter = 0;
        for(int s = 0; s < nSpecies; s++) {
            // is the species an output for any reactions?
            if(end_ix[s]) {

                // find reactions with this species as output
                // add equivalent output_cube data to new_input
                for(int a = 0; a < nReacs; a++) {
                    if(s == maxIx[a]) {selection[selCounter] = a; selCounter++;}
                }
                // if the species is an output for a single reaction
                // it's a 1-1 mapping to new_input
                if(selCounter == 1) {
                    for(int b = 0; b < nCond; b++) {
                        new_input[b][s] = output_cube[b][selection[selCounter-1]];
                    }
                    selCounter = 0;
                }
                // else if species is output for > 1
                if(selCounter > 1) {
                    for(i=0; i < nCond; i++) {
                        or_max = NA;
                        curr_max = 0;
                        for(int p=0; p < selCounter; p++) {
                            if(output_cube[i][selection[p]] >= curr_max && output_cube[i][selection[p]] < NA) {
                                or_max = output_cube[i][selection[p]];
                                curr_max = output_cube[i][selection[p]];
                            }
                        }
                        new_input[i][s] = or_max;
                    }
                    selCounter = 0;
                }
            }
        }

        // reset the stimuli
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nStimuli; j++) {
                curr_max = valueStimuli[i][j];
                if(new_input[i][indexStimuli[j]] > curr_max && new_input[i][indexStimuli[j]] < NA) {
                    curr_max = new_input[i][indexStimuli[j]];
                }
                new_input[i][indexStimuli[j]] = curr_max;
            }
        }

        // reset the inhibitors
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nInhibitors; j++) {
                if(valueInhibitors[i][j] == 0) {
                    new_input[i][indexInhibitors[j]] = 0;
                }
            }
        }

        // set 'NAs' (2s) to 0
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                if(new_input[i][j] == NA) {new_input[i][j] = 0;}
                if(output_prev[i][j] == NA) {output_prev[i][j] = 0;}
            }
        }

        term_check_1 = 0;
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                diff = fabs((new_input[i][j] - output_prev[i][j]));
                if (diff > test_val){
                    term_check_1 = 1;
                    break;  /*  no need to keep going checking other values if
                                one is greater than test_val */
                }
            }
        }
        /*term_check_1 = !(abs(diff) < test_val);*/
        term_check_2 = (count < (nSpecies * 1.2));
        count++;

    } // end of main loop
    // set non-resolved bits to 2 (NA)
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            if(new_input[i][j] != output_prev[i][j])
                new_input[i][j] = NA;
        }
    }

     PROTECT(simResults = allocMatrix(REALSXP, nCond, nSpecies));
    rans = REAL(simResults);
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            if(new_input[i][j] == NA) rans[i + nCond*j] = NA_REAL;
            else rans[i + nCond*j] = new_input[i][j];
        }
    }


/*     PROTECT(simResults = allocMatrix(REALSXP, nCond, nSignals));
    rans = REAL(simResults);
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSignals; j++) {
            if(new_input[i][indexSignals[j]] == NA) rans[i + nCond*j] = NA_REAL;
            else rans[i + nCond*j] = new_input[i][indexSignals[j]];
        }
    }

*/

    free(maxIx);
    free(indexStimuli);
    free(indexInhibitors);
    free(indexSignals);

    for (i = 0; i < nReacs; i++) {
        free(finalCube[i]);
    }
    free(finalCube);

    for (i = 0; i < nReacs; i++) {
        free(ixNeg[i]);
    }
    free(ixNeg);

    for (i = 0; i < nReacs; i++) {
        free(ignoreCube[i]);
    }
    free(ignoreCube);

    for (i = 0; i < nCond; i++) {
        free(valueInhibitors[i]);
    }
    free(valueInhibitors);

    for (i = 0; i < nCond; i++) {
        free(valueStimuli[i]);
    }
    free(valueStimuli);

    for (i = 0; i < nCube1; i++) {
        free(valueGCube[i]);
    }
    free(valueGCube);

    for (i = 0; i < nCube1; i++) {
        free(valueKCube[i]);
    }
    free(valueKCube);

    for (i = 0; i < nCube1; i++) {
        free(powkn[i]);
    }
    free(powkn);

    for (i = 0; i < nCube1; i++) {
        free(valueNCube[i]);
    }
    free(valueNCube);

    UNPROTECT(1);
    return simResults;

}
