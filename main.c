#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include "libs.h"
#include "strucs.h"

int main(int argc, char* argv[]){
    FILE* infoFile = NULL;
    FILE* contentFile = NULL;
    programInfo programInfo;
    char* contentFileName;
    previouslyKnownValues prevKnoVals;

    /* Matrix which save times involved in the different stages to detect bottlenecks */
    clock_t durationSteps[11][2] = {0};

    /* Row 0 of durationSteps devoted to the total program time */
    durationSteps[0][0] = clock();
    /* Row 1 of durationSteps devoted to loading libraries */
    durationSteps[1][0] = clock();

    if (argc != 2) {
        //infoFile = openFile(argv[1]);
        infoFile = openFile("CY62167_090nm_TS.jl");
        //infoFile = openFile("CY62167_130nm_TS.jl");

        if (infoFile != NULL) {
            contentFileName = readJLFile(infoFile, &programInfo, &prevKnoVals);
            contentFile = openFile(contentFileName);
        } else {
            printf("File name not found. \n");
            return EXIT_FAILURE;
        }
    } else {
        return EXIT_FAILURE;
    }

    /* Loading a library with specific function */
    printf("Loading library with specific functions\n");
    durationSteps[1][1] = clock();
    /* Row 2 of durationSteps devoted to loading data */
    printf("Loading data\n");
    durationSteps[2][0] = clock();

    /* Now, the matrix with data and that showing the patterns must be consistent: */
    /* Content matrix dimensions */
    matrixInt322DStruct content;
    content.rows = 0;
    content.cols = 0;

    int commas = 0, c;
    while((c = fgetc(contentFile)) != '\n'){
        if(c == ',')
            commas++;
    }
    programInfo.nRoundsInPattern = (commas+1)/3;
    content.cols = programInfo.nRoundsInPattern * 3;
    rewind(contentFile);

    /* When we found a "\n" the rows can be incremented */
    while ((c=fgetc(contentFile)) != EOF) {
        if(c == '\n'){
            content.rows++;
        }
    }

    /* Now, with all the information, we need to start again the file reading to get the data */
    rewind(contentFile);

    int a, b;
    content.data = (int32_t**)malloc(content.rows*sizeof(int32_t *));
    for (a = 0; a < content.rows; a++)
        content.data[a] = (int32_t*)malloc(content.cols*sizeof(int32_t));

    /* Obtaining the file data */
    for (a = 0; a < content.rows; a++){
        for (b = 0; b < content.cols; b++){
            char* line = malloc(sizeof(char)*11);
            getfield(contentFile, line);
            content.data[a][b] = (int32_t)strtol(line, NULL, 16);
            line = NULL;
            free(line);
        }
    }

    fclose(contentFile);

    durationSteps[2][1] = clock();
    printf("Finished the DATA load.\n");

    /* Verify everything is correct */
    /* The address bits must be an integer and less than 32 */
    bool isAddressSizeGood = false;
    if ((isalnum(programInfo.nBitsAddress) == 0) && (programInfo.nBitsAddress > 0) && (programInfo.nBitsAddress < 32)) {
        isAddressSizeGood = true;
    } else {
        printf("The variable nBitsAddress badly defined in the datafile. Skipping analysis. \n");
        EXIT_FAILURE;
    }

    /* This number determinates the rank of the values in the histogram */
    int32_t LN = pow(2,(programInfo.nBitsAddress-programInfo.nBits4Blocks))-1;

    /* A small dictionary to save and quickly access some values */
    vectorIntStruct nAddressesInRound;
    nAddressesInRound.length = programInfo.nRoundsInPattern;
    nAddressesInRound.data = calloc(nAddressesInRound.length, sizeof(int));

    /* The following value is just to initialize some variables.
     * Also used later to determine Shallowness of XORing Rule. */
    int32_t NDVTop = round(0.5*content.rows*(content.rows-1));

    /* The variables are initialized to be used as backups of some critical variables. */
    matrixUint322DStruct xordvbackup;
    xordvbackup.rows = NDVTop;
    xordvbackup.cols = programInfo.nRoundsInPattern;
    xordvbackup.data = calloc(xordvbackup.rows, sizeof(uint32_t*));
    for(a = 0; a < xordvbackup.rows; a++){
    	xordvbackup.data[a] = calloc(xordvbackup.cols, sizeof(uint32_t));
    }

    matrixInt322DStruct posdvbackup;
    posdvbackup.rows = NDVTop;
    posdvbackup.cols = programInfo.nRoundsInPattern;
    posdvbackup.data = calloc(posdvbackup.rows, sizeof(int32_t*));
    for(a = 0; a < posdvbackup.rows; a++){
    	posdvbackup.data[a] = calloc(posdvbackup.cols, sizeof(int32_t));
    }

    matrixUint323DStruct xordvmatrixbackup;
    xordvmatrixbackup.rows = content.rows;
    xordvmatrixbackup.cols = content.rows;
    xordvmatrixbackup.dims = programInfo.nRoundsInPattern;
    xordvmatrixbackup.data = calloc(content.rows, sizeof(uint32_t**));
    for(a = 0; a < content.rows; a++){
    	xordvmatrixbackup.data[a] = calloc(content.rows, sizeof(uint32_t*));
        for(b = 0; b < content.rows; b++){
            xordvmatrixbackup.data[a][b] = calloc(programInfo.nRoundsInPattern, sizeof(uint32_t));
        }
    }
    matrixInt323DStruct posdvmatrixbackup;
    posdvmatrixbackup.rows = content.rows;
    posdvmatrixbackup.cols = content.rows;
    posdvmatrixbackup.dims = programInfo.nRoundsInPattern;
    posdvmatrixbackup.data = calloc(content.rows, sizeof(int32_t**));
    for(a = 0; a < content.rows; a++){
        posdvmatrixbackup.data[a] = calloc(content.rows, sizeof(int32_t*));
        for(b = 0; b < content.rows; b++){
            posdvmatrixbackup.data[a][b] = calloc(programInfo.nRoundsInPattern, sizeof(int32_t));
        }
    }

    /* Check if the dimensions of the data are correct */
    bool isConsistentPatternRawData = false;
    if (content.cols != (3*programInfo.nRoundsInPattern)){
        printf("The number of rounds differs from pattern matrix to raw data matrix. Fix it. \n");
        EXIT_FAILURE;
    } else {
        isConsistentPatternRawData = true;
    }

    /* If both conditions are accomplished, we can start. If they aren't, skip. */
    if (!isConsistentPatternRawData && !isAddressSizeGood){
    	EXIT_FAILURE;
    }
    printf("Data seem to be well formatted. Starting analysis.\n");

    /* First step. It is necessary to extract excessive repetitions in the
     * histogram associated with every test. */
    /* We create a very large matrix used to create the histograms. One matrix
     * for the XOR operation, the other for the positive */
    printf("Creating Histograms...\n");
    /* Row 3 devoted to creating histograms */
    durationSteps[3][0] = clock();

    /* Thus, the first column of the histogram contains the number, the other
     * the number of occurrences in each test. */
    b = 1;
    matrixInt322DStruct xordvhistogram;
    xordvhistogram.rows = LN;
    xordvhistogram.cols = programInfo.nRoundsInPattern + 1;
    xordvhistogram.data = calloc(LN, sizeof(int32_t *));
    for(a = 0; a < xordvhistogram.rows; a++){
    	xordvhistogram.data[a] = calloc(xordvhistogram.cols, sizeof(int32_t));
        xordvhistogram.data[a][0] = b;
        b++;
    }

    /* The PS DV histogram is identically created. */
    b = 1;
    matrixInt322DStruct posdvhistogram;
    posdvhistogram.rows = LN;
    posdvhistogram.cols = programInfo.nRoundsInPattern + 1;
    posdvhistogram.data = calloc(posdvhistogram.rows, sizeof(int32_t *));
    for(a = 0; a < posdvhistogram.rows; a++){
    	posdvhistogram.data[a] = calloc(posdvhistogram.cols, sizeof(int32_t));
        posdvhistogram.data[a][0] = b;
        b++;
    }

    int ktest, elems, k;
    int32_t ndv;
    for (ktest = 1; ktest < programInfo.nRoundsInPattern+1; ++ktest) {
        elems = content.rows;
        ndv = (0.5*content.rows*(content.rows-1));
        /* First, address vector must be isolated. It is always in the col. 3*ktest-3.
         * Also, dummy elements (0xFFFFFFFF) must be removed. I use a special function. */
        uint32_t* columnvector = (uint32_t*)malloc(sizeof(uint32_t)*(elems));
        for (k = 0; k < elems; ++k) {
            columnvector[k] = content.data[k][3*ktest-3];
        }
        uint32_t* addresses = extract_addressvector(columnvector, &elems);
        columnvector = NULL;

        uint32_t* RWcyclesVector = calloc(elems, sizeof(uint32_t));
        for (k = 0; k < elems; k++) {
        	RWcyclesVector[k] = content.data[k][3*ktest-1];
        }
        nAddressesInRound.data[ktest-1] = elems;

        /* Now, let us create the associated DV vectors: XOR and Positive Subtraction in matrix format and as a vector. */
        matrixUint322DStruct xordvmatrix = create_DVmatrix(addresses, elems, "xor", RWcyclesVector, programInfo.nBits4Blocks, LN);
        vectorUint32Struct xordvvector = create_DVvector(&xordvmatrix);
        matrixUint322DStruct posdvmatrix = create_DVmatrix(addresses, elems, "pos", RWcyclesVector, programInfo.nBits4Blocks, LN);
        vectorUint32Struct posdvvector = create_DVvector(&posdvmatrix);

        free(addresses);
        free(RWcyclesVector);

        vectorInt32Struct xordvhistogramVector = create_histogram(&xordvvector, LN);
        vectorInt32Struct posdvhistogramVector = create_histogram(&posdvvector, LN);

        for(a = 0; a < xordvhistogram.rows; a++){
            xordvhistogram.data[a][ktest] = xordvhistogramVector.data[a+1];
        }
        for(a = 0; a < posdvhistogram.rows; a++){
            posdvhistogram.data[a][ktest] = posdvhistogramVector.data[a+1];
        }
        free(xordvhistogramVector.data);
        free(posdvhistogramVector.data);

        /* Other variables must be saved as well. But the procedure is a bit different
         * as the elements do not have identical size. */
        copyOfMatrix(xordvmatrix.data, xordvmatrixbackup.data, elems, ktest);
        copyOfMatrix(posdvmatrix.data, posdvmatrixbackup.data, elems, ktest);
        copyOfVector(xordvvector.data, xordvbackup.data, ndv, ktest);
        copyOfVector(posdvvector.data, posdvbackup.data, ndv, ktest);
        free(xordvvector.data);
        free(posdvvector.data);
        free(xordvmatrix.data);
        free(posdvmatrix.data);
    }

    printf("Completed creation of partial histograms.\n");
    printf("Creating the combined histogram.\n");
    matrixUint322DStruct totalDVHistogram;
    totalDVHistogram.rows = LN;
    totalDVHistogram.cols = 3;
    totalDVHistogram.data = calloc(LN, sizeof(uint32_t*));
    b = 1;
	for (a = 0; a < totalDVHistogram.rows; a++){
		totalDVHistogram.data[a] = calloc(totalDVHistogram.cols, sizeof(uint32_t));
        totalDVHistogram.data[a][0] = b;
        b++;
	}

    if (programInfo.nRoundsInPattern == 1) { /* Only one experiment, total its not necessary. */
        liberaUint(totalDVHistogram.data, LN);
    }

    vectorIntStruct arrayMaxXORValue;
    arrayMaxXORValue.length = programInfo.nRoundsInPattern;
    arrayMaxXORValue.data = calloc(arrayMaxXORValue.length, sizeof(int));
    vectorIntStruct arrayMaxPOSValue;
    arrayMaxPOSValue.length = programInfo.nRoundsInPattern;
    arrayMaxPOSValue.data = calloc(arrayMaxPOSValue.length, sizeof(int));

	if (programInfo.nRoundsInPattern > 1){
		/* First column for LN, Second column for XOR, Third for subtraction.*/
		int kvalues;
		for (kvalues = 0; kvalues < LN; kvalues++){
            /* Sum of all the columns on the histogram of a selected row */
			int32_t sumXor = 0, sumPos = 0;
			for (b = 1; b < programInfo.nRoundsInPattern+1; b++) {
                if (xordvhistogram.data[kvalues][b] > arrayMaxXORValue.data[b-1]) {
                    arrayMaxXORValue.data[b-1] = xordvhistogram.data[kvalues][b];
                }
                if (posdvhistogram.data[kvalues][b] > arrayMaxPOSValue.data[b-1]) {
                    arrayMaxPOSValue.data[b-1] = posdvhistogram.data[kvalues][b];
                }
				sumXor += xordvhistogram.data[kvalues][b];
				sumPos += posdvhistogram.data[kvalues][b];
			}
			totalDVHistogram.data[kvalues][1] = sumXor;
			totalDVHistogram.data[kvalues][2] = sumPos;
			sumXor = 0;
			sumPos = 0;

		}
		printf("Done.\n");
    }

    durationSteps[3][1] = clock();

    printf("Counting repetitions in partial histograms...\n");
    /* Row 5 devoted to counting events */
    durationSteps[4][0] = clock();

    /* We need to get the maximum value for the vector creation  */
    int maxXORValue = 0, maxPOSValue = 0;
    for (int i = 0; i < programInfo.nRoundsInPattern; i++) {
        if (arrayMaxXORValue.data[i] > maxXORValue) {
            maxXORValue = arrayMaxXORValue.data[i];
        }
    }
    for (int i = 0; i < programInfo.nRoundsInPattern; i++) {
        if (arrayMaxPOSValue.data[i] > maxPOSValue) {
            maxPOSValue = arrayMaxPOSValue.data[i];
        }
    }

    matrixInt322DStruct xordvrepetitions;
    xordvrepetitions.rows = maxXORValue + 1;
    xordvrepetitions.cols = programInfo.nRoundsInPattern + 1;
    xordvrepetitions.data = calloc(xordvrepetitions.rows, sizeof(int32_t*));
    for(a = 0; a < xordvrepetitions.rows; a++){
        xordvrepetitions.data[a] = calloc(xordvrepetitions.cols, sizeof(int32_t));
    }

    matrixInt322DStruct posdvrepetitions;
    posdvrepetitions.rows = maxPOSValue + 1;
    posdvrepetitions.cols = programInfo.nRoundsInPattern + 1;
    posdvrepetitions.data = calloc(posdvrepetitions.rows, sizeof(int32_t*));
    for(a = 0; a < posdvrepetitions.rows; a++){
        posdvrepetitions.data[a] = calloc(posdvrepetitions.cols, sizeof(int32_t));
    }

	for (int i = 1; i < xordvrepetitions.rows; i++){
		xordvrepetitions.data[i][0] = i;
	}
	for (int i = 1; i < posdvrepetitions.rows; i++){
		posdvrepetitions.data[i][0] = i;
	}

    /* Let us evaluate the repetitions in the histogram. A little confusing, perhaps,
     * but true. Every column xordvhistogram[:,ktest+1] is an independent histogram
     * and counts return the number of times elements between 0 and the maximum
     * value in the histogram appears. Everything is saved in a temporary variable
     * to save some additional zeros in the matrix of the repetitions. */
    for (ktest = 1; ktest < programInfo.nRoundsInPattern+1; ktest++){
        vectorInt32Struct xordvrepetitionstmp;
        vectorInt32Struct posdvrepetitionstmp;

        xordvrepetitionstmp = countsInt(&xordvhistogram, ktest, arrayMaxXORValue.data[ktest-1]);
        posdvrepetitionstmp = countsInt(&posdvhistogram, ktest, arrayMaxPOSValue.data[ktest-1]);

        free(xordvrepetitionstmp.data);
        free(posdvrepetitionstmp.data);
    }

    printf("Completed task\n.");

    int maxTotalXORValue = 0, maxTotalPOSValue = 0;
    for (a = 0; a < totalDVHistogram.rows; a++) {
        if (totalDVHistogram.data[a][1] > maxTotalXORValue) {
            maxTotalXORValue = totalDVHistogram.data[a][1];
        }
        if (totalDVHistogram.data[a][2] > maxTotalPOSValue) {
            maxTotalPOSValue = totalDVHistogram.data[a][2];
        }
    }

    matrixInt322DStruct xorDVtotalrepetitions;
    xorDVtotalrepetitions.rows = maxTotalXORValue + 1;
    xorDVtotalrepetitions.cols = 2;
	xorDVtotalrepetitions.data = calloc(xorDVtotalrepetitions.rows, sizeof(int32_t *));
	for (a = 0; a < xorDVtotalrepetitions.rows; a++){
		xorDVtotalrepetitions.data[a] = calloc(xorDVtotalrepetitions.cols, sizeof(int32_t));
	}

    matrixInt322DStruct posDVtotalrepetitions;
    posDVtotalrepetitions.rows = maxTotalPOSValue + 1;
    posDVtotalrepetitions.cols = 2;
    posDVtotalrepetitions.data = calloc(posDVtotalrepetitions.rows, sizeof(int32_t *));
	for (a = 0; a < posDVtotalrepetitions.rows; a++){
		posDVtotalrepetitions.data[a] = calloc(posDVtotalrepetitions.cols, sizeof(int32_t));
	}

    if (programInfo.nRoundsInPattern > 1){
        printf("Checking Total Histogram...\n");
        /* The first column is filled with numbers from 0 to max */
        for (int i = 0; i < xorDVtotalrepetitions.rows; i++){
			xorDVtotalrepetitions.data[i][0] = i;
		}

        vectorInt32Struct xordvtotaltemp = countsUint(&totalDVHistogram, 1, maxTotalXORValue);
        for (a = 0; a < xorDVtotalrepetitions.rows; a++) {
            xorDVtotalrepetitions.data[a][1] = xordvtotaltemp.data[a];
        }
        free(xordvtotaltemp.data);

        /* The first column is filled with numbers from 0 to max */
		for (int i = 0; i < posDVtotalrepetitions.rows; i++){
			posDVtotalrepetitions.data[i][0] = i;
		}

        vectorInt32Struct posdvtotaltemp = countsUint(&totalDVHistogram, 2, maxTotalPOSValue);
        for (a = 0; a < posDVtotalrepetitions.rows; a++) {
            posDVtotalrepetitions.data[a][1] = posdvtotaltemp.data[a];
        }
        free(posdvtotaltemp.data);

        printf("Created statistics for combinations.\n");
    }

    durationSteps[4][1] = clock();
    printf("STARTING TO LOOK UP ANOMALOUSLY REPEATED VALUES...\n");
    durationSteps[5][0] = clock();
    int XORANOMALS = 0;

    matrixInt322DStruct XORextracted_values00 = extractAnomalDVSelfConsistency("xor", &xorDVtotalrepetitions, &totalDVHistogram, &nAddressesInRound, &xordvhistogram, &xordvmatrixbackup, &XORANOMALS);
    matrixInt322DStruct POSextracted_values00 = extractAnomalDVSelfConsistency("pos", &posDVtotalrepetitions, &totalDVHistogram, &nAddressesInRound, &posdvhistogram, &posdvmatrixbackup, &XORANOMALS);

    durationSteps[5][1] = clock();

    /* Thus, selfconsistency is completed. Now, it is time to apply the TRACE rule. */
    printf("APPLYING THE TRACE RULE FOR XOR DV SET...\n");
    vectorInt32Struct XORextracted_TraceRule = traceRule(&xorDVtotalrepetitions, &totalDVHistogram, &XORextracted_values00, LN);

    matrixInt322DStruct XORextracted_values01;
    XORextracted_values01.rows = XORextracted_values00.rows + XORextracted_TraceRule.length;
    XORextracted_values01.cols = 2;
    XORextracted_values01.data = calloc(XORextracted_values01.rows, sizeof(int32_t*));
    for (a = 0; a < XORextracted_values01.rows; a++) {
        XORextracted_values01.data[a] = calloc(XORextracted_values01.cols, sizeof(int32_t));
    }

    int index = 0;
    for (a = 0; a < XORextracted_values00.rows; a++) {
        for (b = 0; b < 2; b++) {
            XORextracted_values01.data[a][b] = XORextracted_values00.data[a][b];
        }
        index++;
    }
    for (a = 0; a < XORextracted_TraceRule.length; a++) {
        XORextracted_values01.data[index][0] = XORextracted_TraceRule.data[a];
        XORextracted_values01.data[index][1] = totalDVHistogram.data[XORextracted_TraceRule.data[a]-1][1];
        index++;
    }

    matrixInt322DStruct POSextracted_values01;
    POSextracted_values01.rows = POSextracted_values00.rows;
    POSextracted_values01.cols = 2;
    POSextracted_values01.data = calloc(POSextracted_values01.rows, sizeof(int32_t*));
    for (a = 0; a < POSextracted_values01.rows; a++) {
        POSextracted_values01.data[a] = calloc(2, sizeof(int32_t));
    }
    for (a = 0; a < POSextracted_values01.rows; a++) {
        for (b = 0; b < 2; b++) {
            POSextracted_values01.data[a][b] = POSextracted_values00.data[a][b];
        }
    }

    /** USING DV VALUES INTRODUCED BY USER. **/
    printf("\nRECYCLING PREVIOUSLY KNOWN DV ELEMENTS...");

    vectorInt32Struct UnknowXORValues = setdiff(&prevKnoVals.XORextracted_values, &XORextracted_values01, 0);
    matrixInt322DStruct opXOR;
    opXOR.rows = UnknowXORValues.length;
    opXOR.cols = 2;
    opXOR.data = calloc(opXOR.rows, sizeof(int32_t*));
    for (int i = 0; i < opXOR.rows; i++) {
        opXOR.data[i] = calloc(opXOR.cols, sizeof(int32_t));
    }
    for (int i = 0; i < UnknowXORValues.length; i++) {
        opXOR.data[i][0] = totalDVHistogram.data[UnknowXORValues.data[i]][0];
        opXOR.data[i][1] = totalDVHistogram.data[UnknowXORValues.data[i]][1];
    }
    matrixInt322DStruct XORextracted_values02 = vcat(&XORextracted_values01, &opXOR);

    vectorInt32Struct UnknowPOSValues = setdiff(&prevKnoVals.POSextracted_values, &POSextracted_values01, 0);
    matrixInt322DStruct opPOS;
    opPOS.rows = UnknowPOSValues.length;
    opPOS.cols = 2;
    opPOS.data = calloc(opPOS.rows, sizeof(int32_t*));
    for (int i = 0; i < opPOS.rows; i++) {
        opPOS.data[i] = calloc(opPOS.cols, sizeof(int32_t));
    }
    for (int i = 0; i < UnknowPOSValues.length; i++) {
        opPOS.data[i][0] = totalDVHistogram.data[UnknowPOSValues.data[i]][0];
        opPOS.data[i][1] = totalDVHistogram.data[UnknowPOSValues.data[i]][2];
    }
    matrixInt322DStruct POSextracted_values02 = vcat(&POSextracted_values01, &opPOS);

    int nUnknownXORValues = UnknowXORValues.length, nUnknownPOSValues = UnknowPOSValues.length;

    printf("\n\tXOR: ");
    if (nUnknownXORValues > 0){
        printf("%d unknown value", nUnknownXORValues);
        if (nUnknownXORValues != 1){
            printf("s");
            }
        printf(" added.");
    } else {
        printf("Nothing to include.");
    }
    printf("\n\tPositive Subtraction: ");
    if (nUnknownPOSValues > 0) {
        printf("%d unknown value", nUnknownPOSValues);
        if (nUnknownPOSValues != 1) {
            printf("s");
        }
        printf(" added.");
    } else {
        printf("Nothing to include.\n");
    }
    printf(" ");

    /* Time for Exploring MCUs */
    printf("\nSEARCHING INSIDE MCUs (FIRST PASS)...");

    matrixInt322DStruct XORextracted_values03;
    matrixInt322DStruct POSextracted_values03;

    vectorInt32Struct discoveredXORDvs;
    vectorInt32Struct discoveredPOSDvs;

    extractAnomalDVfromClusters(&content, &XORextracted_values02, &POSextracted_values02, &xorDVtotalrepetitions, &posDVtotalrepetitions, &totalDVHistogram, &xordvmatrixbackup, &posdvmatrixbackup, &nAddressesInRound, LN, &discoveredXORDvs, &discoveredPOSDvs, &XORextracted_values03, &POSextracted_values03);

    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");

    /* XORing RULE */
    printf("\n\nAPPLYING XORING RULE...");
    vectorInt32Struct XORDVfromXORing;

    matrixInt322DStruct XORextracted_values04 = criticalXORvaluesFromXORingRule(&XORextracted_values03, &xorDVtotalrepetitions, &totalDVHistogram, &XORDVfromXORing);

    /** Time for Exploring MCUs **/
    printf("\nSEARCHING INSIDE MCUs (SECOND PASS)...");

    matrixInt322DStruct XORextracted_values05;
    matrixInt322DStruct POSextracted_values05;
    vectorInt32Struct discoveredXORDvs2;
    discoveredXORDvs2.length = 0;
    vectorInt32Struct discoveredPOSDvs2;
    discoveredPOSDvs2.length = 0;

    extractAnomalDVfromClusters(&content, &XORextracted_values04, &POSextracted_values03, &xorDVtotalrepetitions, &posDVtotalrepetitions, &totalDVHistogram, &xordvmatrixbackup, &posdvmatrixbackup, &nAddressesInRound, LN, &discoveredXORDvs2, &discoveredPOSDvs2, &XORextracted_values05, &POSextracted_values05);

    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");
    durationSteps[0][1] = clock();

    printf(" ");

    printf("\nMAKING THE FINAL ORGANIZATION OF ADDRESSES FOR PRESENTATION...");
    matrixInt323DStruct cmbRes = propose_MCUs("cmb", &xordvmatrixbackup, &posdvmatrixbackup, &XORextracted_values05, &POSextracted_values03, &nAddressesInRound);

    matrixInt322DStruct finalResult = condensate_summary(&cmbRes, &nAddressesInRound);

    printf("\tEnded.");

    printf("*******************************************************************");


    printf("\nPRESENTING RESULTS:\n");
    printf("\nWARNING: This program is not prepared to deal with Experiments with MBUs!!!!!!!!");


    FILE* fileResults = NULL;
    fileResults = fopen("Results.txt","w");
    if(fileResults == NULL){
    	printf("\n Error creating the finalResult archive for QT");
    }

    matrixInt322DStruct mbuSummary = locate_mbus(&content, programInfo.pattern, programInfo.nRoundsInPattern, programInfo.nWordWidth);

    if((mbuSummary.rows * mbuSummary.cols) > 0){
        printf("\n\tThe following MBUs have been observed...\n");
        printf("\n THE PRESENCE OF MBUs MAKES THE FOLLOWING RESULTS INACCURATE.");
        printf("\n\tFor your information, the expected average numbers of several SBUs in a cell is ");
        printf("\n");

        for(int i = 0; i < programInfo.nRoundsInPattern; i++){
            printf("Test %d: %d \t", i+1, 7*nAddressesInRound.data[i]*(nAddressesInRound.data[i]-1)/16);
        }
        printf("\n");
    } else {
    	printf("\n\t---> Fortunately, your data lacks MBUs and following results are not commited.\n");
    }

    matrixInt323DStruct hexCMBRES;
    hexCMBRES = index2address(&cmbRes, &content);

    printf("\nKINDS OF EVENTS");
    int i,j,z;

	char *hex2 = calloc(255, sizeof(char));
    char *prt = calloc(255, sizeof(char));

    int32_t test = programInfo.nRoundsInPattern;
    char *num = calloc(255, sizeof(char));
    sprintf(num, "%d", test);
    fputs(num, fileResults);
    fputs("\n", fileResults);
    free(num);
    int32_t filas = hexCMBRES.rows;
    char *num1 = calloc(255, sizeof(char));
    sprintf(num1, "%d", test);
    fputs(num1, fileResults);
    fputs("\n", fileResults);
    free(num1);

    for(i = 0; i< programInfo.nRoundsInPattern;i++){
    	printf("\nTest %d: %d SBUs", i+1, finalResult.data[0][i+1]);
        sprintf(prt, "%d", finalResult.data[0][i+1]);
        fputs(prt, fileResults);
        fputs("\n", fileResults);
    	for(j = 1; j < finalResult.rows; j++){
            if(finalResult.data[j][i+1] !=0){
                int32_t print;
    			printf("\n \t%d %d-bit MCU", finalResult.data[j][i+1], j+1);
                char *num = calloc(1, sizeof(char));
                sprintf(num, "%d", j+1);
                fputs(num, fileResults);
                fputs(" ", fileResults);
				print = finalResult.data[j][i+1];
				sprintf(hex2, "%d", print);
    			fputs(hex2, fileResults);

    			if(finalResult.data[j][i+1]!=1){
    				printf("s");
    			}
			}else{
                    fputs("-", fileResults);
					printf("			");

            }
            fputs("\n", fileResults);
    	}

    	int32_t sum = 0;
    	for(z= 0 ;z < finalResult.rows;z++){
    		sum+= (finalResult.data[z][0]*finalResult.data[z][i+1]);
    	}
    	printf("\n Affected Addresses: %d",sum);
    	sprintf(hex2, "%d", sum);
        fputs("s", fileResults);
    	fputs(hex2, fileResults);
        fputs("\n", fileResults);
        fputs("*\n", fileResults);
    }

    fputs("*****\n", fileResults);
    printf("\n\n*******************************************************************");
    printf("\nSUMMARY OF THE PROPOSED ANOMALOUS DV VALUES:");
    printf("\nSelf-consistence after XORing: ");

    for(i = 0; i < XORextracted_values00.rows; i++){
    	printf("\n 0x%x", XORextracted_values00.data[i][0]);
    }
    printf("\n---------------------------");

    printf("\n Self-consistence after Positive Subtraction: ");
    for(i = 0; i < POSextracted_values00.rows; i++){
    	printf("\n%d", POSextracted_values00.data[i][0]);
    	printf(" (0x%x )\n", POSextracted_values00.data[i][0]);
    }
     printf("\n---------------------------");

    printf("\n Trace Rule: ");
    for(i = 0; i < XORextracted_TraceRule.length; i++){
    	printf("\n 0x%x", XORextracted_TraceRule.data[i]);
    }
    printf("\n---------------------------");
    printf("\nAdded anomalous DV values: ");
    printf("\n\tXOR operation:");

    for (i = 0; i < UnknowXORValues.length; i++) {
        printf("\n0x%x", UnknowXORValues.data[i]);
    }

    printf("\n ");
    printf("\n\tPOS. SUB. operation:");

    for (i = 0; i < UnknowPOSValues.length; i++) {
        printf("\n0x%x", UnknowPOSValues.data[i]);
    }

    printf("\n--------------------------");
    printf("\nMCU for XORing (First Pass): ");

    for(i = 0; i < discoveredXORDvs.length; i++){
    	k = discoveredXORDvs.data[i];
    	printf("\n0x%x", k);
    }

    printf("\n--------------------------");
    printf("\nMCU for Pos. Sub. (First Pass):");
    for(i = 0; i < discoveredPOSDvs.length; i++){
    	k = discoveredPOSDvs.data[i];
    	printf("\n %d (0x%x)", discoveredPOSDvs.data[i],k);
    }

    printf("\n--------------------------");
    printf("\nXORING RULE: ");
    for(i = 0; i < XORDVfromXORing.length; i++){
    	printf("\n0x%x", XORDVfromXORing.data[i]);
    }

    printf("\n--------------------------");
    printf("\nMCU for XORing (Second Pass): ");
    for(i = 0; i < discoveredXORDvs2.length; i++){
    	k = discoveredXORDvs2.data[i];
    	printf("\n0x%x", k);
    }

    printf("\n--------------------------");
    printf("\nMCU for Pos. Sub. (Second Pass): ");
    for(i = 0; i < discoveredPOSDvs2.length; i++){
    	k = discoveredPOSDvs2.data[i];
    	printf("\n %d (0x%x)", discoveredPOSDvs2.data[i], k);
    }

    printf("\n-------------------------");
    printf("\n-------------------------");

    printf("\nElapsed Time: %f s\n", ((float)(durationSteps[0][1] - durationSteps[0][0]) / 1000000.0F ));

    time_t tiempo = time(0);
    struct tm *tlocal = localtime(&tiempo);
    char output[128];
    char archivo [128];
    strftime(output, 128,"%Y-%m-%dT%H:%M:%S", tlocal);
    char* str = "AffectedAddresses-";
    strcpy(archivo, str);
    strcat(archivo, output);

    FILE* FileForAddresses = fopen(archivo, "w");
    if (FileForAddresses == NULL) {
        printf("File not created. \n");
    }

    fputs("ADDRESSES INVOLVED IN MCUs\n", FileForAddresses);
    int rows = 0, cols = 0;
    for(ktest = 0; ktest < programInfo.nRoundsInPattern; ktest++){
        fputs( "\n\tTest Number: ", FileForAddresses);
        char *test = calloc(3, sizeof(char));
        sprintf(test, "%d", ktest+1);
        fputs(test, FileForAddresses);
        free(test);
        fputs("\n\t------------------\n", FileForAddresses);
        for( rows = 0; rows < hexCMBRES.rows; rows++){
            if(hexCMBRES.data[rows][0][ktest] != 0){
                for(cols = 0; cols < hexCMBRES.cols;cols++){
                    int32_t toBePrinted;
                    toBePrinted = (hexCMBRES.data[rows][cols][ktest]);
                    if(toBePrinted != 0){
                        int32_t storeWord;
                        storeWord = content.data[cmbRes.data[rows][cols][ktest]][3*ktest+1];
                        fputs( "\t0x", FileForAddresses);
                        char *hex = calloc(255, sizeof(char));
                        sprintf(hex, "%.*x",6, toBePrinted);
                        fputs(hex, FileForAddresses);
                        char *hex1 = calloc(255, sizeof(char));
                        sprintf(hex1, " (0x%.*x)",2, storeWord);
                        fputs(hex1, FileForAddresses);
                        free(hex);
                        free(hex1);
                    }else{
                        fputs("\t, ", FileForAddresses);
                    }
                }
                fputs("\n", FileForAddresses);
            }else{
                break;
            }
        }
        fputs("\n------------------\n", FileForAddresses);
    }

    for(ktest = 0; ktest < programInfo.nRoundsInPattern; ktest++){
        char *test = calloc(3, sizeof(char));
        sprintf(test, "%d", ktest+1);
        fputs(test, fileResults);
        free(test);
        fputs("\n", fileResults);
        for(rows = 0; rows < hexCMBRES.rows; rows++){
            if(hexCMBRES.data[rows][0][ktest] != 0){
                for(cols = 0; cols < hexCMBRES.cols;cols++){
                    int32_t toBePrinted;
                    toBePrinted = (hexCMBRES.data[rows][cols][ktest]);
                    if(toBePrinted != 0){
                        int32_t storeWord;
                        storeWord = content.data[cmbRes.data[rows][cols][ktest]][3*ktest+1];
                        fputs( "\t0x", fileResults);
                        char *hex = calloc(255, sizeof(char));
                        sprintf(hex, "%.*x", 6, toBePrinted);
                        fputs(hex, fileResults);
                        free(hex);
                    }else{
                        fputs("\t- ", fileResults);
                    }
                }
                fputs("\n", fileResults);
            }else{
                break;
            }
        }
    }
    fclose(fileResults);
    fclose(FileForAddresses);
    return 0;
}
