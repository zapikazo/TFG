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

//TimeMessages POR DEFINIR

int main(){

    bool** v1 = calloc(3, sizeof(bool*));
    for (int i = 0; i < 3; i++) {
        v1[i] = calloc(3, sizeof(bool));
        v1[i][0] = true;
        v1[i][1] = true;
        v1[i][2] = false;
    }
    int32_t* v2 = calloc(5, sizeof(int32_t));
    int32_t* v3 = calloc(3, sizeof(int32_t));

    
    
    
   // int32_t* lol  = create_histogram(v2, 5, 10);
        
    for (int i = 0; i < 10; i++) {
        //for (int j = 0; j < 3; j++) {
           // printf("%d", lol[i]);
        //}
        //printf("\n");
    }

    /* Matrix which save times involved in the different stages to detect bottlenecks */
    float durationSteps[11][2] = {0};
    
    /* Row 0 of durationSteps devoted to the total program time */
    durationSteps[0][0] = time(NULL);
    /* Row 1 of durationSteps devoted to loading libraries */
    durationSteps[1][0] = time(NULL);
    /* Loading a library with specific function */
    printf("Loading library with specific functions\n");
    durationSteps[1][1] = time(NULL);
    /* Row 2 of durationSteps devoted to loading data */
    printf("Loading data\n");
    durationSteps[2][0] = time(NULL);
    
    /* File opening */
    FILE* file;
    file = fopen("CY62167_090nm_TS.csv", "r");
    if (file == NULL) {
        printf("ERROR OPENING FILE.\n");
        return 0;
    }

   /* const Pattern = [
                     1 0x00
                     2 0x55
                     3 0xFF
                     ];*/

    //OBTENER EL TAMANIO DE LOS DATOS, ARCHIVO JL
    int nBitsAddress = 0;
    int nBits4Blocks = 0;
    int NWordWidth = 0;
    int nRoundsInPattern;
    int commas = 0, c;
    
    /* Now, the matrix with data and that showing the patterns must be consistent: */
    /* Content matrix dimensions */
    int rawDataMatrixNRows = 0;
    int rawDataMatrixNCols;
    
    while((c = fgetc(file)) != '\n'){
        if(c == ',')
            commas++;
    }
    nRoundsInPattern = (commas+1)/3;
    rawDataMatrixNCols = nRoundsInPattern * 3;
    rewind(file);
    while((fgetc(file)) != ','){
        nBitsAddress++;
    }
    //nBitsAddress = (nBitsAddress-2)*4;
    //SEGUN EL JL SON 21, necesitamos la apertura de dos ficheros
    nBitsAddress = 21;
    nBits4Blocks = 1;
    NWordWidth = 8;
    // Leer jl para obtener patterns
    int32_t** pattern = calloc(nRoundsInPattern, sizeof(int32_t*));
    for (int i = 0; i < nRoundsInPattern; i++) {
        pattern[i] = calloc(2, sizeof(int32_t));
        pattern[i][0] = i+1;
    }
    pattern[0][1] = 0x00;
    pattern[1][1] = 0x55;
    pattern[2][3] = 0xFF;

    /* When we found a "\n" the rows can be incremented */
    while ((c=fgetc(file)) != EOF) {
        if(c == '\n'){
            rawDataMatrixNRows++;
        }
    }
    /* Now, with all the information, we need to start again the file reading to get the data */
    rewind(file);
    
    int a, b;
    int32_t **content;
    content = (int32_t **)malloc(rawDataMatrixNRows*sizeof(int32_t *));
    for (a = 0; a < rawDataMatrixNRows; a++)
        content[a] = (int32_t *)malloc(rawDataMatrixNCols*sizeof(int32_t));
    
    /* Obtaining the file data */
    for (a = 0; a < rawDataMatrixNRows; a++){
        for (b = 0; b < rawDataMatrixNCols; b++){
            char* line = malloc(sizeof(char)*11);
            getfield(file, line);
            content[a][b] = (int32_t)strtol(line, NULL, 16);
            line = NULL;
            free(line);
        }
    }
   /* for (a = 0; a < rawDataMatrixNRows; a++){
        printf("\n");
        for (b = 0; b < rawDataMatrixNCols; b++){
            printf(" %#x", content[a][b] );
        }
    }*/

    durationSteps[2][1] = time(NULL);
    printf("Finished the DATA load.\n");

    /* Verify everything is correct */
    /* The address bits must be an integer and less than 32 */
    bool isAddressSizeGood = false;
    if ((isalnum(nBitsAddress) == 0) && (nBitsAddress > 0) && (nBitsAddress < 32)) {
        isAddressSizeGood = true;
    } else {
        printf("Upsss: The variable nBitsAddress badly defined in the datafile. Skipping analysis.\n");
        EXIT_FAILURE;
    }

    int32_t LN = pow(2,(nBitsAddress-nBits4Blocks))-1;

    /* A small dictionary to save and quickly access some values: Addresses in rounds */
    int* nAddressesInRound;
    nAddressesInRound = calloc(nRoundsInPattern, sizeof(int));
    
    /* I will use the following value just to initialize some variables.
     * Also used later to determine Shallowness of XORing Rule. */
    int32_t NDVTop = round(0.5*rawDataMatrixNRows*(rawDataMatrixNRows-1));
    
    /* The variables are initialized to be used as backups of some critical variables. */
    //xordvbackup
    uint32_t **xordvbackup = (uint32_t**)malloc(NDVTop*sizeof(uint32_t*));
    for(a = 0; a < NDVTop; a++){
    	xordvbackup[a] = (uint32_t*)malloc(nRoundsInPattern*sizeof(uint32_t));
    }
    //posdvbackup
    int32_t **posdvbackup = (int32_t**)malloc(NDVTop*sizeof(int32_t*));
    for(a = 0; a < NDVTop; a++){
    	posdvbackup[a] = (int32_t*)malloc(nRoundsInPattern*sizeof(int32_t));
    }
    
    //xordvmatrixbackup
    uint32_t*** xordvmatrixbackup = calloc(rawDataMatrixNRows, sizeof(uint32_t**));
    for(a = 0; a < rawDataMatrixNRows; a++){
    	xordvmatrixbackup[a] = calloc(rawDataMatrixNRows, sizeof(uint32_t*));
        for(b = 0; b < rawDataMatrixNRows; b++){
            xordvmatrixbackup[a][b]= calloc(nRoundsInPattern, sizeof(uint32_t));
        }
    }

    //posdvmatrixbackup
    int32_t*** posdvmatrixbackup = calloc(rawDataMatrixNRows, sizeof(int32_t**));
    for(a = 0; a < rawDataMatrixNRows; a++){
        posdvmatrixbackup[a] = calloc(rawDataMatrixNRows, sizeof(int32_t*));
        for(b = 0; b < rawDataMatrixNRows; b++){
            posdvmatrixbackup[a][b] = calloc(nRoundsInPattern, sizeof(int32_t));
        }
    }

    bool isConsistentPatternRawData = false;
    if (rawDataMatrixNCols != (3*nRoundsInPattern)){
        printf("Mmmmmm: The number of rounds differs from pattern matrix to raw data matrix. Fix it. Bye.\n");
        EXIT_FAILURE;
    } else {
        isConsistentPatternRawData = true;
    }
    
    /*If both conditions are accomplished, we can start. If they aren't, skip. */
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
    durationSteps[3][0]=time(NULL);

    /* Thus, the first column of the histogram contains the number, the other
     * the number of occurrences in each test. */
    b = 1;
    int32_t **xordvhistogram = (int32_t**)calloc(LN, sizeof(int32_t *));
    for(a = 0; a < LN; a++){
    	xordvhistogram[a] = (int32_t*)calloc((nRoundsInPattern+1), sizeof(int32_t));
        xordvhistogram[a][0] = b;
        b++;
    }

    /* The PS DV histogram is identically created. */
    b = 1;
    int32_t **posdvhistogram = (int32_t**)calloc(LN, sizeof(int32_t *));
    for(a = 0; a < LN; a++){
    	posdvhistogram[a] = (int32_t*)calloc((nRoundsInPattern+1), sizeof(int32_t));
        posdvhistogram[a][0] = b;
        b++;
    }
	
    int ktest, elems, k;
    int32_t ndv;
    for (ktest = 1; ktest < nRoundsInPattern+1; ++ktest) {
        elems = rawDataMatrixNRows;
        ndv = (0.5*rawDataMatrixNRows*(rawDataMatrixNRows-1));
        /* First, address vector must be isolated. It is always in the col. 3*ktest-3.
         * Also, dummy elements (0xFFFFFFFF) must be removed. I use a special function. */
        uint32_t* columnvector = (uint32_t*)malloc(sizeof(uint32_t)*(elems));
        for (k = 0; k < elems; ++k) {
            columnvector[k] = content[k][3*ktest-3];
        }
        uint32_t* addresses = extract_addressvector(columnvector, &elems);
        columnvector = NULL;
        free(columnvector);

        uint32_t* RWcyclesVector = calloc(elems, sizeof(uint32_t));
        for (k = 0; k < elems; k++) {
        	RWcyclesVector[k] = content[k][3*ktest-1];
        }
        nAddressesInRound[ktest-1] = elems;

        /* Now, let us create the associated DV vectors: XOR and Positive Subtraction
         * in matrix format and as a vector. */
        uint32_t** xordvmatrix;
        uint32_t** posdvmatrix;
        uint32_t* xordvvector;
        uint32_t* posdvvector;
        int xordvvectortam = 0;
        int posdvvectortam = 0;
        xordvmatrix = create_DVmatrix(addresses, elems, "xor", RWcyclesVector, nBits4Blocks, LN);
        xordvvector = create_DVvector(xordvmatrix, elems, &xordvvectortam);
        posdvmatrix = create_DVmatrix(addresses, elems, "pos", RWcyclesVector, nBits4Blocks, LN);
        posdvvector = create_DVvector(posdvmatrix, elems, &posdvvectortam);
        
        free(addresses);
        free(RWcyclesVector);
        
        //Rellenar la columna del exp con for
        int32_t* xordvhistogramVector;
        int32_t* posdvhistogramVector;
        
        xordvhistogramVector = create_histogram(xordvvector, xordvvectortam, LN);
        posdvhistogramVector = create_histogram(posdvvector, posdvvectortam, LN);
        
       /* printf("xordvhistogramVector");
        for (a = 1; a < 7; a++) {
            printf("%d ", xordvhistogramVector[a]);
        }
        printf("\n");*/
        
        //Meter valores del vector en la matriz de histogramas, en la posición del cero está el 1
        for(a = 0; a < LN; a++){
            xordvhistogram[a][ktest] = xordvhistogramVector[a+1];
        }
        for(a = 0; a < LN; a++){
            posdvhistogram[a][ktest] = posdvhistogramVector[a+1];
        }
        free(xordvhistogramVector);
        free(posdvhistogramVector);
        
        /* Other variables must be saved as well. But the procedure is a bit different
         * as the elements do not have identical size. */
        copyOfMatrix(xordvmatrix, xordvmatrixbackup, elems, ktest);
        copyOfMatrix(posdvmatrix, posdvmatrixbackup, elems, ktest);
        copyOfVector(xordvvector, xordvbackup, ndv, ktest);
        copyOfVector(posdvvector, posdvbackup, ndv, ktest);
        free(xordvvector);
        free(posdvvector);
        free(xordvmatrix);
        free(posdvmatrix);
    }
    
    printf("Completed creation of partial histograms.\n");
    printf("Creating the combined histogram.\n");
	uint32_t** totalDVHistogram = calloc(LN, sizeof(uint32_t*));
    b = 1;
	for (a = 0; a < LN; a++){
		totalDVHistogram[a] = calloc(3, sizeof(uint32_t));
        totalDVHistogram[a][0] = b;
        b++;
	}
    if (nRoundsInPattern == 1) { // Único experimento no es necesario el total. Esto no sé si
        liberaUint(totalDVHistogram, LN);
    }
    
    int* arrayMaxXORValue = calloc(nRoundsInPattern, sizeof(int));
    int* arrayMaxPOSValue = calloc(nRoundsInPattern, sizeof(int));
    
	if (nRoundsInPattern > 1){
		/* First column for LN, Second column for XOR, Third for subtraction.*/
		int kvalues,valor=0;
		for (kvalues = 0; kvalues < LN; kvalues++){
            /* Sum of all the columns on the histogram of a selected row */
			int32_t sumXor = 0, sumPos = 0;
			for (b = 1; b < nRoundsInPattern+1; b++) {
                if (xordvhistogram[kvalues][b] > arrayMaxXORValue[b-1]) {
                    arrayMaxXORValue[b-1] = xordvhistogram[kvalues][b];
                }
                if (posdvhistogram[kvalues][b] > arrayMaxPOSValue[b-1]) {
                    arrayMaxPOSValue[b-1] = posdvhistogram[kvalues][b];
                }
				sumXor += xordvhistogram[kvalues][b];
				sumPos += posdvhistogram[kvalues][b];
			}
			totalDVHistogram[kvalues][1] = sumXor;
			totalDVHistogram[kvalues][2] = sumPos;
			sumXor = 0;
			sumPos = 0;

		}
		printf("Done.\n");
    }
    
    /*for (int i = 0; i < nRoundsInPattern; i++) {
        printf("%d ", arrayMaxXORValue[i]);
    }
    printf("\n");
    for (int i = 0; i < nRoundsInPattern; i++) {
        printf("%d ", arrayMaxPOSValue[i]);
    }
    printf("\n");*/

    durationSteps[3][1] = time(NULL);
            
    printf("Counting repetitions in partial histograms...\n");
    /* Row 5 devoted to counting events */
    durationSteps[4][0] = time(NULL);
    
    // Obtención de valor mayor de repeticiones para la creación de los vectores
    int maxXORValue = 0, maxPOSValue = 0;
    for (int i = 0; i < nRoundsInPattern; i++) {
        if (arrayMaxXORValue[i] > maxXORValue) {
            maxXORValue = arrayMaxXORValue[i];
        }
    }
    for (int i = 0; i < nRoundsInPattern; i++) {
        if (arrayMaxPOSValue[i] > maxPOSValue) {
            maxPOSValue = arrayMaxPOSValue[i];
        }
    }

    //xordvrepetitions
    int32_t** xordvrepetitions = (int32_t **)calloc((maxXORValue+1), sizeof(int32_t*));
    for(a = 0; a < maxXORValue+1; a++){
        xordvrepetitions[a] = (int32_t *)calloc((nRoundsInPattern+1), sizeof(int32_t));
    }

    //posdvrepetitions
    int32_t** posdvrepetitions = (int32_t **)calloc((maxPOSValue+1), sizeof(int32_t*));
    for(a = 0; a < maxPOSValue+1; a++){
        posdvrepetitions[a] = (int32_t *)calloc((nRoundsInPattern+1), sizeof(int32_t));
    }

	for (int i = 1; i < maxXORValue+1; i++){
		xordvrepetitions[i][0] = i;
	}
	for (int i = 1; i < maxPOSValue+1; i++){
		posdvrepetitions[i][0] = i;
	}

    /* Let us evaluate the repetitions in the histogram. A little confusing, perhaps,
     * but true. Every column xordvhistogram[:,ktest+1] is an independent histogram
     * and counts return the number of times elements between 0 and the maximum
     * value in the histogram appears. Everything is saved in a temporary variable
     * to save some additional zeros in the matrix of the repetitions. */
            
    for (ktest = 1; ktest < nRoundsInPattern+1; ktest++){
        int i;
        int32_t* xordvrepetitionstmp;
        int32_t* posdvrepetitionstmp;

        xordvrepetitionstmp = counts(xordvhistogram, LN, ktest, arrayMaxXORValue[ktest-1]);
        posdvrepetitionstmp = counts(posdvhistogram, LN, ktest, arrayMaxPOSValue[ktest-1]);

        free(xordvrepetitionstmp);
        free(posdvrepetitionstmp);
        
        
 /*       for (i = 0; i < maxXORValue+1; i++) {
            printf("%d ", xordvrepetitions[i][ktest]);
        }
        printf("\n");
        for (i = 0; i < maxPOSValue+1; i++) {
            printf("%d ", posdvrepetitions[i][ktest]);
        }
        printf("\n");*/
        
    }
    
    printf("Completed task\n.");
    
    int maxTotalXORValue = 0, maxTotalPOSValue = 0;
    for (a = 0; a < LN; a++) {
        if (totalDVHistogram[a][1] > maxTotalXORValue) {
            maxTotalXORValue = totalDVHistogram[a][1];
        }
        if (totalDVHistogram[a][2] > maxTotalPOSValue) {
            maxTotalPOSValue = totalDVHistogram[a][2];
        }
    }
    
	//xorDVtotalrepetitions 
	int32_t** xorDVtotalrepetitions = (int32_t**)calloc((maxTotalXORValue + 1), sizeof(int32_t *));
	for (a = 0; a < maxTotalXORValue + 1; a++){
		xorDVtotalrepetitions[a] = (int32_t*)calloc(2, sizeof(int32_t));
	}
	//posDVtotalrepetitions
	int32_t** posDVtotalrepetitions = (int32_t**)calloc((maxTotalPOSValue + 1), sizeof(int32_t*));
	for (a = 0; a < maxTotalPOSValue + 1; a++){
		posDVtotalrepetitions[a] = (int32_t*)calloc(2, sizeof(int32_t));
	}
    if (nRoundsInPattern > 1){
        printf("Checking Total Histogram...\n");
        // Se rellena la primera columna con valores del 0 hasta el máximo
        for (int i = 0; i < maxTotalXORValue+1; i++){
			xorDVtotalrepetitions[i][0] = i;
		}

        int32_t* xordvtotaltemp = counts(totalDVHistogram, LN, 1, maxTotalXORValue);
        for (a = 0; a < maxTotalXORValue+1; a++) {
            xorDVtotalrepetitions[a][1] = xordvtotaltemp[a];
        }
        free(xordvtotaltemp);

        // Se rellena la primera columna con valores del 0 hasta el máximo
		for (int i = 0; i < maxTotalPOSValue+1; i++){
			posDVtotalrepetitions[i][0] = i;
		}

        int32_t* posdvtotaltemp = counts(totalDVHistogram, LN, 2, maxTotalPOSValue);
        for (a = 0; a < maxTotalPOSValue+1; a++) {
            posDVtotalrepetitions[a][1] = posdvtotaltemp[a];
        }
        free(posdvtotaltemp);

        printf("Created statistics for combinations.\n");
    }
                    
    durationSteps[4][1] = time(NULL);
    printf("STARTING TO LOOK UP ANOMALOUSLY REPEATED VALUES...\n");
    durationSteps[5][0] = time(NULL);
    int XORextracted_values00Lenght = 0, POSextracted_values00Lenght = 0;
    int XORANOMALS = 0;

    int32_t** XORextracted_values00 = extractAnomalDVSelfConsistency("xor", xorDVtotalrepetitions, (maxTotalXORValue + 1), totalDVHistogram, LN, LN, nAddressesInRound, nRoundsInPattern, xordvhistogram, xordvmatrixbackup, rawDataMatrixNRows, &XORANOMALS, &XORextracted_values00Lenght);
    int32_t** POSextracted_values00 = extractAnomalDVSelfConsistency("pos", posDVtotalrepetitions, (maxTotalPOSValue + 1), totalDVHistogram, LN, LN, nAddressesInRound, nRoundsInPattern, posdvhistogram, posdvmatrixbackup, rawDataMatrixNRows, &XORANOMALS, &POSextracted_values00Lenght);
    
    durationSteps[5][1] = time(NULL);
    
    /* Thus, selfconsistency is completed. Now, it is time to apply the TRACE rule. */
    printf("APPLYING THE TRACE RULE FOR XOR DV SET...\n");
    
    int traceRuleLong = 0;
    int32_t* XORextracted_TraceRule = traceRule(xorDVtotalrepetitions, (maxTotalXORValue + 1), totalDVHistogram, XORextracted_values00, XORextracted_values00Lenght, LN, &traceRuleLong);

    //XORextracted_values00 y el pos, necesitamos guardarlos, no se borran
    int32_t XORextracted_values01Lenght = XORextracted_values00Lenght + traceRuleLong;
    int32_t** XORextracted_values01 = calloc(XORextracted_values01Lenght, sizeof(int32_t*));
    for (a = 0; a < XORextracted_values01Lenght; a++) {
        XORextracted_values01[a] = calloc(2, sizeof(int32_t));
    }
    //XORextracted_values01=vcat(XORextracted_values00, TotalDVhistogram[XORextracted_TraceRule, 1:2])
    // Esto es el vcat, función?
    int index = 0;
    for (a = 0; a < XORextracted_values00Lenght; a++) {
        for (b = 0; b < 2; b++) {
            XORextracted_values01[a][b] = XORextracted_values00[a][b];
        }
        index++;
    }
    for (a = 0; a < traceRuleLong; a++) {
        XORextracted_values01[index][0] = XORextracted_TraceRule[a];
        XORextracted_values01[index][1] = totalDVHistogram[XORextracted_TraceRule[a]-1][1];
        index++;
    }
    

     /*for (a = 0; a < XORextracted_values01Lenght; a++) {
         for (b = 0; b < 2; b++) {
             printf("%d ",XORextracted_values01[a][b]);
         }
         printf("\n");
     }
     printf("\n");*/
    
    //POSextracted_values01=POSextracted_values00, copia simple de matrices, hacer función?
    int32_t POSextracted_values01Lenght =  POSextracted_values00Lenght;
    int32_t** POSextracted_values01 = calloc(POSextracted_values01Lenght, sizeof(int32_t*));
    for (a = 0; a < POSextracted_values01Lenght; a++) {
        POSextracted_values01[a] = calloc(2, sizeof(int32_t));
    }
    for (a = 0; a < POSextracted_values01Lenght; a++) {
        for (b = 0; b < 2; b++) {
            POSextracted_values01[a][b] = POSextracted_values00[a][b];
        }
    }
    
    
    /** USING DV VALUES INTRODUCED BY USER. **/
   // printf("\nRECYCLING PREVIOUSLY KNOWN DV ELEMENTS...");
    

 /*   printf("\nRECYCLING PREVIOUSLY KNOWN DV ELEMENTS...");
    //UnknowXORValues = setdiff(PreviouslyKnownXOR, XORextracted_values01[:,1])
    //XORextracted_values02=vcat(XORextracted_values01, TotalDVhistogram[UnknowXORValues, 1:2])
    //UnknowPOSValues = setdiff(PreviouslyKnownPOS, POSextracted_values01[:,1])
    //POSextracted_values02=vcat(POSextracted_values01, TotalDVhistogram[UnknowPOSValues, [1,3]])
    int nUnknownXORValues, nUnknownPOSValues;
    //NUnknownXORValues = length(UnknowXORValues)
    //NUnknownPOSValues = length(UnknowPOSValues)
    
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
    printf(" ");*/

// Time for Exploring MCUs
    printf("\nSEARCHING INSIDE MCUs (FIRST PASS)...");

    int32_t** XORextracted_values03 = NULL;
    int XORextracted_values03Length = 0;
    int32_t** POSextracted_values03 = NULL;
    int POSextracted_values03Length = 0;
    /*XORextracted_values03 = calloc(XORextracted_values01Lenght, sizeof(int32_t*));
    for (a = 0; a < XORextracted_values01Lenght; a++) {
        XORextracted_values03[a] = calloc(2, sizeof(int32_t));
    }

    POSextracted_values03 = calloc(POSextracted_values01Lenght, sizeof(int32_t*));
    for (a = 0; a < POSextracted_values01Lenght; a++) {
        POSextracted_values03[a] = calloc(2, sizeof(int32_t));
    }*/
    int32_t* discoveredXORDvs = NULL;
    //int32_t* discoveredXORDvs = calloc(rawDataMatrixNRows, sizeof(int32_t));
    int discoveredXORDVsLength = 0;
    int32_t* discoveredPOSDvs = NULL;
    //int32_t* discoveredPOSDvs = calloc(rawDataMatrixNRows, sizeof(int32_t));
    int discoveredPOSDVsLength = 0;

    extractAnomalDVfromClusters(content, rawDataMatrixNRows, rawDataMatrixNCols, XORextracted_values01, XORextracted_values01Lenght, 2, POSextracted_values01, POSextracted_values01Lenght, 2, xorDVtotalrepetitions, (maxTotalXORValue + 1), posDVtotalrepetitions, totalDVHistogram, xordvmatrixbackup, rawDataMatrixNRows, rawDataMatrixNRows, nRoundsInPattern, posdvmatrixbackup, nAddressesInRound, LN, &discoveredXORDvs, &discoveredXORDVsLength, &discoveredPOSDvs, &discoveredPOSDVsLength, &XORextracted_values03, &XORextracted_values03Length, &POSextracted_values03, &POSextracted_values03Length);

   /* printf("\n");
    printf("\n Discoveredxordvs");
    for (a = 0; a < discoveredXORDVsLength; a++) {
        printf("%d ", discoveredXORDvs[a]);
    }
    printf("\n");

  /*  printf("\n Discoveredposdvs");
    for (a = 0; a < discoveredPOSDVsLength; a++) {
        printf("%d ", discoveredPOSDvs[a]);
    }
    printf("\n");

    printf("\n xorextractedvalues03");
    for (a = 0; a < XORextracted_values03Length; a++) {
        for (b = 0; b < 2; b++) {
            printf("%d ", XORextracted_values03[a][b]);
        }
        printf("\n");
    }
    printf("\n");

    printf("\n posextractedvalues03");
    for (a = 0; a < POSextracted_values03Length; a++) {
        for (b = 0; b < 2; b++) {
            printf("%d ", POSextracted_values03[a][b]);
        }
        printf("\n");
    }
    printf("\n");*/


    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");
    //POSDVsMCU1=[];
    //POSextracted_values03=POSextracted_values02;
// XORing RULE
    printf("\n\nAPPLYING XORING RULE...");
    int32_t** XORextracted_values04 =  NULL;
    int XORextracted_values04Length = 0;
    int32_t* XORDVfromXORing = NULL;
    int XORDVfromXORingLength = 0;

    criticalXORvaluesFromXORingRule(XORextracted_values03, XORextracted_values03Length, xorDVtotalrepetitions, (maxTotalXORValue + 1), totalDVHistogram, LN, &XORextracted_values04, &XORextracted_values04Length, &XORDVfromXORing, &XORDVfromXORingLength);
    //POSextracted_values04=POSextracted_values03

   /* printf("\n");
    for (int i = 0; i < XORextracted_values04Length; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%d ", XORextracted_values04[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < XORDVfromXORingLength; i++) {
        printf("%d ", XORDVfromXORing[i]);
    }

    printf("\n");*/


/** Time for Exploring MCUs **/
    printf("\nSEARCHING INSIDE MCUs (SECOND PASS)...");
    /*XORDVsMCU2, POSDVsMCU2, XORextracted_values05,
    POSextracted_values05=ExtractAnomalDVfromClusters(  Content[:,1:3:end],
                                                      XORextracted_values04,
                                                      POSextracted_values04,
                                                      xorDVtotalrepetitions,
                                                      posDVtotalrepetitions,
                                                      TotalDVhistogram,
                                                      xordvmatrixbackup,
                                                      posdvmatrixbackup,
                                                      LN,
                                                      RandomnessThreshold)*/

    int32_t** XORextracted_values05 = NULL;
    int XORextracted_values05Length = 0;
    int32_t** POSextracted_values05 = NULL;
    int POSextracted_values05Length = 0;
    int32_t* discoveredXORDvs2 = NULL;
    int discoveredXORDVs2Length = 0;
    int32_t* discoveredPOSDvs2 = NULL;
    int discoveredPOSDVs2Length = 0;
    extractAnomalDVfromClusters(content, rawDataMatrixNRows, rawDataMatrixNCols, XORextracted_values04, XORextracted_values04Length, 2, POSextracted_values03, POSextracted_values03Length, 2, xorDVtotalrepetitions, (maxTotalXORValue + 1), posDVtotalrepetitions, totalDVHistogram, xordvmatrixbackup, rawDataMatrixNRows, rawDataMatrixNRows, nRoundsInPattern, posdvmatrixbackup, nAddressesInRound, LN, &discoveredXORDvs2, &discoveredXORDVs2Length, &discoveredPOSDvs2, &discoveredPOSDVs2Length, &XORextracted_values05, &XORextracted_values05Length, &POSextracted_values05, &POSextracted_values05Length);


   /*  printf("\n");
     printf("\n Discoveredxordvs");
     for (a = 0; a < discoveredXORDVs2Length; a++) {
     printf("%d ", discoveredXORDvs2[a]);
     }
     printf("\n");

       printf("\n Discoveredposdvs");
     for (a = 0; a < discoveredPOSDVs2Length; a++) {
     printf("%d ", discoveredPOSDvs[a]);
     }
     printf("\n");

     printf("\n xorextractedvalues03");
     for (a = 0; a < XORextracted_values05Length; a++) {
     for (b = 0; b < 2; b++) {
     printf("%d ", XORextracted_values05[a][b]);
     }
     printf("\n");
     }
     printf("\n");

     printf("\n posextractedvalues03");
     for (a = 0; a < POSextracted_values05Length; a++) {
     for (b = 0; b < 2; b++) {
     printf("%d ", POSextracted_values05[a][b]);
     }
     printf("\n");
     }
     printf("\n");*/



    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");
    //POSDVsMCU2=[];
    //POSextracted_values05=POSextracted_values04;
    durationSteps[0][1]=time(NULL);
    
    printf(" ");
    
    printf("\nMAKING THE FINAL ORGANIZATION OF ADDRESSES FOR PRESENTATION...");
    int cmbResRows, cmbResCols, cmbResDims;
    int32_t*** cmbRes = propose_MCUs("cmb", xordvmatrixbackup, rawDataMatrixNRows, rawDataMatrixNRows, nRoundsInPattern, posdvmatrixbackup, XORextracted_values05, XORextracted_values05Length, POSextracted_values03, POSextracted_values03Length, nAddressesInRound, &cmbResRows, &cmbResCols, &cmbResDims);

    int finalResultRows = 0, finalResultCols = 0;
    int32_t** finalResult = condensate_summary(cmbRes, cmbResRows, cmbResCols, cmbResDims, nAddressesInRound, &finalResultRows, &finalResultCols);

  /*  for (a = 0; a < cmbResRows; a++) {
        for (b = 0; b < cmbResCols; b++) {
            printf("%d ", cmbRes[a][b][0]);
        }
        printf("\n");
    }
    printf("\n");

    for (a = 0; a < cmbResRows; a++) {
        for (b = 0; b < cmbResCols; b++) {
            printf("%d ", cmbRes[a][b][1]);
        }
        printf("\n");
    }
    printf("\n");

    for (a = 0; a < cmbResRows; a++) {
        for (b = 0; b < cmbResCols; b++) {
            printf("%d ", cmbRes[a][b][2]);
        }
        printf("\n");
    }
    printf("\n");

    for (a = 0; a < finalResultRows; a++) {
        for (b = 0; b < finalResultCols; b++) {
            printf("%d ", finalResult[a][b]);
        }
        printf("\n");
    }
    printf("\n");*/


    printf("\tEnded.");
    
    printf("*******************************************************************");
    printf("*******************************************************************");


    printf("\nPRESENTING RESULTS:\n");
    printf("\nWARNING: This program is not prepared to deal with Experiments with MBUs!!!!!!!!");

    int mbuSummaryRows = 0, mbuSummaryCols = 0;
    int32_t** mbuSummary = locate_mbus(content, rawDataMatrixNRows, rawDataMatrixNCols, pattern, nRoundsInPattern, NWordWidth, &mbuSummaryRows, &mbuSummaryCols);


    for (a = 0; a < mbuSummaryRows; a++) {
        for (b = 0; b < mbuSummaryCols; b++) {
            printf("%d ", mbuSummary[a][b]);
        }
        printf("\n");
    }
    printf("\n");

    printf("\n EXIT!");

    if((mbuSummaryRows*mbuSummaryCols)>0)
        printf("\n\tThe following MBUs have been observed...\n");

	/*for (a = 0; a < mbuSummaryRows; a++) {
		printf("\t\t Address %d in Test %d (%#): bitflips in the word.\n",mbuSummary[a,1],mbuSummary[a,0],mbuSummary[a,2],mbuSummary[a,3]);
	}*/
        /*for kmbu=1:length(MBUSummary[:,1])
            print("\t\t Address ",  MBUSummary[kmbu, 2], " in Test ", MBUSummary[kmbu, 1],
                  " (0x", hex(MBUSummary[kmbu, 3],6)"): ", MBUSummary[kmbu, 4], " bitflips in the word.\n")
            end
            println("\nTHE PRESENCE OF MBUs MAKES THE FOLLOWING RESULTS INACCURATE.\n")
            println("\tFor you information, the expected average numbers of several SBUs in a cell is ")
            print("\t\t")
            for k=1:NRoundsInPattern
                print("Test ", k, ": ", 7*NAddressesInRound[k]*(NAddressesInRound[k]-1)/16,"\t")
                end
                println("\n")

                else

                    println("\t---> Fortunately, your data lacks MBUs and following results are not commited.\n")

                    end*/




    return 0;
}
