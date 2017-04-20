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

    
    v2[0] = 2;
    v2[1] = 1;
    v2[2] = 2;
    v2[3] = 3;
    v2[4] = 2;
    
    v3[0] = 0;
    v3[1] = 2;
    v3[2] = 1;
    
    int result = 0;

    printf("%d", lookForElem(v2, 5, 0));
    printf("%d", lookForElem(v2, 5, 1));
    printf("%d", lookForElem(v2, 5, 2));
    printf("%d", lookForElem(v2, 5, 3));
    printf("%d", lookForElem(v2, 5, 4));
    printf("\n");
    
    
   // int32_t* lol  = create_histogram(v2, 5, 10);
        
    for (int i = 0; i < 10; i++) {
        //for (int j = 0; j < 3; j++) {
           // printf("%d", lol[i]);
        //}
        //printf("\n");
    }

    /* Matrix which save times involved in the different stages to detect bottlenecks */
    float durationSteps[11][2] = {0}; //EL VALOR A COPIAR AQUI F-1 y C-1
    
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
    
    //ABRIR FICHERO
    FILE* file;
    file = fopen("CY62167_090nm_TS.csv", "r");
    if (file == NULL) {
        printf("ERROR OPENING FILE.\n");
        return 0;
    }
    //OBTENER EL TAMANIO DE LOS DATOS
    int nBitsAddress = 0;
    int nBits4Blocks = 0;
    int nRoundsInPattern; // NUM DE EXPERIMENTOS
    int commas = 0, c;
    
    /* Now, the matrix with data and that showing the patterns must be consistent: */
    //Dimensiones de matriz content
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


    while ((c=fgetc(file)) != EOF) {
        if(c == '\n'){
            rawDataMatrixNRows++;
        }
    }
    //LEER EL FICHERO PARA GUARDARLO EN LA MATRIZ
    rewind(file);
    
    int a, b;
    int32_t **content;
    content = (int32_t **)malloc(rawDataMatrixNRows*sizeof(int32_t *));
    for (a = 0; a < rawDataMatrixNRows; a++)
        content[a] = (int32_t *)malloc(rawDataMatrixNCols*sizeof(int32_t));
    
    // Rellenamos la matriz content con los datos de fichero
    for (a = 0; a < rawDataMatrixNRows; a++){
        for (b = 0; b < rawDataMatrixNCols; b++){
            char* line = malloc(sizeof(char)*11);
            getfield(file, line);
            content[a][b] = (int32_t)strtol(line, NULL, 16);
            line = NULL;
            free(line);
        }
    }
    for (a = 0; a < rawDataMatrixNRows; a++){
        printf("\n");
        for (b = 0; b < rawDataMatrixNCols; b++){
            printf(" %#x", content[a][b] );
        }
    }

    
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
     * the number of occurrences in each test (Test K --> Col. K +1). */
    b = 1;
    int32_t **xordvhistogram = (int32_t**)calloc(LN, sizeof(int32_t *));
    for(a = 0; a < LN; a++){
    	xordvhistogram[a] = (int32_t*)calloc((nRoundsInPattern+1), sizeof(int32_t));
        //la columna cero ponemos numeros del 1 al LN
        xordvhistogram[a][0] = b;
        b++;
    }

    /* The PS DV histogram is identically created. */
    b = 1;
    int32_t **posdvhistogram = (int32_t**)calloc(LN, sizeof(int32_t *));
    for(a = 0; a < LN; a++){
    	posdvhistogram[a] = (int32_t*)calloc((nRoundsInPattern+1), sizeof(int32_t));
        //la columna cero ponemos numeros del 1 al LN
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

        
        //Se debe convertir de decimal a hexadecimal el vector RWCycles?
        uint32_t* RWcyclesVector = calloc(elems, sizeof(uint32_t));
        for (k = 0; k < elems; k++) {
        	RWcyclesVector[k] = content[k][3*ktest-1];
            printf("%#x\n", RWcyclesVector[k]);
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

        //Por aquí asegurarse del pos, 
        
        free(addresses);
        free(RWcyclesVector);
        
        //Rellenar la columna del exp con for
        int32_t* xordvhistogramVector;
        int32_t* posdvhistogramVector;
        
        xordvhistogramVector = create_histogram(xordvvector, xordvvectortam, LN);
        posdvhistogramVector = create_histogram(posdvvector, posdvvectortam, LN);
        
        /*printf("xordvhistogramVector");
        for (a = 0; a < LN; a++) {
            printf("%d ", xordvhistogramVector[a]);
        }
        printf("\n");*/
        
        //Meter valores del vector en la matriz de histogramas
        for(a = 0; a < LN; a++){
            xordvhistogram[a][ktest+1] = xordvhistogramVector[a];
        }
        for(a = 0; a < LN; a++){
            posdvhistogram[a][ktest+1] = posdvhistogramVector[a];
        }
        free(xordvhistogramVector);
        free(posdvhistogramVector);
        
        
        for (int i = 0; i < 4; i++) {
            printf("%d \n", xordvhistogram[0][i]);
        }
        for (int i = 0; i < 4; i++) {
            printf("%d \n", xordvhistogram[1][i]);
        }
        for (int i = 0; i < 4; i++) {
            printf("%d \n", xordvhistogram[2][i]);
        }
        
       /* for (a = 0; a < LN; a++) {
            for (b = 0; nRoundsInPattern+1; b++) {
                printf(" %#x", xordvhistogram[a][b] );
            }
        }*/
        
        /* Other variables must be saved as well. But the procedure is a bit different
         * as the elements do not have identical size. */
        copyOfMatrix(xordvmatrix, xordvmatrixbackup, elems, ktest);
        copyOfMatrix(posdvmatrix, posdvmatrixbackup, elems, ktest);
        copyOfVector(xordvvector, xordvbackup, ndv, ktest);
        copyOfVector(posdvvector, posdvbackup, ndv, ktest); //WARNING
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
	int maxXORValue = 0, maxPOSValue = 0;
    if (nRoundsInPattern == 1) { // Único experimento no es necesario el total.
        liberaUint(totalDVHistogram, LN);
    }
	if (nRoundsInPattern > 1){
		/* PRIMERA RANGO LN, Second column for XOR, Third for subtraction.*/
		int kvalues;
		for (kvalues = 0; kvalues < LN; kvalues++){
			// SUMA DE TODAS LAS COLUMNAS DE LOS HISTOGRAMAS QUE CORRESPONDAN A UN LN
			int32_t sumXor = 0, sumPos = 0;
			for (b = 1; b < nRoundsInPattern+1; b++) {
				sumXor += xordvhistogram[kvalues][b];
				sumPos += posdvhistogram[kvalues][b];
			}
			if (sumXor > maxXORValue) {
				maxXORValue = sumXor;
			}
			if (sumPos > maxPOSValue) {
				maxPOSValue = sumPos;
			}
			totalDVHistogram[kvalues][1] = sumXor;
			totalDVHistogram[kvalues][2] = sumPos;
			sumXor = 0;
			sumPos = 0;
		}
		printf("Done.\n");
    }
    
    durationSteps[3][1] = time(NULL);
            
    printf("Counting repetitions in partial histograms...\n");
    /* Row 5 devoted to counting events */
    durationSteps[4][0] = time(NULL);
    
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
            
    for (ktest = 1; ktest < nRoundsInPattern+1; ++ktest){
        int32_t xordvrepetitionstmp = countsOfElems(xordvhistogram, maxXORValue, ktest, LN);
		int32_t posdvrepetitionstmp = countsOfElems(posdvhistogram, maxPOSValue, ktest, LN);
        
        xordvrepetitions[0][ktest] = xordvrepetitionstmp;
        posdvrepetitions[0][ktest]= posdvrepetitionstmp;
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
        for (int i = 1; i < maxTotalXORValue+1; i++){
			xorDVtotalrepetitions[i][0] = i;
		}
        int32_t xorCounts = countsOfElems((int32_t**)totalDVHistogram, maxTotalXORValue, 1, LN);
        for (a = 0; a < maxTotalXORValue; a++) {
            xorDVtotalrepetitions[a][1] = xorCounts;
        }

        // Se rellena la primera columna con valores del 0 hasta el máximo
		for (int i = 1; i < maxTotalPOSValue+1; i++){
			posDVtotalrepetitions[i][0] = i;
		}
        // Se rellena la segunda columna con la aparición de esos valores en el histograma
        int32_t posCounts = countsOfElems((int32_t**)totalDVHistogram, maxTotalPOSValue, 2, LN);
        for (a = 0; a < maxTotalPOSValue; a++) {
            posDVtotalrepetitions[a][1] = posCounts;
        }
        
        printf("Created statistics for combinations.\n");
    }
                    
    durationSteps[4][1] = time(NULL);
    printf("STARTING TO LOOK UP ANOMALOUSLY REPEATED VALUES...\n");
    durationSteps[5][0] = time(NULL);
    int XORextracted_values00Lenght = 0, POSextracted_values00Lenght = 0;
    
    //int32_t** XORextracted_values00 = extractAnomalDVSelfConsistency("xor", xorDVtotalrepetitions, (maxTotalXORValue + 1), totalDVHistogram, LN, LN, nRoundsInPattern, xordvhistogram, xordvmatrixbackup, rawDataMatrixNRows, &XORextracted_values00Lenght);
    //int32_t** POSextracted_values00 = extractAnomalDVSelfConsistency("pos", posDVtotalrepetitions, (maxTotalPOSValue + 1), totalDVHistogram, LN, LN, nRoundsInPattern, posdvhistogram, posdvmatrixbackup, rawDataMatrixNRows, &POSextracted_values00Lenght);

    
    durationSteps[5][1] = time(NULL);
    
    /* Thus, selfconsistency is completed. Now, it is time to apply the TRACE rule. */
    printf("APPLYING THE TRACE RULE FOR XOR DV SET...\n");
    
    int traceRuleLong = 0;
    //int32_t* XORextracted_TraceRule = traceRule(xorDVtotalrepetitions, (maxXORValue + 1), totalDVHistogram, XORextracted_values00, XORextracted_values00Lenght, LN, &traceRuleLong);
    
    //XORextracted_values01=vcat(XORextracted_values00, TotalDVhistogram[XORextracted_TraceRule, 1:2])
    //POSextracted_values01=POSextracted_values00
    
    
    
    /** USING DV VALUES INTRODUCED BY USER. **/
    
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
  /*  XORDVsMCU1, POSDVsMCU1, XORextracted_values03,
    POSextracted_values03=ExtractAnomalDVfromClusters(  Content[:,1:3:end],
                                                      XORextracted_values02,
                                                      POSextracted_values02,
                                                      xorDVtotalrepetitions,
                                                      posDVtotalrepetitions,
                                                      TotalDVhistogram,
                                                      xordvmatrixbackup,
                                                      posdvmatrixbackup,
                                                      LN,
                                                      RandomnessThreshold)*/
    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");
    //POSDVsMCU1=[];
    //POSextracted_values03=POSextracted_values02;
// XORing RULE
    printf("\n\nAPPLYING XORING RULE...");
    /*XORDVfromXORing,  XORextracted_values04=CriticalXORvaluesFromXORingRule( XORextracted_values03,
                                                                            xorDVtotalrepetitions,
                                                                            TotalDVhistogram, LN,RandomnessThreshold)*/
    //POSextracted_values04=POSextracted_values03

// Time for Exploring MCUs
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
    printf("\n\tWARNING: MCUs IN POSITIVE SUBTRACION IN QUARANTINE.\n");
    //POSDVsMCU2=[];
    //POSextracted_values05=POSextracted_values04;
    durationSteps[0][1]=time(NULL);
    
    printf(" ");
    
    printf("\nMAKING THE FINAL ORGANIZATION OF ADDRESSES FOR PRESENTATION...");
   /* CMBRES,XORRES, POSRES=propose_MCUs(xordvmatrixbackup,
                                       posdvmatrixbackup,
                                       XORextracted_values05[:,1],
                                       POSextracted_values05[:,1],
                                       NAddressesInRound,
                                       UnrealMCUsize)*/
    //FinalResult = condensate_summary(CMBRES, Content)
    printf("\tEnded.");
    
    printf("*******************************************************************");
    printf("*******************************************************************");








/*Para matrices bidimensionales
 *
 *   for(x=0; x<5; x++)
 *   	for(y=0; y<5; y++)
    	free(matrix[x][y]);
     for(z=0; z<5; z++)
    	free(matrix[z]);
  	 free(matrix)
 *
 * Para matrices tridimensionales
 *
 *   for(x=0; x<5; x++)
    	free(matrix[x]);
  	 free(matrix)
 *
 *
 *
 *
 *
 *
 *
 *
 */
    
    //00:30
    
    return 0;
}
