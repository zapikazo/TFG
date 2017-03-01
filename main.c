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
    
    
   /* tAddress* columnVector = malloc(sizeof(tAddress)*17);
    
    
    for (i = 0; i < elems; i++){
     printf("%s\n", columnVector[i]);
     }
    //flipped_bits(columnVector[0], columnVector[2], sizeof(tAddress));
    */
    
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
    file = fopen("CY62167_090nm_TSS.csv", "r");
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
    nBitsAddress = (nBitsAddress-2)*4;
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
    /*for (a = 0; a < rawDataMatrixNRows; a++){
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
    uint32_t ***xordvmatrixbackup =(uint32_t***)malloc((rawDataMatrixNRows)*sizeof(uint32_t**));
    for(a = 0; a < (rawDataMatrixNRows); a++)
    	xordvmatrixbackup[a] =(uint32_t**)malloc((rawDataMatrixNRows)* sizeof(uint32_t*));
    for(a = 0; a < (rawDataMatrixNRows+); a++)
     for(b = 0; b < nRoundsInPattern; b++)
    	 xordvmatrixbackup[a][b]= (uint32_t*)malloc(nRoundsInPattern*sizeof(uint32_t));

    //posdvmatrixbackup
    int32_t ***posdvmatrixbackup = (int32_t***)malloc((rawDataMatrixNRows)* sizeof(int32_t**));
    for(a = 0; a < (rawDataMatrixNRows); a++)
    	posdvmatrixbackup[a] = (int32_t**)malloc((rawDataMatrixNRows)* sizeof(int32_t*));
    for(a = 0; a < (rawDataMatrixNRows); a++)
     for(b = 0; b < nRoundsInPattern; b++)
    	 posdvmatrixbackup[a][b] = (int32_t*)malloc(nRoundsInPattern*sizeof(int32_t));

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
    uint32_t **xordvhistogram = (uint32_t**)calloc(LN, sizeof(uint32_t *));
    for(a = 0; a < LN; a++){
    	xordvhistogram[a] = (uint32_t*)calloc((nRoundsInPattern+1), sizeof(uint32_t));
    }
    /* The PS DV histogram is identically created. */
    int32_t **posdvhistogram = (int32_t**)calloc(LN, sizeof(int32_t *));
    for(a = 0; a < LN; a++){
    	posdvhistogram[a] = (int32_t*)calloc((nRoundsInPattern+1), sizeof(int32_t));
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

        uint32_t* RWcyclesVector = (uint32_t*)malloc(sizeof(uint32_t)*(elems));
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
        xordvmatrix = create_DVmatrix(addresses, elems, "xor", RWcyclesVector, nBits4Blocks, (2^20-1));
        xordvvector = create_DVvector(xordvmatrix, elems);
        posdvmatrix = create_DVmatrix(addresses, elems, "pos", RWcyclesVector, nBits4Blocks, (2^20-1));
        posdvvector = create_DVvector(posdvmatrix, elems);
        
        //Rellenar la columna del exp con for
        int32_t* xordvhistogramVector;
        int32_t* posdvhistogramVector;
        
        xordvhistogramVector = create_histogram(xordvvector, ndv, LN);
        posdvhistogramVector = create_histogram(posdvvector, ndv, LN);
        
        //Meter valores del vector en la matriz de histogramas
        for(a = 0; a < LN; a++){
            xordvhistogram[a][ktest] = xordvhistogramVector[a];
        }
        for(a = 0; a < LN; a++){
            posdvhistogram[a][ktest] = posdvhistogramVector[a];
        }
        
        /* Other variables must be saved as well. But the procedure is a bit different
         * as the elements do not have identical size. */
        copyOfMatrix(xordvmatrix, xordvmatrixbackup, elems, ktest);
        copyOfMatrix(posdvmatrix, posdvmatrixbackup, elems, ktest);
        copyOfVector(xordvvector, xordvbackup, ndv, ktest);
        copyOfVector(posdvvector, posdvbackup, ndv, ktest); //WARNING
    }
    
    printf("Completed creation of partial histograms.\n");
    printf("Creating the combined histogram.\n");
    uint32_t** totalDVHistogram = (uint32_t**)malloc(LN*sizeof(uint32_t*));
    for(a = 0; a < LN; a++){
        totalDVHistogram[a] = (uint32_t *)malloc(3*sizeof(uint32_t));
    }
    int maxXORValue = 0, maxPOSValue = 0;

    if (nRoundsInPattern > 1){
        /* PRIMERA RANGO LN, Second column for XOR, Third for subtraction.*/
        int kvalues;
        for (kvalues = 0; kvalues < LN; ++kvalues){
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
    int32_t** xordvrepetitions = (int32_t **)malloc((maxXORValue+1)*sizeof(int32_t *));
    for(a = 0; a < maxXORValue+1; a++){
        xordvrepetitions[a] = (int32_t *)malloc((nRoundsInPattern+1)*sizeof(int32_t));
    }
    b = 0;
    for (a = 0; a < maxXORValue+1; a++) {
        xordvrepetitions[a][0] = b;
        b++;
    }    
    //posdvrepetitions
    int32_t** posdvrepetitions = (int32_t **)malloc((maxPOSValue+1)*sizeof(int32_t *));
    for(a = 0; a < maxPOSValue+1; a++){
        posdvrepetitions[a] = (int32_t *)malloc((nRoundsInPattern+1)*sizeof(int32_t));
    }
    b = 0;
    for (a = 0; a < maxPOSValue+1; a++) {
        posdvrepetitions[a][0] = b;
        b++;
    }
            
    //xordvrepetitions[:,1]=collect(0:1:length(xordvrepetitions[:,1])-1)
    //posdvrepetitions[:,1]=collect(0:1:length(posdvrepetitions[:,1])-1)
            
    /* Let us evaluate the repetitions in the histogram. A little confusing, perhaps,
     * but true. Every column xordvhistogram[:,ktest+1] is an independent histogram
     * and counts return the number of times elements between 0 and the maximum
     * value in the histogram appears. Everything is saved in a temporary variable
     * to save some additional zeros in the matrix of the repetitions. */
            
    for (ktest = 1; ktest < nRoundsInPattern+1; ++ktest){
        //xordvrepetitionstmp=counts(xordvhistogram[:,ktest+1], 0:maximum(xordvhistogram[:,ktest+1]))
        //posdvrepetitionstmp=counts(posdvhistogram[:,ktest+1], 0:maximum(posdvhistogram[:,ktest+1]))
                
        //xordvrepetitions[1:length(xordvrepetitionstmp),ktest+1]=xordvrepetitionstmp
        //posdvrepetitions[1:length(posdvrepetitionstmp),ktest+1]=posdvrepetitionstmp
    }
                
    printf("Completed task\n.");
    //xorDVtotalrepetitions
    int32_t** xorDVtotalrepetitions = (int32_t**)malloc((maxXORValue+1)*sizeof(int32_t *));
    for(a = 0; a < maxXORValue+1; a++){
        xorDVtotalrepetitions[a] = (int32_t*)malloc(2*sizeof(int32_t));
    }
    //posDVtotalrepetitions
    int32_t** posDVtotalrepetitions = (int32_t**)malloc((maxPOSValue+1)*sizeof(int32_t*));
    for(a = 0; a < maxPOSValue+1; a++){
        posDVtotalrepetitions[a] = (int32_t*)malloc(2*sizeof(int32_t));
    }
    
    if (nRoundsInPattern > 1){
        
        printf("Checking Total Histogram...\n");

                    //xorDVtotalrepetitions[:,1]=collect(0:1:length(xorDVtotalrepetitions[:,1])-1)
                    //xorDVtotalrepetitions[:,2]=counts(TotalDVhistogram[:,2], 0:maximum(TotalDVhistogram[:,2]))
                    //posDVtotalrepetitions[:,1]=collect(0:1:length(posDVtotalrepetitions[:,1])-1)
                    //posDVtotalrepetitions[:,2]=counts(TotalDVhistogram[:,3], 0:maximum(TotalDVhistogram[:,3]))
        printf("Created statistics for combinations.\n");
    }
                    
    durationSteps[4][1] = time(NULL);
    printf("STARTING TO LOOK UP ANOMALOUSLY REPEATED VALUES...\n");
    durationSteps[5][0] = time(NULL);
    //XORextracted_values00, POSextracted_values00 =
    extractAnomalDVSelfConsistency(xordvrepetitions, posDVtotalrepetitions, totalDVHistogram, LN);
    
    durationSteps[5][1] = time(NULL);
    
    /* Thus, selfconsistency is completed. Now, it is time to apply the TRACE rule. */
    printf("APPLYING THE TRACE RULE FOR XOR DV SET...\n");
   // XORextracted_TraceRule=TraceRule(xorDVtotalrepetitions, TotalDVhistogram,XORextracted_values00, LN, RandomnessThreshold)
    //XORextracted_values01=vcat(XORextracted_values00, TotalDVhistogram[XORextracted_TraceRule, 1:2])
    //POSextracted_values01=POSextracted_values00
    
    
    
    /** USING DV VALUES INTRODUCED BY USER. **/
    



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
    
    //22:00
    
    return 0;
}
