#ifndef libs_h
#define libs_h

#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/**** CONSTANTS ****/

/* The following constant is extremely important: it determines the
 * expected number values that separates randomness from causality. We set
 * this value to 0.05 (1 chance out of 20) */
#define randomnessThreshold 0.05 //double

/* A value decided by the researcher and supposed to be higher than the
 * largest expected MCU */
#define UnrealMCUsize 30

/* The second part of the XORing rule (DV XOR CritXOR = CritXOR) has demonstrated
 * to be a bottleneck and it is better to use it only if the highest NDV is low
 * enough. The following constant can be used to determine the size of the DV
 * sets from which the rule cannot be applied and restricted only to know critical
 * know XOR values */
#define NDVforShallowXORrule 1E5

/* Sometimes, the second pass of the MCU rule leads to crush. In this case,
 * swap this constant to false */
#define DoSecondPassMCURule false

#define NMaxAddressesInEvent 200


/*
 * This function allows getting characters from the file, transforming them into fields
 */
char* getfield(FILE* file, char* taken){
    int c = 0, elems = 0;
    char* ch = malloc(sizeof(char)); // En visual da excepción, se arregla con calloc y tam = 11
    c = getc(file);
    while ((c != ',') && (c != '\n')) {
        elems++;
        ch[0] = (char)c;
        strcat(taken, &ch[0]);
        c = getc(file);
    }
    ch = NULL;
    free(ch);
    taken = realloc(taken, elems);
    return taken;
}

/*
 * This function allows separating actual values of addresses in tests from
 * those added to create a matrix. As these elements are equal to 0xFFFFFFFF,
 * are easy to found. First of all, we verify if the element is present. If not,
 * the function returns the whole column. If so, it looks for the first occurrence
 * of 0xFFFFFFFF and returns a vector starting from 1 to the element before the
 * first occurrence of the anomalous value.
 */
uint32_t* extract_addressvector(uint32_t* columnVector, int* elems){
    int i;
    // Comparing last position of the array with 0xFFFFFFFF == 4294967295
    if(columnVector[*elems-1] == 0xFFFFFFFF){
        i = 0;
        uint32_t* newColumnVector;
        newColumnVector = (uint32_t*)malloc(sizeof(uint32_t)*(*elems));
        while(columnVector[i] != 4294967295){
            newColumnVector[i] = columnVector[i];
            i++;
        }
        *elems = i;
        //Reallocating memory
        newColumnVector = (uint32_t*)realloc(newColumnVector, sizeof(uint32_t)*i);
        return newColumnVector;
    } else {
        return columnVector;
    }
}

/*
 * This function returns a vector with variable lenght indicating the positions
 * of the bits where a known word differs from the pattern
 */
int32_t* flipped_bits(int32_t word, int32_t pattern, int wordWidth){
    int k, nFound = 0, count = 0;
    int32_t* result = calloc(wordWidth, sizeof(int32_t));
    for (k = 0; k < wordWidth; k++) {
        result[k] = 1;
        result[k] *= (wordWidth+1);
    }
    int32_t bitflips = word ^ pattern;
    for(k = 0; k < wordWidth-1; k++){
        if ((bitflips/2) == 1) {
            nFound++;
            result[nFound] = k;
        }
        bitflips = bitflips >> 1;
    }
    //return result[1:findlast(result.!=(wordwidth+1))]
    // Recorremos el vector desde la última posición para encontrar el primero desde el final,
    // que cumple result != wordWith+1, realloc de result, para devolver vector con nuevo tamaño
    k = wordWidth;
    while (result[k] == wordWidth+1) {
        count++;
    }
    result = realloc(result, (wordWidth-count)*sizeof(int32_t));
    return result;
}

uint32_t** Transposed_matrix(uint32_t** A, int NumRow, int NumCol){
	uint32_t** trans;
	int a = 0;
	trans = (int32_t **)malloc(NumRow*sizeof(int32_t *));
	for (a = 0; a < NumRow; a++)
		trans[a] = (int32_t *)malloc(NumCol*sizeof(int32_t));

	for (int a = 0; a < NumRow; a++){
		for (int b = 0; b < NumCol; b++){
			trans[a][b] = A[b][a];
		}
	}
	//mostrar transpuesta
	/*for (int a = 0; a < NumRow; a++){
		printf("\n");
		for (int b = 0; b < NumCol; b++){
			printf(" %x", trans[a][b]);
		}
	}*/

	return trans;
}
//Funcion de relleno de fila superior y columna izquierda vector, elemens y matriz para modificarla


/*
 * Function to create the DV sets associated with the XOR and positive operations.
 * The idea is simple: The address vector is replicated to create a square matrix
 * traspose it, and operating. Every element in elements above the main
 * diagonal is an element of the DV.
 */
uint32_t** create_DVmatrix(uint32_t* addresses, int elems, char* op, uint32_t* RWcycles, int nbits4blocks, long int ln){
    int32_t** DVmatrix = (int32_t **)malloc((elems+2)*sizeof(int32_t *));
    int i, j, cont = 0, ntriu = 0;
    for (i = 0; i < elems+2; i++){
        DVmatrix[i] = (int32_t*)malloc((elems+2)*sizeof(int32_t));
    }
    /* Número de elementos que debería haber en la triangular superior */
    for(i = 0; i < elems; i++){
        ntriu += i;
    }
    //Rellenar fila superior y columna izquierda
    for (j = 1; j < elems+1; j++){
        DVmatrix[0][j] = addresses[j-1];
    }
    for (i = 1; i < elems+1; i++){
        DVmatrix[i][0] = addresses[i-1];
    }

    if (strcmp(op, "xor") == 0){
        for (i = 1; i < elems+1; i++){
            for (j = 2+cont; j < elems+1; j++){
                DVmatrix[i][j] = DVmatrix[0][j] ^ DVmatrix[i][0];
            }
            cont++;
        }
    } else if (strcmp(op, "pos") == 0){
        /* Subtracting is a problem with unsigned integers. It is better to convert
         * them into signed format and later return to the original format. */
        for (i = 1; i < elems+1; i++){
            for (j = 2+cont; j < elems+1; j++){
                DVmatrix[i][j] = abs(DVmatrix[0][j] - DVmatrix[i][0]);
            }
            cont++;
        }
    } else {
        printf("\n A problem creating the DV set. Operation not recognized\n");
        EXIT_FAILURE;
    }

    /* Now, we will incorporate the information about the R-W cycles in the same
     * way. We create a matrix with the vector, traspose it and look for coincidence
     * element to element. Later, the matrix with True/False elements is converted
     * into integer and multiplied element to element with the DVmatrix. */
    //DVTMP1 = repmat(RWcycles, 1, N);
    
    //CoincidenceMatrix = round(Int32, triu(DVTMP1.==DVTMP1'))
    //DVmatrix = DVmatrix.*CoincidenceMatrix
                                          
    /* The same for the block division. We multiply the address vector by 2^bits4blocks/LN
     * and round to the floor integer. Traspose and multiply. */
    
    int32_t lnn = pow(2,24)-1;
    uint32_t* nuevo = malloc(sizeof(uint32_t)*elems);
    for (i = 0; i < elems; i++) {
        nuevo[i] = floor(addresses[i]/lnn);
    }
    
    //DVTMP1 = repmat(addressblocks, 1, N);
                                          
    //CoincidenceMatrix = round(Int32, triu(DVTMP1.==DVTMP1'))
                                                                                
    //DVmatrix = DVmatrix.*CoincidenceMatrix
    
    
    /* both return. Let us change the format a bit... */
    if (strcmp(op, "xor") == 0) {
       // DVmatrix = round(UInt32, DVmatrix)
    }
    return (uint32_t**)DVmatrix;
}

/*
 * This function creates the DVvector cotaining the elements of the triu in
 * the DVmatrix. elems is the vector's size.
 * Returns that vector.
 */
uint32_t* create_DVvector(uint32_t** DVmatrix, int elems){
    int ndv = 0.5*elems*(elems-1);
    uint32_t* DVvector = (uint32_t*)calloc(ndv, sizeof(uint32_t));

    /* The elements of the matrix are recovered. Elements are added fixing the
     * column and going down. */
    int index = 0, cont = 0;
    int kcol, krow;
    for (krow = 1; krow < elems+1; krow++){
        for (kcol = 2+cont; kcol < elems+1; kcol++){
            DVvector[index] = DVmatrix[krow][kcol];
            index++;
        }
        cont++;
    }
    return DVvector;
}

/* Copia de matrices */
void copyOfMatrix(uint32_t** matrix, uint32_t*** matrixbackup, int elems, int ktest){
    int i, j;
    for(i = 0; i < elems+2; i++){ // Se copian con los valores de las direcciones?
        for (j = 0; j < elems+2; j++) {
            matrixbackup[i][j][ktest-1] = matrix[i][j];
        }
    }
}

/* Copia de vector a matriz */
void copyOfVector(uint32_t* vector, uint32_t** vectorbackup, int ndv, int ktest){
    int i;
    for(i = 0; i < ndv; i++){
        vectorbackup[i][ktest-1] = vector[i];
    }
}

/*
 * This function returns an array which has on it the count of each element on the histogram
 */
int32_t* countsOfElems(int32_t** histogram, int maxValue, int col, int LN){
    int32_t* repetitionstmp = (int32_t*)calloc(maxValue, sizeof(int32_t));
    int i;
    for (i = 0; i < LN; i++) {
        if (histogram[i][col] > 0 && histogram[i][col] < maxValue) {
            repetitionstmp[histogram[i][col]]++;
        }
    }
    return repetitionstmp;
}

/*
 * This function looks for the number of appearances of an element
 * vector[position]. nvd is the vector's size.
 * Return the number of appearances of that element.
 */
int lookForElem(uint32_t* vector, int ndv, int position){
    int sum = 0, i;
    for (i = 0; i < ndv; i++) {
        if (vector[position] == vector[i]) {
            sum++;
        }
    }
    return sum;
}

/* 
 * The use of this function is quite simple. The DV is checked step by step
 * and every time that an element K appears, histogram(k) increases 1.
 * LN determines the length of the set of possible values.
 */
int32_t* create_histogram(uint32_t* vector, int ndv, long int ln){
    int32_t* histogram = calloc(ln, sizeof(int32_t));
    int k, sum = 0;
    for (k = 0; k < ndv; k++) {
        sum = lookForElem(vector, ndv, k);
        histogram[k] = sum;
    }
    return histogram;
}
int calculaFatorial(int num)
{
    return ((num <= 1) ? 1 : (num * calculaFatorial(num - 1)));
}

long long binomial(int n, int p)
{
    return (calculaFatorial(n) / (calculaFatorial(p)*calculaFatorial(n - p)));
}


/* This function just implements the theoretical expresions to determine expected
 * repetitions.  Variable names are meaningful. */
double ExpectedRepetitions(int m, int ndv, long int LN, char* operation){
  double result = 0,logresulttmp = 0,i;

  if (strcmp(operation, "pos") == 0){
    //Difference of values
    if(fabs(ndv/LN) < 0.1){
        logresulttmp = (float)((1 - m)*log(LN) + m*log(2) - log(m + 1) + log(binomial((int)ndv, m)));
      result = exp(logresulttmp);
      result = result*(1+(3*(m)-2*(int64_t)(ndv))/(LN+1))*(1+2*(ndv-m)/(LN-2)/(m+2));
    }else{
        //          OJOOOOOOOOOOOOOOOOOOOOOOOO  EL COLLECT...
        /*pk=2/LN/(LN+1)*(LN+1-collect(1:1:LN));
      result = sum(binomial(BigInt(ndv), m)*(pk.^m).*(1-pk).^(ndv-m));*/

        for (i = 0; i < LN; i++){
            double pk = 2 / LN / (LN + 1)*(LN + 1 - i);
            result += (binomial((double)(ndv), m)*(pow(pk, m))*pow((1 - pk), (ndv - m)));
        }
    }
  }
  else if (strcmp(operation, "xor") == 0){
    // XOR
      logresulttmp = (float)(log(binomial((double)(ndv), m)) + (ndv - m)*log(LN - 1) - (ndv - 1)*log(LN));
    result = exp(logresulttmp);
  }else{
    //Badly specified function.
    printf("\n\tAn unknown operation to calculate expected repetitions. Returning a strange value\n");
    result = -1;
  }
  return result;
}

int32_t ExcessiveRepetitions(int32_t** RepInHistVector, int RepInHistVectorCol, long int LN, char* operation, int threshold){

    // This function allows calculating the lowest number of repetitions from which
    // the expected number is below threshold.
    // RepInHistVector is a column vector containing the statistics.It is assumed
    // that the first element correspond to 0 repetitions.

    int32_t nthreshold = 0, ndv = 0;
    int32_t nmax = sizeof(RepInHistVector) / sizeof(int32_t);       //si no le pasamos el tamaño

    int32_t* OccurrenceIndex = (int32_t*)calloc(sizeof(int32_t),(nmax - 1));
    for (int i = 0; i < nmax; i++)
        OccurrenceIndex[i] = i;

    for (int i = 0; i < nmax; i++)
        ndv += ((int32_t)RepInHistVector[i][RepInHistVectorCol] * (int32_t)OccurrenceIndex[i]);

    //Inner product and later addition
    //LN = sum(RepInHistVector);
    //just a trick to reobtain the LN value from the data.
    //REMOVED: LN included as an function argument since it is not easily deducible
    // in the case of dividing the memory into blocks.
    int k0 = 0;
    if (strcmp(operation, "pos") == 0){
        if ((ExcessiveRepetitions(RepInHistVector, RepInHistVectorCol, LN, "xor", threshold) - 1) > 0)
            k0 = ExcessiveRepetitions(RepInHistVector, RepInHistVectorCol, LN, "xor", threshold) - 1;
    }
    // This is a trick to speed up "pos" calculations. Instead of starting to
    // iterate from 0, we calculate the threshold for XOR, which is always the
    // lowest and it is really easy to determine.

    //for(int k = k0; k < nmax; k++){   // mejor con un while, porque en cuanto encuentre el valor, se sale
    int k = k0;
    while ((nthreshold != k) && (k < nmax)){
        if (ExpectedRepetitions(k, LN, ndv, operation) < threshold){
            nthreshold = k;
        }
        k++;
    }
    return nthreshold;
}

/*
 * This function lists the element anomalously repeated in the histogram.
 * occurrences is a column vector [V(0)...V(K)] where V(K) indicates the number
 * of elements that are repeated K times.
 * nthreshold indicates the limit value from which repetitions should not
 * be expected in only-SBU experiments.
 */
int32_t** find_anomalies_histogram(uint32_t** histogram, int32_t histogramLenght, int histogramCol, int32_t** occurrences, int32_t occurrencesLenght, int occurrencesCol, int32_t nthreshold, int* n_anomalous_values){ //n_anomalous_values por referencia porque necesitamos el valor
    /* The total number of elements anomalously repeated. We use +1 since the first
     * element of the vector is related to 0 repetitions. 
     */
    int i, index = 0, k, k2;
    int32_t nAnomalies = 0;
    for (i = nthreshold; i < occurrencesLenght; i++) { // nthreshold+1 o nthreshold??????
        nAnomalies += occurrences[i][occurrencesCol];
    }
    /*The number of repetitions of the most repeated element. Coincides roughly with the vector length. */
    int32_t highestAnomaly = occurrencesLenght - 1;
    //initialize the value to return
    int32_t** anomalies = (int32_t**)calloc(nAnomalies, sizeof(int32_t*));
    for (i = 0; i < nAnomalies; i++) {
        anomalies[i] = (int32_t*)calloc(2, sizeof(int32_t));
    }
    /* the following variable will be used to indicate the number of elements different
     * from 0, above nthreshold, in the occurrences vector. */
    *n_anomalous_values = 0;
    for (k = highestAnomaly; k > nthreshold; k--) { // Los dos +1 o no????
        /* I've prefered to implement this solution instead of using the native function
         * "find()" since I believe this solution is faster as it includes breaks once
         * the number of anomalous occurrences are achieved. */
        if (occurrences[k][occurrencesCol] != 0){
            int n_occurrence_value = 0;
            int occurrence_value = k-1;
            *n_anomalous_values += 1;
            for (k2 = 1; k2 < histogramLenght; k2++) {
                if (histogram[k2][histogramCol] == occurrence_value){
                    anomalies[index][0] = k2;
                    anomalies[index][1] = occurrence_value;
                    index++;
                    n_occurrence_value++;
                    if (n_occurrence_value == occurrences[k][occurrencesCol]) {
                        break;
                    }
                }
            }
        }
    }
    return anomalies;
/*Example of return:
 #Anomalies=[
 #49152  13
 #  6   6
 #49158   4
 #1974137   4
 #1974143   4
 #2023289   4
 #2023295   4
 #57344    3
 #893437   3]
 #n_anomalous_values=4
 #Obviously, since {13, 6, 4, 3} are four posibilities.*/
}

/*
 # Now, a simple function to extract a subset with the most repeated values in
 # the anomalous patterns. Take into account that we look at the number of repeated
 # values and not to the number of elements. For instance, if we have the following
 # anomalous_repetitions vector:
 #
 #   0x0000c000  13
 #   0x00000006  6
 #   0x0000c006  4
 #   0x001e1f79  4
 #   0x001e1f7f  4
 #   0x001edf79  4
 #   0x001edf7f  4
 #   0x0000e000  3
 #   0x000da1fd  3
 #
 # And we specify n=1, the first row will be returned. If, n=2, rows 1 & 2
 # but, if n=3, the function returns rows 1->7 as they contain 3 different
 # repetition values. If n = 4 or more, the function returns the input matrix
 # unchanged.
 */
int32_t** extract_some_critical_values(int32_t** anomalous_repetitions, int32_t anomalous_repetitionsLenght, int* n_anomalous_repetitions, int n){
    int observed_values = 1, nrow = 1, k;
    if (*n_anomalous_repetitions <= n) {
        int32_t** extracted_values = calloc(anomalous_repetitionsLenght, sizeof(int32_t*));
        for (k = 0; k < anomalous_repetitionsLenght; k++) {
            extracted_values[k] = calloc(2, sizeof(int32_t));
        }
        //extracted_values = anomalous_repetitions; //Copia de la matriz en otra matriz
        for (k = 0; k < anomalous_repetitionsLenght; k++) {
            extracted_values[k][0] = anomalous_repetitions[k][0];
            extracted_values[k][1] = anomalous_repetitions[k][1];
        }
        return extracted_values;
    } else {
        for (k = 1; k < anomalous_repetitionsLenght; k++){
            if (anomalous_repetitions[k][1] != anomalous_repetitions[k-1][1]) {
                if (observed_values == n) {
                    nrow = k - 1;
                }
                observed_values++;
            }
        }
        int32_t** extracted_values = calloc(nrow, sizeof(int32_t*));
        for (k = 0; k < nrow; k++) {
            extracted_values[k] = calloc(2, sizeof(int32_t));
        }
        //extracted_values = anomalous_repetitions[1:nrow,:];
        for (k = 0; k < nrow; k++) {
            extracted_values[k][0] = anomalous_repetitions[k][0];
            extracted_values[k][1] = anomalous_repetitions[k][1];
        }
        return extracted_values;
    }
}

/*
 #This function marks the elements in the DV matrix that are equal to one of
 # the elements of the critical_values vector. The procedure is quite simple to
 # understand. First of all, a boolean matrix with the size of the DV matrix
 # is created to initialize the result to return.
 # Later, a matrix with DV size and with all the elements equal to one of
 # the critical values is created
 #
 #           critical_values[kvalue,1]*ones(UInt32, size(davmatrix))
 #
 # and compared element to element with the DV matrix:
 #
 # round(Bool, davmatrix.==critical_values[kvalue,1]*ones(UInt32, size(davmatrix)))
 #
 # and ORED with the preliminary result value. This recursive solution is
 # implemented in a loop and tried to be parallelized with the instruction @simd
 # If this works, only God knows.
 */
bool** marking_addresses(int32_t** davmatrix, int32_t davmatrixRows, int32_t davmatrixCols, int32_t** critical_values, int32_t critical_valuesLenght){
    int i, x, y;
    int32_t value;
    bool** result = calloc(davmatrixRows, sizeof(bool*));
    for (i = 0; i < davmatrixRows; i++) {
        result[i] = calloc(davmatrixCols, sizeof(bool));
    }
    /*for value in critical_values[:,1]
        result = result | (round(Bool, davmatrix.==value*ones(UInt32, size(davmatrix))))
    end*/
    for (i = 0; i < critical_valuesLenght; i++) {
        value = critical_values[i][0];
        //result = result | (round(Bool, davmatrix.==value*ones(UInt32, size(davmatrix))))
        for (x = 0; x < davmatrixRows; x++) { // Coste muy alto porque hay que recorrer la matriz entera para cada elemento
            for (y = 0; y < davmatrixCols; y++) {
                if (davmatrix[x][y] == value) {
                    result[x][y] = true;
                }
            }
        }
    }
    return result;
}

/*
 * Esta función, crea un vector temporal de tamaño thesummaryCols, para hacer la unión entre las filas addressRowMin y 
 * addressRowMax, colocando los elementos en orden mientras se insertan.
 */
int32_t* unionAgrupate_mcus(int32_t** thesummary, int32_t thesummaryCols, int32_t addressRowMin, int32_t addressRowMax){
    int i, j, x;
    int32_t* thesummarytemp = calloc(thesummaryCols, sizeof(int32_t)); // Necesitamos un vector para
    // guardar los elementos de las dos filas, para poder hacer la unión
    j = 0; // Indice de thesummarytemp
    for (i = 0; i < thesummaryCols; i++) {
        if (thesummary[addressRowMin][i] != 0) {
            if (thesummarytemp[0] == 0) { // Esta vacío
                thesummarytemp[j] = thesummary[addressRowMin][i];
                j++;
            } else if(thesummarytemp[j-1] < thesummary[addressRowMin][i]){ // Si el valor último menor
                thesummarytemp[j] = thesummary[addressRowMin][i];
                j++;
            } else if(thesummarytemp[j-1] > thesummary[addressRowMin][i]){ // Si el valor último mayor, se abre hueco
                x = j - 1;
                while (thesummarytemp[x] > thesummary[addressRowMin][i]) {
                    thesummarytemp[x+1] = thesummarytemp[x]; // Se desplaza uno
                    x--;
                }
                thesummarytemp[x+1] = thesummary[addressRowMin][i]; // Se copia en su hueco correspondiente
            }
        }
    }for (i = 0; i < thesummaryCols; i++) {
        if (thesummary[addressRowMax][i] != 0) {
            if (thesummarytemp[0] == 0) { // Esta vacío
                thesummarytemp[j] = thesummary[addressRowMax][i];
                j++;
            } else if(thesummarytemp[j-1] < thesummary[addressRowMax][i]){ // Si el valor último menor
                thesummarytemp[j] = thesummary[addressRowMax][i];
                j++;
            } else if(thesummarytemp[j-1] > thesummary[addressRowMax][i]){ // Si el valor último mayor, se abre hueco
                x = j - 1;
                while (thesummarytemp[x] > thesummary[addressRowMax][i]) {
                    thesummarytemp[x+1] = thesummarytemp[x]; // Se desplaza uno
                    x--;
                }
                thesummarytemp[x+1] = thesummary[addressRowMax][i]; // Se copia en su hueco correspondiente
            }
        }
    }
    return thesummarytemp;
}

/*
 # keep in mind that davmatrix is a boolean matrix containing FALSE if two elements
 # are not involved in the same MCU and TRUE if so.
 
 # First of all, let us detect the elements different from FALSE. It is easy
 # with the function FIND. The problem is that Julia considers matrix as vectors
 # following the column order. Therefore, it is necessary to make some calculations
 # to place the TRUE values in the matrix. Let us remember that, in Julia,
 # the element [row, col] in a matrix is the element (col-1)*NROws+row in the
 # equivalent vector, with NROws the total number of columns.
 */

int32_t** agrupate_mcus(bool** davmatrix, int32_t davmatrixRows, int32_t davmatrixCols, int* nTotalEvents, int* largestMCU){
    //Let us locate the elements
    //Con FIND tenemos que encontrar los indices cuyos valores no sean cero
    int count = 0, newCount = 0, i, j, k, x, xx;
    int32_t* nonzerovectorelementstmp = calloc(davmatrixRows*davmatrixCols, sizeof(int32_t));
    for (i = 0; i < davmatrixRows; i++) {
        for (j = 0; j < davmatrixCols; j++) {
            nonzerovectorelementstmp[i * davmatrixCols + j] = -1; //metemos -1 de entrada para ahorrar una vuelta
            if (davmatrix [i][j] != false) {
                count++; //Sabemos cuantos elementos, pero alguno puede estar repetido
                nonzerovectorelementstmp[i * davmatrixCols + j] = i * davmatrixCols + j;
            }
        }
    }
    // Tenemos que volver a recorrer para encontrar los -1 que son vacíos
    int32_t* nonzerovectorelements = calloc(count, sizeof(int32_t)); // Luego redimensionamos
    for (i = 0; i < count; i++) {
        if (nonzerovectorelementstmp[i] != -1) {
            nonzerovectorelements[newCount] = nonzerovectorelementstmp[i];
            newCount++; // número de elementos reales
        }
    }
    nonzerovectorelementstmp = NULL;
    free(nonzerovectorelementstmp);
    nonzerovectorelements = realloc(nonzerovectorelements, newCount);
    //Now, transform the elements in pairs creating a matrix with 2 columns.
    int32_t** relatedpairs = calloc(newCount, sizeof(int32_t*));
    for (i = 0; i < newCount; i++) {
        relatedpairs[i] = calloc(2, sizeof(int32_t));
    }
    
/*    @simd for k1 = 1:length(nonzerovectorelements)
# The first operation calculates the column
        @inbounds relatedpairs[k1,:] = [rem(nonzerovectorelements[k1], NRows) div(nonzerovectorelements[k1],NRows)+1]
# REM means "remainder of division", DIV integer division.
        end*/
    for (i = 0; i < newCount; i++) {
        relatedpairs[i][0] = nonzerovectorelements[i] / davmatrixRows;
        relatedpairs[i][1] = (nonzerovectorelements[i] / davmatrixRows) + 1;
    }
    /* Good, we have created a matrix containing the related pairs. Now, we must
     # group them in larger events.
     
     # First step: A matrix initilizing the events. Perhaps too large, but it is better
     # to have spare space.
     */
    //thesummary = zeros(Int32, length(relatedpairs[:,1]), minimum([length(relatedpairs[:,1]),nmaxaddressesinevent]))
    int32_t** thesummary = calloc(newCount, sizeof(int32_t*));
    int32_t thesummaryRows = newCount, thesummaryCols;
    if (newCount < NMaxAddressesInEvent) { // Mínimo es newCount
        thesummaryCols = newCount;
        for (i = 0; i < newCount; i++) {
            thesummary[i] = calloc(newCount, sizeof(int32_t));
        }
    } else { // Mínimo es NMaxAddressesInEvent
        thesummaryCols = NMaxAddressesInEvent;
        for (i = 0; i < newCount; i++) {
            thesummary[i] = calloc(NMaxAddressesInEvent, sizeof(int32_t));
        }
    }
    //Now, we save the first pair in the first column:
    *nTotalEvents = 0;
    int32_t firstAddress, secondAddress, firstAddressRow = 0, secondAddressRow = 0;
    for (k = 0; k < newCount; k++) {
        //There are several cases. Let us list them in order.
        //Positions are extrated from the pairs.
        firstAddress = relatedpairs[k][0];
        secondAddress = relatedpairs[k][1];
        //TENEMOS QUE MIRAR SI firstAddress y secondAddress están en thesummary
        bool firstAddressFound = false, secondAddressFound = false;
        for (i = 0; i < thesummaryRows; i++) {
            for (j = 0; j < thesummaryCols; j++) {
                if (thesummary[i][j] == firstAddress) {
                    firstAddressFound = true;
                    firstAddressRow = i;
                }
                if (thesummary[i][j] == secondAddress) {
                    secondAddressFound = true;
                    secondAddressRow = i;
                }
            }
        }
        /*
         # Case 1:
         # Both addresses are not present in thesummary. They are placed at a new_row
         # at NTotalEvents
         */
        if (!firstAddressFound && !secondAddressFound){
            thesummary[*nTotalEvents][0] = relatedpairs[k][0];
            thesummary[*nTotalEvents][1] = relatedpairs[k][1];
            *nTotalEvents ++;
        }
        /*
         # Case 2:
         # FirstAddress is present, but not SecondAddress
         */
        if (firstAddressFound && !secondAddressFound){
            //First of all, let us locate the row where the FirstAddress is.
            //Now, it is necessary to locate the first free column to put the new address.
            x = 0;
            while ((thesummary[firstAddressRow][x] != 0) && (x < thesummaryCols)) {
                x++;
            }
            //Fixed problem. We do know exactly the position to put the new value.
            if (x < thesummaryCols) {
                thesummary[firstAddressRow][x] = secondAddress;
            }
        }
        /*
         # Case 3:
         # FirstAddress absent, Second Address present. similar to previous one.
         */
        if (!firstAddressFound && secondAddressFound) {
            //First of all, let us locate the row where the Second Address is.
            //Now, it is necessary to locate the first free column to put the new address.
            x = 0;
            while ((thesummary[secondAddressRow][x] != 0) && (x < thesummaryCols)) {
                x++;
            }
            //Fixed problem. We do know exactly the position to put the new value.
            if (x < thesummaryCols) {
                thesummary[secondAddressRow][x] = firstAddress;
            }
        }
        /*
         #Case 4:
         # Both addresses had been previously detected. Two subcases appear: They
         # are included in the same row (MCU) so the program must skip this pair or
         # are included in differen rows. Therefore, both rows must be carefully merged.
         */
        if (firstAddressFound && secondAddressFound) {
            /*
             # If both values are identical, the pair must be skipped and the program
             # continue. If different, the rows must be merged.
             */
            if (firstAddressRow != secondAddressRow){
                //First step: It is necessary to range the numbers in increasing order.
                int32_t addressRowMin = 0, addressRowMax = 0;
                if (firstAddressRow < secondAddressRow) {
                    addressRowMin = firstAddressRow;
                    addressRowMax = secondAddressRow;
                } else {
                    addressRowMin = secondAddressRow;
                    addressRowMax = firstAddressRow;
                }
                //Now, let us locate elements in RowMin different from 0
                int32_t* thesummarytemp = unionAgrupate_mcus(thesummary, thesummaryCols, addressRowMin, addressRowMax);
                //RawCombinedAddresses =union(thesummary[AddressRowMin, 1:findfirst(thesummary[AddressRowMin,:].==0)-1], thesummary[AddressRowMax, 1:findfirst(thesummary[AddressRowMax,:].==0)-1]);
                /*
                 * But this is a column vector where the elements are not repeated but
                 * not disposed in increasing order. Let us solve this:
                 */
                //RawCombinedAddresses = sort(vec(RawCombinedAddresses))'
                /*
                 * Vec vectorizes the elements of the matrix.
                 * Time to put the elements in the row.
                 */
                //thesummarytemp[AddressRowMin,1:length(RawCombinedAddresses)]=RawCombinedAddresses;
                for (i = 0; i < thesummaryCols; i++) {
                    thesummary[addressRowMin][i] = thesummarytemp[i]; // Se copia el vector de la unión en la matriz
                }
                //Also, as two events have been merged, the number of total events is lower:
                *nTotalEvents--;
            }
        }
    }
    /*
     #Really good. Now, we will separate the cells with information from those with
     #zeros.
     # The number of rows is easy to calculate: NTotalEvents. Concerning the other
     # element:
     */
    int sum = 0;
    for (i = 3; i < thesummaryRows; i++) {
        for (j = 0; j < thesummaryCols; j++) {
            sum += thesummary[i][j];
        }
        if (sum == 0) {
            *largestMCU = i - 1;
            break;
        }
    }
    int32_t** result = calloc(*nTotalEvents, sizeof(int32_t*));
    for (i = 0; i < *nTotalEvents; i++) {
        result[i] = calloc(*largestMCU, sizeof(int32_t));
    }
    for (x = 0; x < *nTotalEvents; x++) {
        for (xx = 0; xx < *largestMCU; xx++) {
            result[x][xx] = thesummary[x][xx];
        }
    }
    return result;
}

/*
 *This is a simple function to cut files and rows of redundant zeros in 1 or 2
 *dimensions matrix.
 */
/*int32_t** CutZerosFromArray(int32_t*** matrix, int32_t matrixRows, int32_t matrixCols, int32_t matrixDim){
    if (matrixCols == 2) {
        // Upss, this is a vector. Be careful.
        7
    } else if (matrixCols > 2){
        
    } else {
        
    }
}*/
/*
function CutZerosFromArray(A)

if ((Ndimension==2)&(1 in Adimensions))
# Upss, this is a vector. Be careful.
if (Adimensions[1]==1)
result =A[findfirst(A.!=0):findlast(A.!=0)]'
else
result =A[findfirst(A.!=0):findlast(A.!=0)]
end

return result


elseif ((Ndimension==2)&!(1 in Adimensions))
# this is a classical matrix.
firstrow=1;
lastrow =Adimensions[1];
firstcol=1;
lastcol=Adimensions[2];
for krow = 1:Adimensions[1]
if (length(find(A[krow,:].!=0))!=0)
break;
end
firstrow +=1;
end
for krow = Adimensions[1]:-1:1
if (length(find(A[krow,:].!=0))!=0)
break;
end
lastrow -=1;
end
for kcol = 1:Adimensions[2]
if (length(find(A[:,kcol].!=0))!=0)
break;
end
firstcol +=1;
end
for kcol = Adimensions[2]:-1:1
if (length(find(A[:,kcol].!=0))!=0)
break;
end
lastcol -=1;
end

return result = A[firstrow:lastrow, firstcol:lastcol]

else

println("\tMatriz dimension different from 1 or 2. Exiting.")
return A
end

end
*/

/*
 # An alternative implementation of the trace rule.
 
 ### Anomal XOR values is used to avoid the repetition of elements. It is just
 ### XORExtractedValues[1,:]
 */
void traceRule(int32_t** xorDVtotalrepetitions, uint32_t** totalDVhistogram, int32_t** AnomalXORvalues, int32_t LN){
    	// Anomal XOR values is used to avoid the repetition of elements. It is just
	// XORExtractedValues[1,:]
	int32_t nBits = (log(LN + 1) / log(2));

	// Elements with 1-trace
	printf("\n\tInvestigating Trace = 1: ");
	//CandidatesT1 = 2.^collect(0:Nbits-1)
	int32_t CandidatesT1 = 0;
	for (int i = 0; i < nBits - 1; i++){
		CandidatesT1 += pow(2, i);
	}

	int32_t NT1Threshold = ExcessiveRepetitions(xorDVtotalrepetitions, 1, LN, "xor", randomnessThreshold);

	//SelectedT1 =
}

// Esta función depende de la operación
int32_t** extractAnomalDVSelfConsistency(char* op, int32_t** opDVtotalrepetitions, int32_t opDVtotalrepetitionsRows, uint32_t** totalDVhistogram, int32_t histogramLenght, long int LN, int nRoundsInPattern, int32_t** opdvhistogram, int32_t*** opdvmatrixbackup, int opdvmatrixbackupRows){
    int32_t opNthreshold;
    int* n_anomalous_values = malloc(sizeof(int));
    int32_t** testOPa;
    int32_t** oPExtracted_values;
    bool oPSelfConsistence = true;
    int n_anomalous_repetitions = 1, i, j, test;
    printf("\n\tDetermining the threshold for repetition excess:\n");
    if (strcmp(op, "xor") == 0) {
        printf("\n\t\tXOR operation...");
        opNthreshold = ExcessiveRepetitions(opDVtotalrepetitions, 1, LN, "xor", randomnessThreshold);
        testOPa = find_anomalies_histogram(totalDVhistogram, histogramLenght, 1, opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, opNthreshold, n_anomalous_values);
        
        printf("\n\tXOR operation:\n");
        oPExtracted_values = extract_some_critical_values(testOPa, (sizeof(testOPa)/sizeof(int32_t))/2, n_anomalous_values, 1);
        
        while (oPSelfConsistence) {
            printf("\t\tStep %d ", n_anomalous_repetitions);
            oPExtracted_values = extract_some_critical_values(testOPa, (sizeof(testOPa)/sizeof(int32_t))/2, n_anomalous_values, n_anomalous_repetitions);
            int nOPextrValues = (sizeof(testOPa)/sizeof(int32_t))/2;
            int32_t** opPartRepsSelfCons = calloc(nOPextrValues, sizeof(int32_t*));
            for (i = 0; i < nOPextrValues; i++) {
                opPartRepsSelfCons[i] = calloc(nRoundsInPattern, sizeof(int32_t));
            }
            for (i = 0; i < nOPextrValues; i++) {
                for (j = 0; j < nRoundsInPattern; j++) {
                    opPartRepsSelfCons[i][j] = opdvhistogram[oPExtracted_values[i][0]][j+1]; //lolllz
                }
            }
            for (test = 0; test < nRoundsInPattern; test++) {
                printf("Test %d", test);
                int32_t** opdvmatrix = calloc(opdvmatrixbackupRows, sizeof(int32_t*));
                for (i = 0; i < opdvmatrixbackupRows; i++) {
                    opdvmatrix[i] = calloc(opdvmatrixbackupRows, sizeof(int32_t));
                }
                for (i = 0; i < opdvmatrixbackupRows; i++) {
                    for (j = 0; j < opdvmatrixbackupRows; j++) {
                        opdvmatrix[i][j] = opdvmatrixbackup[i][j][test];
                    }
                }
                bool** opMarkedPairs = marking_addresses(opdvmatrix, opdvmatrixbackupRows, opdvmatrixbackupRows, oPExtracted_values, nOPextrValues);
                int32_t** opProposedMcus = agrupate_mcus(opMarkedPairs, opdvmatrixbackupRows, opdvmatrixbackupRows);
                // LargestMCUSize = length(XORproposed_MCUs[1,:])
                int continuation;
                // Continuation=(length(find(XORPartRepsSelfCons[:,ktest].<=LargestMCUSize))==0)
                if (continuation == 0) {
                    oPSelfConsistence = false;
                }
            }
            if (oPSelfConsistence) {
                n_anomalous_repetitions++;
                printf("\n");
                if (n_anomalous_repetitions > *n_anomalous_values) {
                    oPSelfConsistence = false;
                    printf("\n\t\tNo more anomalous elements to check. Exiting\n");
                } else {
                    printf("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n");
                    n_anomalous_repetitions--;
                    oPExtracted_values = extract_some_critical_values(testOPa, (sizeof(testOPa)/sizeof(int32_t))/2, *n_anomalous_values, n_anomalous_repetitions);
                }
            }
        }
        
    } else if(strcmp(op, "pos") == 0){
        printf("\n\t\tPOS operation (It can take a long if there are too many addresses)...\n");
        opNthreshold = ExcessiveRepetitions(opDVtotalrepetitions, 1, LN, "pos", randomnessThreshold);
        testOPa = find_anomalies_histogram(totalDVhistogram, histogramLenght, 2, opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, opNthreshold, n_anomalous_values);
        
        printf("\n\tPOS operation:\n");
        oPExtracted_values = extract_some_critical_values(testOPa, (int32_t)(sizeof(testOPa)/sizeof(int32_t)), *n_anomalous_values, 1);
        while (oPSelfConsistence) {
            printf("\t\tStep %d ", n_anomalous_repetitions);
            oPExtracted_values = extract_some_critical_values(testOPa, (sizeof(testOPa)/sizeof(int32_t)), *n_anomalous_values, n_anomalous_repetitions);
            int nOPextrValues = (sizeof(testOPa)/sizeof(int32_t))/2;
            int32_t** opPartRepsSelfCons = calloc(nOPextrValues, sizeof(int32_t*));
            for (i = 0; i < nOPextrValues; i++) {
                opPartRepsSelfCons[i] = calloc(nRoundsInPattern, sizeof(int32_t));
            }
            for (i = 0; i < nOPextrValues; i++) {
                for (j = 0; j < nRoundsInPattern; j++) {
                    opPartRepsSelfCons[i][j] = opdvhistogram[oPExtracted_values[i][0]][j+1]; //lolllz
                }
            }
            for (test = 0; test < nRoundsInPattern; test++) {
                printf("Test %d", test);
                int32_t** opdvmatrix = calloc(opdvmatrixbackupRows, sizeof(int32_t*));
                for (i = 0; i < opdvmatrixbackupRows; i++) {
                    opdvmatrix[i] = calloc(opdvmatrixbackupRows, sizeof(int32_t));
                }
                for (i = 0; i < opdvmatrixbackupRows; i++) {
                    for (j = 0; j < opdvmatrixbackupRows; j++) {
                        opdvmatrix[i][j] = opdvmatrixbackup[i][j][test];
                    }
                }
                bool** opMarkedPairs = marking_addresses(opdvmatrix, opdvmatrixbackupRows, opdvmatrixbackupRows, oPExtracted_values, nOPextrValues);
                int32_t** opProposedMcus = agrupate_mcus(opMarkedPairs, opdvmatrixbackupRows, opdvmatrixbackupRows);
                // LargestMCUSize = length(XORproposed_MCUs[1,:])
                int continuation;
                // Continuation=(length(find(XORPartRepsSelfCons[:,ktest].<=LargestMCUSize))==0)
                if (continuation == 0) {
                    oPSelfConsistence = false;
                }
            }
            if (oPSelfConsistence) {
                n_anomalous_repetitions++;
                printf("\n");
                if (n_anomalous_repetitions > *n_anomalous_values) {
                    oPSelfConsistence = false;
                    printf("\n\t\tNo more anomalous elements to check. Exiting\n");
                } else {
                    printf("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n");
                    n_anomalous_repetitions--;
                    oPExtracted_values = extract_some_critical_values(testOPa, (sizeof(testOPa)/sizeof(int32_t))/2, *n_anomalous_values, n_anomalous_repetitions);
                }
            }
        }
        
    }
    return oPExtracted_values;
}




/*




LargestMCUSize = length(XORproposed_MCUs[1,:])
Continuation=(length(find(XORPartRepsSelfCons[:,ktest].<=LargestMCUSize))==0)
if (!Continuation)
XORSelfConsistence = false
end

end

if (XORSelfConsistence)
n_anomalous_repetitions +=1;
print("\n")
if (n_anomalous_repetitions>testXORb)
XORSelfConsistence = false
print("\n\t\tNo more anomalous elements to check. Exiting\n")
end
else
print("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n")
n_anomalous_repetitions -=1
XORextracted_values=extract_some_critical_values(testXORa, testXORb, n_anomalous_repetitions)
end

end


POSextracted_values=extract_some_critical_values(testPOSa, testPOSb, 1)

while (POSSelfConsistence)
print("\t\tStep ", n_anomalous_repetitions, ": ")

POSextracted_values=extract_some_critical_values(testPOSa, testPOSb, n_anomalous_repetitions)
NPOSextrValues = length(POSextracted_values[:,1])
POSPartRepsSelfCons = zeros(Int32, NPOSextrValues, NRoundsInPattern)

for kval=1:NPOSextrValues
POSPartRepsSelfCons[kval,:]=posdvhistogram[POSextracted_values[kval,1], 2:end]
end

#n_anomalous_repetitions +=1

for ktest = 1:NRoundsInPattern
print("Test ", ktest, " ")

posdvmatrix = posdvmatrixbackup[1:NAddressesInRound[ktest],1:NAddressesInRound[ktest],ktest]
POSmarked_pairs = marking_addresses(posdvmatrix, POSextracted_values)
POSproposed_MCUs = agrupate_mcus(POSmarked_pairs)

LargestMCUSize = length(POSproposed_MCUs[1,:])

Continuation=(length(find(POSPartRepsSelfCons[:,ktest].<=LargestMCUSize))==0)
if (!Continuation)
POSSelfConsistence = false
end

end

if (POSSelfConsistence)
print("\n")
n_anomalous_repetitions +=1;
if (n_anomalous_repetitions>testXORb)
POSSelfConsistence = false
print("\n\t\tNo more anomalous elements to check. Exiting.\n")
end
else
print("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n")
n_anomalous_repetitions -=1
POSextracted_values=extract_some_critical_values(testPOSa, testPOSb, n_anomalous_repetitions)
end
end

return XORextracted_values, POSextracted_values

end
*/

int32_t*** index2address(int32_t*** indexMatrix, int indexMatrixRows, int indexMatrixCols, int indexMatrixDims, int32_t** addressMatrix){
    int i, j, z;
    int32_t*** result = calloc(indexMatrixRows, sizeof(int32_t**));
    for (i = 0; i < indexMatrixRows; i++) {
        result[i] = calloc(indexMatrixCols, sizeof(int32_t*));
        for (j = 0; j < indexMatrixDims; j++) {
            result[i][j] = calloc(indexMatrixDims, sizeof(int32_t));
        }
    }
    for (z = 0; z < indexMatrixDims; z++) {
        for (i = 0; i < indexMatrixRows; i++) {
            for (j = 0; j < indexMatrixCols; j++) {
                if (indexMatrix[i][j][z] == 0) {
                    break;
                } else {
                    result[i][j][z] = addressMatrix[indexMatrix[i][j][z]][z];
                }
            }
        }
    }
    return result;
}

/*
 # This function investigates occurrences of several bitflips in an only one
 # address. It returns a matrix of variable number of rows but with four collumns.
 # The information is:
 # First Column: Experiment in which it was observed.
 # Second column: The position of the address in the address vector.
 # Third Column: The exact address.
 # Fourth Column: The number of flipped bits.
 */
void locate_mbus(int32_t** content, int32_t contentRows, int32_t contentCols, int32_t** pattern, int datawidth){
    int32_t nRounds = contentCols/3;
    int i, j, nDetected = 0;
    uint32_t** summary = calloc(3*contentRows, sizeof(uint32_t*));
    for (i = 0; i < 3*contentRows; i++) {
        summary[i] = calloc(4, sizeof(uint32_t));
    }
    for (i = 0; i < nRounds; i++) {
        for (j = 0; j < contentRows; j++) {
            if (content[j][3*i-2] == 0xFFFFFFFF) {
                break;
            } else {
                int32_t* corruptedBits;
                corruptedBits = flipped_bits(content[contentRows][3*i-1], pattern[i][1], datawidth);
                //if (length(corruptedBits) > 1){ //Lenght no vale
                    // Summary[NDetected,:]=[ktest kaddress Content[kaddress, (3ktest-2)] length(CorruptedBits)]
                    nDetected++;
                //}
            }
        }
    }
    //return cutZerosFromArray(summary);
}

int32_t** condensate_summary(int32_t*** summary3D, int summary3DRows, int summary3DCols, int summary3DDims, int32_t** content, int contentRows){
    int i, k;
    int32_t** condensateSummary = calloc(summary3DCols, sizeof(int32_t*));
    for (i = 0; i < summary3DDims+1; i++) {
        condensateSummary[i] = calloc(summary3DDims+1, sizeof(int32_t));
        // CondensateSummary[:,1]=round(Int32,linspace(1,LargestMCU, LargestMCU))
        // condensateSummary[i][0] =
    }
    for (i = 1; i < summary3DDims+1; i++) {
        // elems tendrá la longitud del nuevo vector
        int elems = contentRows;
        uint32_t* columnvector = (uint32_t*)malloc(sizeof(uint32_t)*(elems));
        for (k = 0; k < elems; ++k) {
            columnvector[k] = content[k][3*i-3];
        }
        uint32_t* addresses = extract_addressvector(columnvector, &elems);
        //Ahora para cada test: Summary = Summary3D[:,:,ktest]
        // for kmcu = LargestMCU:-1:2
        for (k = summary3DCols; k > 1; k--) {
            //NEvents = length(find(Summary[:,kmcu]))
            //CondensateSummary[kmcu, ktest+1] = NEvents - sum(CondensateSummary[kmcu+1:end, ktest+1])
        }
        // Apparently, the following solution is a bit faster.
        //CondensateSummary[1, ktest+1] = NAddressesInRound - length(find(Summary))
    }
    return condensateSummary;
}


int32_t** propose_MCUs(char* op, uint32_t*** XORDVmatrix3D, int numRows, int numCols, int dim, int32_t*** POSDVmatrix3D, int32_t** anomalXOR, int32_t anomalXORlength, int32_t** anomalPOS, int32_t anomalPOSlength, int* nAddressesInRound){

	int i, j, z, ktest;
	numRows = ceil(numRows / 2 - 1);
	if (strcmp(op, "xor") == 0){

		int32_t*** XOR_summary = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			XOR_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				XOR_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		int32_t*** xordvmatrix = calloc(nAddressesInRound[ktest], sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			xordvmatrix[i] = calloc(nAddressesInRound[ktest], sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				xordvmatrix[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < dim; i++){

			for (j = 0; j < numRows; i++) {
				for (z = 0; z < numRows; i++) {
					xordvmatrix[j][z] = XORDVmatrix3D[j][z][ktest];
				}
			}


			bool** XORmarked_pairs = marking_addresses(xordvmatrix, numRows, numCols, anomalXOR, anomalXORlength);
			//The following results only concern an operation.Just present for illustrative
			//operations.In practical use, only CMB should be used.
			int32_t** XORproposed_MCUs = agrupate_mcus(XORmarked_pairs, numRows, numCols);


			//Let us save everything.
			for (i = 0; i < numRows; i++) {
				for (j = 0; j < numCols; i++) {
					XOR_summary[i][j][ktest] = XORproposed_MCUs[i][j];
				}
			}
			// However, as there are too many zeros, it is interesting to reshape the matrix.
			// Thus, some memory is saved.But we must be careful with the different dimensions
			// of every partial XY matrix in XYZ arrays.A simple but unelegant way is :
			int32_t*** finalXOR_summary = calloc(numRows, sizeof(int32_t**));
			for (i = 0; i < numRows; i++) {
				finalXOR_summary[i] = calloc(numCols, sizeof(int32_t*));
				for (j = 0; j < dim; j++) {
					finalXOR_summary[i][j] = calloc(dim, sizeof(int32_t));
				}
			}
			//return finalXOR_summary = CutZerosFromMCUsummary(XOR_summary, numRows, numCols, dim);	//Falta terminar hasta que sepamos que valores pasamos
		}
	}
	else if (strcmp(op, "pos") == 0){

		int32_t*** POS_summary = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			POS_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				POS_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		int32_t*** posdvmatrix = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			posdvmatrix[i] = calloc(numCols, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				posdvmatrix[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < dim; i++){

			for (j = 0; j < numRows; i++) {
				for (z = 0; z < numRows; i++) {
					posdvmatrix[j][z] = POSDVmatrix3D[j][z][ktest];
				}
			}


			bool** POSmarked_pairs = marking_addresses(posdvmatrix, numRows, numCols, anomalPOS, anomalPOSlength);
			//The following results only concern an operation.Just present for illustrative
			//operations.In practical use, only CMB should be used.
			int32_t** POSproposed_MCUs = agrupate_mcus(POSmarked_pairs, numRows, numCols);


			//Let us save everything.
			for (i = 0; i < numRows; i++) {
				for (j = 0; j < numCols; i++) {
					POS_summary[i][j][ktest] = POSproposed_MCUs[i][j];
				}
			}
			// However, as there are too many zeros, it is interesting to reshape the matrix.
			// Thus, some memory is saved.But we must be careful with the different dimensions
			// of every partial XY matrix in XYZ arrays.A simple but unelegant way is :
			int32_t*** finalPOS_summary = calloc(numRows, sizeof(int32_t**));
			for (i = 0; i < numRows; i++) {
				finalPOS_summary[i] = calloc(numCols, sizeof(int32_t*));
				for (j = 0; j < dim; j++) {
					finalPOS_summary[i][j] = calloc(dim, sizeof(int32_t));
				}
			}
			//return finalPOS_summary = CutZerosFromMCUsummary(POS_summary, numRows, numCols, dim);	//Falta terminar hasta que sepamos que valores pasamos
		}
	}
	else{//op == CMB

		int32_t*** CMB_summary = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			CMB_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				CMB_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}
		int32_t*** xordvmatrix = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			xordvmatrix[i] = calloc(numCols, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				xordvmatrix[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		int32_t*** posdvmatrix = calloc(numRows, sizeof(int32_t**));
		for (i = 0; i < numRows; i++) {
			posdvmatrix[i] = calloc(numCols, sizeof(int32_t*));
			for (j = 0; j < dim; j++) {
				posdvmatrix[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < dim; i++){

			for (j = 0; j < numRows; i++) {
				for (z = 0; z < numRows; i++) {
					xordvmatrix[j][z] = XORDVmatrix3D[j][z][ktest];
				}
			}
			for (j = 0; j < numRows; i++) {
				for (z = 0; z < numRows; i++) {
					posdvmatrix[j][z] = POSDVmatrix3D[j][z][ktest];
				}
			}

			bool** XORmarked_pairs = marking_addresses(xordvmatrix, numRows, numCols, anomalXOR, anomalXORlength);
			bool** POSmarked_pairs = marking_addresses(posdvmatrix, numRows, numCols, anomalPOS, anomalPOSlength);

			bool** XOR_POS_marked_pairs = calloc(numRows, sizeof(bool*));
			for (i = 0; i < numRows; i++) {
				XOR_POS_marked_pairs[i] = calloc(numCols, sizeof(bool));
			}

			for (j = 0; j < numRows; i++) {
				for (z = 0; z < numRows; i++) {
					XOR_POS_marked_pairs[j][z] = XORmarked_pairs[j][z] | POSmarked_pairs[j][z];
				}
			}

			int32_t** CMBproposed_MCUs = agrupate_mcus(XOR_POS_marked_pairs, numRows, numCols);

			//Let us save everything.
			for (i = 0; i < numRows; i++) {
				for (j = 0; j < numCols; i++) {
					CMB_summary[i][j][ktest] = CMBproposed_MCUs[i][j];
				}
			}
			// However, as there are too many zeros, it is interesting to reshape the matrix.
			// Thus, some memory is saved.But we must be careful with the different dimensions
			// of every partial XY matrix in XYZ arrays.A simple but unelegant way is :
			int32_t*** finalCMB_summary = calloc(numRows, sizeof(int32_t**));
			for (i = 0; i < numRows; i++) {
				finalCMB_summary[i] = calloc(numCols, sizeof(int32_t*));
				for (j = 0; j < dim; j++) {
					finalCMB_summary[i][j] = calloc(dim, sizeof(int32_t));
				}
			}
			//return finalCMB_summary = CutZerosFromMCUsummary(CMB_summary, numRows, numCols, dim);	//Falta terminar hasta que sepamos que valores pasamos
		}
	}
}

int32_t** CriticalXORValuesFromClusters(int32_t*** PrelMCUSummary,int numRows, int numCols, int nRounds, int32_t** AddressMatrix, int32_t** XORextracted_values, int32_t**xorDVtotalrepetitions, int32_t** DVHistogram, long int LN){


	int32_t* NewCandidates;
	int ktest,kRow,kAddress1,kAddress2,i,j;

	int32_t*** index1 = calloc(numRows, sizeof(int32_t**));
	for (i = 0; i < numRows; i++) {
		index1[i] = calloc(numCols, sizeof(int32_t*));
		for (j = 0; j < nRounds; j++) {
			index1[i][j] = calloc(nRounds, sizeof(int32_t));
		}
	}


	for (ktest = 0; ktest < nRounds; ktest++){

		for (kRow = 0; kRow < numRows; kRow++){
			for (kAddress1 = 1; kAddress1 < numCols; kAddress1++){
				index1 = PrelMCUSummary[kRow][kAddress1][ktest];
			}


		}
		
	}

}

#endif
