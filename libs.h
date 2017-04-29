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
 * Función de apertura de fichero
 */


//libera memoria de una matriz**[][]
void liberaInt(int32_t** matrix, int rows){
    int i;
    for (i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}
void liberaUint(uint32_t** matrix, int rows){
    int i;
    for (i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}
//libera memoria matriz***[][][]
void libera3D(int32_t*** matrix, int row,int col){
    int i, j;
    for (int i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}


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

void sort(int32_t* A, int n){
	int min, i, j, aux;
	for (i = 0; i < n - 1; i++){
		min = i;
		for (j = i + 1; j < n; j++)
			if (A[min] > A[j])
				min = j;

		aux = A[min];
		A[min] = A[i];
		A[i] = aux;
	}
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
	if (columnVector[*elems - 1] == 0xFFFFFFFF){
		i = 0;
		uint32_t* newColumnVector;
		newColumnVector = (uint32_t*)malloc(sizeof(uint32_t)*(*elems));
		while (columnVector[i] != 4294967295){
			newColumnVector[i] = columnVector[i];
			i++;
		}
		*elems = i;
		//Reallocating memory
		newColumnVector = (uint32_t*)realloc(newColumnVector, sizeof(uint32_t)*i);
		return newColumnVector;
	}
	else {
		return columnVector;
	}
}

/*
* This function returns a vector with variable lenght indicating the positions
* of the bits where a known word differs from the pattern
*/
int32_t* flipped_bits(int32_t word, int32_t pattern, int wordWidth, int32_t* corruptedBitsLength){
	int k, nFound = 0, count = 0;
	int32_t* result = calloc(wordWidth, sizeof(int32_t));
	for (k = 0; k < wordWidth; k++) {
		result[k] = 1;
		result[k] *= (wordWidth + 1);
	}


	int32_t bitflips = word ^ pattern;
	for (k = 0; k < wordWidth; k++){
		if ((bitflips % 2) == 1) {
            result[nFound] = k;
            nFound++;
        }
		bitflips = bitflips >> 1;
	}
	//return result[1:findlast(result.!=(wordwidth+1))]
	// Recorremos el vector desde la última posición para encontrar el primero desde el final,
	// que cumple result != wordWith+1, realloc de result, para devolver vector con nuevo tamaño

	//como la matriz esta rellena inicialmente de wordwidth+1 = 9 y y los elementos encontrados
	//se introducen desde la posicion[0] hasta nFounds, no hace falta recorrer el vector.
	//se devuelven los nFound primeros del vector
	//si nFound = 7-> se devuelve result[0] hasta result[8]

	/*k = nFound-1;
	while ((k > -1) || (result[k] == wordWidth + 1)) {
		count++;
        k--;
	}*/
    if (nFound > 0) {
        result = realloc(result, nFound*sizeof(int32_t));
    }
    *corruptedBitsLength = nFound;
	return result;
}

uint32_t** transposed_matrix(uint32_t** A, int NumRow, int NumCol){
	uint32_t** trans;
	int a = 0;
	trans = malloc(NumRow*sizeof(uint32_t *));
	for (a = 0; a < NumRow; a++)
		trans[a] = malloc(NumCol*sizeof(uint32_t));

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
    int i, j;
    // DVmatrix será la matriz a devolver, tamaño = elems * elems
    uint32_t** DVmatrix = calloc((elems),sizeof(uint32_t *));
    for (i = 0; i < elems; i++){
        DVmatrix[i] = calloc((elems),sizeof(uint32_t));
    }
    
    // DVmatrixOP es una matriz para realizar las operaciones
    uint32_t** DVmatrixOP = calloc((elems),sizeof(uint32_t *));
    for (i = 0; i < elems; i++){
        DVmatrixOP[i] = calloc((elems),sizeof(uint32_t));
    }
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            DVmatrixOP[i][j] = addresses[i];
        }
    }
    
    // DVmatrixTrans matriz transpuesta para realizar las operaciones
    uint32_t** DVmatrixTrans = transposed_matrix(DVmatrixOP, elems, elems);

	if (strcmp(op, "xor") == 0){
		for (i = 0; i < elems; i++){
            for (j = 0; j < elems; j++){
				DVmatrix[i][j] = DVmatrixOP[i][j] ^ DVmatrixTrans[i][j];
			}
        }
	}
	else if (strcmp(op, "pos") == 0){
		for (i = 0; i < elems; i++){
            for (j = 0; j < elems; j++){
                 DVmatrix[i][j] = abs(DVmatrixOP[i][j] - DVmatrixTrans[i][j]);
			}
		}
	}
	else {
		printf("\n A problem creating the DV set. Operation not recognized\n");
		EXIT_FAILURE;
	}
    
    // Se deben eliminar los elementos de la triangular inferior
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            if ((i == j) || (j < i)) {
                DVmatrix[i][j] = 0;
            }
        }
    }

	/* Now, we will incorporate the information about the R-W cycles in the same
	* way. We create a matrix with the vector, traspose it and look for coincidence
	* element to element. Later, the matrix with True/False elements is converted
	* into integer and multiplied element to element with the DVmatrix. */
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            DVmatrixOP[i][j] = RWcycles[i];
        }
    }
    
    DVmatrixTrans = transposed_matrix(DVmatrixOP, elems, elems);
    // coincidence es una matriz para descubrir los elementos que coinciden en las dos matrices de operación
    uint32_t** coincidence = calloc((elems),sizeof(uint32_t *));
    for (i = 0; i < elems; i++){
        coincidence[i] = calloc((elems),sizeof(uint32_t));
    }
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            if (DVmatrixOP[i][j] == DVmatrixTrans[i][j]) {
                coincidence[i][j] = 1;
            } else {
                coincidence[i][j] = 0;
            }
        }
    }
    
    //Se deben eliminar los elementos de la triangular inferior
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            if ((i == j) || (j < i)) {
                coincidence[i][j] = 0;
            }
        }
    }

    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            DVmatrix[i][j] = DVmatrix[i][j] * coincidence[i][j];
        }
    }

	/* The same for the block division. We multiply the address vector by 2^bits4blocks/LN
	* and round to the floor integer. Traspose and multiply. */
	int32_t* addressblocks = calloc(elems, sizeof(int32_t));
	for (i = 0; i < elems; i++) {
		addressblocks[i] = floor(addresses[i] / ln);
	}
    
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            DVmatrixOP[i][j] = addressblocks[i];
        }
    }
    
    DVmatrixTrans = transposed_matrix(DVmatrixOP, elems, elems);
    
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            if (DVmatrixOP[i][j] == DVmatrixTrans[i][j]) {
                coincidence[i][j] = 1;
            } else {
                coincidence[i][j] = 0;
            }
        }
    }
    
    //Se deben eliminar los elementos de la triangular inferior
    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            if ((i == j) || (j < i)) {
                DVmatrix[i][j] = 0;
            }
        }
    }

    for (i = 0; i < elems; i++) {
        for (j = 0; j < elems; j++) {
            DVmatrix[i][j] = DVmatrix[i][j] * coincidence[i][j];
        }
    }

    liberaUint(coincidence, elems);
    liberaUint(DVmatrixOP, elems);
    liberaUint(DVmatrixTrans, elems);
	return DVmatrix;
}

/*
* This function creates the DVvector cotaining the elements of the triu in
* the DVmatrix. elems is the vector's size.
* Returns that vector.
*/
uint32_t* create_DVvector(uint32_t** DVmatrix, int elems, int* newTam){
	int ndv = 0.5*elems*(elems - 1);
    // Vector con los elementos de la triangular superior de la matriz seleccionada
    uint32_t* DVvector = (uint32_t*)calloc(ndv, sizeof(uint32_t));

	/* The elements of the matrix are recovered. Elements are added fixing the
	* column and going down. */
	int index = 0, cont = 0;
	int kcol, krow;
    
    for (kcol = 1; kcol < elems; kcol++) {
        for (krow = 0; krow < kcol; krow++) {
            if (DVmatrix[krow][kcol] != 0) {
                DVvector[index] = DVmatrix[krow][kcol];
                index++;
            }
        }
    }
    DVvector = realloc(DVvector, index*sizeof(uint32_t));
    *newTam = index;
	return DVvector;
}

/* Copia de matrices */
uint32_t*** copyOfMatrix(uint32_t** matrix, uint32_t*** matrixbackup, int elems, int ktest){
	int i, j;
	for (i = 0; i < elems; i++){
        for (j = 0; j < elems; j++) {
            matrixbackup[i][j][ktest - 1] = matrix[i][j];
        }
	}
	return matrixbackup;
}

uint32_t** copyOfMatrix3to2(uint32_t*** matrix, uint32_t** matrixbackup, int rows, int cols, int ktest){
	int i, j;
	for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++) {
            matrixbackup[i][j] = matrix[i][j][ktest];
        }
	}
    return matrixbackup;
}

/* Copia de vector a matriz */
void copyOfVector(uint32_t* vector, uint32_t** vectorbackup, int ndv, int ktest){
	int i;
	for (i = 0; i < ndv; i++){
		vectorbackup[i][ktest - 1] = vector[i];
	}
}

/*
* This function returns an array which has on it the count of each element on the histogram
*/
int32_t countsOfElems(int32_t** histogram, int maxValue, int col, long int LN){
	int i;
    int32_t repetitionstmp = 0;
	for (i = 0; i < LN; i++) {
		if (histogram[i][col] > 0 && histogram[i][col] < maxValue) {
			repetitionstmp++;
		}
	}
	return repetitionstmp;
}

int32_t* counts(int32_t** matrix, int matrixLength, int col, int maxValue){
    int i, j, count = 0;
    int32_t* result = calloc(maxValue+1, sizeof(int32_t));
    for (i = 0; i < maxValue+1; i++) {
        for (j = 0; j < matrixLength; j++) {
            if (matrix[j][col] == i) {
                count++;
            }
        }
        result[i] = count;
        count = 0;
    }
    return result;
}

/*
* This function looks for the number of appearances of an element
* vector[position]. nvd is the vector's size.
* Return the number of appearances of that element in the vector.
* Ex: [2, 1, 2, 3, 2] -> 3 1 3 1 3
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
	int k = 0, sum = 0;
	for (k = 0; k < ndv; k++) {
		sum = lookForElem(vector, ndv, k);
		histogram[vector[k]] = sum;
        sum = 0;
	}
	return histogram;
}

long long binomial(int n, int k){
    long long ans = 1;
    int j = 1;
    if (k > (n-k)) {
        k = n - k;
    }
    for(; j<=k; j++, n--){
        if(n % j == 0){
            ans *= n/j;
        } else {
            if(ans % j == 0){
                ans = ans/j * n;
            } else {
                ans = (ans*n) / j;
            }
        }
    }
    return ans;
}


//concatena 2 matrices con mismo numero de columnas.
int32_t** vcat(int32_t ** matrizA, int rowA, int colA, int32_t** matrizB, int rowB, int colB){
	int numRow = rowA + rowB;
	int32_t** matrizC = (int32_t **)malloc((rowA + rowB) *sizeof(int32_t *));
	int i = 0, j = 0;
	for (i = 0; i < numRow; i++)
		matrizC[i] = (int32_t *)malloc(colA*sizeof(int32_t));

	for (i = 0; i < rowA; i++){
		for (j = 0; j < colA; j++) {
			matrizC[i][j] = matrizA[i][j];
		}
	}
	for (i = rowA; i < numRow; i++){
		for (j = 0; j < colB; j++){
			matrizC[i][j] = matrizB[i][j];
		}
	}
	return matrizC;

}



/* This function just implements the theoretical expresions to determine expected
* repetitions.  Variable names are meaningful. */
long double expectedRepetitions(int m, long int LN, int ndv, char* operation){
	double result = 0, logresulttmp = 0, i;

	if (strcmp(operation, "pos") == 0){
		//Difference of values
		if (fabs((double)ndv / (double)LN) < 0.1){
			logresulttmp = (double)((1 - m)*log(LN) + m*log(2) - log(m + 1) + log(binomial((long long int)ndv, m)));
			result = exp(logresulttmp);
			result = result*(1 + (3 * (m)-2 * (int64_t)(ndv)) / (LN + 1))*(1 + 2 * (ndv - m) / (LN - 2) / (m + 2));
		} else{
			//          OJOOOOOOOOOOOOOOOOOOOOOOOO  EL COLLECT...
			/*pk=2/LN/(LN+1)*(LN+1-collect(1:1:LN));
			result = sum(binomial(BigInt(ndv), m)*(pk.^m).*(1-pk).^(ndv-m));*/

			for (i = 0; i < LN; i++){
				double pk = 2 / LN / (LN + 1)*(LN + 1 - i);
				result += (binomial((long long int)(ndv), m)*(pow(pk, m))*pow((1 - pk), (ndv - m)));
			}
		}
	}
	else if (strcmp(operation, "xor") == 0){
		// XOR
		logresulttmp = (double)(log(binomial((long long int)(ndv), m)) + (ndv - m)*log(LN - 1) - (ndv - 1)*log(LN));
        result = (exp(logresulttmp));
        
	}
	else{
		//Badly specified function.
		printf("\n\tAn unknown operation to calculate expected repetitions. Returning a strange value\n");
		result = -1;
	}
	return result;
}

/*
 // This function allows calculating the lowest number of repetitions from which
	// the expected number is below threshold.
	// RepInHistVector is a column vector containing the statistics.It is assumed
	// that the first element correspond to 0 repetitions.
 */
int32_t excessiveRepetitions(int32_t** repInHistVector, int repInHistVectorLength, int repInHistVectorCol, long int LN, char* operation, double threshold){
	int32_t nthreshold = 0, ndv = 0, k, i;
	int32_t* occurrenceIndex = calloc(repInHistVectorLength, sizeof(int32_t));
	for (i = 0; i < repInHistVectorLength; i++)
		occurrenceIndex[i] = i;

	for (i = 0; i < repInHistVectorLength; i++)
		ndv += ((int32_t)repInHistVector[i][repInHistVectorCol] * occurrenceIndex[i]);

	k = 0;
	if (strcmp(operation, "pos") == 0){
        if ((excessiveRepetitions(repInHistVector, repInHistVectorLength, repInHistVectorCol, LN, "xor", threshold) - 1) > 0){
			k = excessiveRepetitions(repInHistVector, repInHistVectorLength, repInHistVectorCol, LN, "xor", threshold) - 1;
        } else {
            k = 0;
        }
    } else {
        k = 0;
    }
	// This is a trick to speed up "pos" calculations. Instead of starting to
	// iterate from 0, we calculate the threshold for XOR, which is always the
	// lowest and it is really easy to determine.

    for (i = k; i < repInHistVectorLength; i++) {
        double o = expectedRepetitions(i, LN, ndv, operation);
        int lol;
        if (expectedRepetitions(i, LN, ndv, operation) < threshold) {
            nthreshold = i;
            break;
        }
    }
    free(occurrenceIndex);
	return nthreshold;
}

/*
* This function lists the element anomalously repeated in the histogram.
* occurrences is a column vector [V(0)...V(K)] where V(K) indicates the number
* of elements that are repeated K times.
* nthreshold indicates the limit value from which repetitions should not
* be expected in only-SBU experiments.
*/
int32_t** find_anomalies_histogram(uint32_t** histogram, int32_t histogramLenght, int histogramCol, int32_t** occurrences, int32_t occurrencesLenght, int occurrencesCol, int32_t nthreshold, int* newLenght, int* n_anomalous_values){ //n_anomalous_values por referencia porque necesitamos el valor
	/* The total number of elements anomalously repeated. We use +1 since the first
	* element of the vector is related to 0 repetitions.
	*/
	int i, index = 0, k, k2;
	int32_t nAnomalies = 0;

	for (i = nthreshold; i < occurrencesLenght; i++) {
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
    *newLenght = nAnomalies;
    for (k = highestAnomaly; k > nthreshold-1; k--) {
		/* I've prefered to implement this solution instead of using the native function
		* "find()" since I believe this solution is faster as it includes breaks once
		* the number of anomalous occurrences are achieved. */
		if (occurrences[k][occurrencesCol] != 0){
			int n_occurrence_value = 0;
			int occurrence_value = k;
			*n_anomalous_values += 1;
			for (k2 = 0; k2 < histogramLenght; k2++) {
				if (histogram[k2][histogramCol] == occurrence_value){
					anomalies[index][0] = k2 + 1;
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
int32_t** extract_some_critical_values(int32_t** anomalous_repetitions, int32_t anomalous_repetitionsLenght, int* newLenght, int* n_anomalous_repetitions, int n){
    // No sabemos si n = 1 forever
	int observed_values = 1, nrow = 1, k;
	if (*n_anomalous_repetitions <= n) {
		/*int32_t** extracted_values = calloc(anomalous_repetitionsLenght, sizeof(int32_t*));
		for (k = 0; k < anomalous_repetitionsLenght; k++) {
			extracted_values[k] = calloc(2, sizeof(int32_t));
		}
		//extracted_values = anomalous_repetitions; //Copia de la matriz en otra matriz
		for (k = 0; k < anomalous_repetitionsLenght; k++) {
			extracted_values[k][0] = anomalous_repetitions[k][0];
			extracted_values[k][1] = anomalous_repetitions[k][1];
		}
		return extracted_values;*/
        *newLenght = anomalous_repetitionsLenght;
        return anomalous_repetitions;
	} else {
		for (k = 1; k < anomalous_repetitionsLenght; k++){
			if (anomalous_repetitions[k][1] != anomalous_repetitions[k - 1][1]) {
				if (observed_values == n) {
					nrow = k;
				}
				observed_values++;
			}
		}
        *newLenght = nrow;
		int32_t** extracted_values = calloc(*newLenght, sizeof(int32_t*));
		for (k = 0; k < *newLenght; k++) {
			extracted_values[k] = calloc(2, sizeof(int32_t));
		}
		//extracted_values = anomalous_repetitions[1:nrow,:];
		for (k = 0; k < *newLenght; k++) {
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

	for (i = 0; i < critical_valuesLenght; i++) {
		value = critical_values[i][0];
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
			}
			else if (thesummarytemp[j - 1] < thesummary[addressRowMin][i]){ // Si el valor último menor
				thesummarytemp[j] = thesummary[addressRowMin][i];
				j++;
			}
			else if (thesummarytemp[j - 1] > thesummary[addressRowMin][i]){ // Si el valor último mayor, se abre hueco
				x = j - 1;
				while (thesummarytemp[x] > thesummary[addressRowMin][i]) {
					thesummarytemp[x + 1] = thesummarytemp[x]; // Se desplaza uno
					x--;
				}
				thesummarytemp[x + 1] = thesummary[addressRowMin][i]; // Se copia en su hueco correspondiente
			}
		}
	}
    for (i = 0; i < thesummaryCols; i++) {
		if (thesummary[addressRowMax][i] != 0) {
			if (thesummarytemp[0] == 0) { // Esta vacío
				thesummarytemp[j] = thesummary[addressRowMax][i];
				j++;
			}
			else if (thesummarytemp[j - 1] < thesummary[addressRowMax][i]){ // Si el valor último menor
				thesummarytemp[j] = thesummary[addressRowMax][i];
				j++;
			}
			else if (thesummarytemp[j - 1] > thesummary[addressRowMax][i]){ // Si el valor último mayor, se abre hueco
				x = j - 1;
				while (thesummarytemp[x] > thesummary[addressRowMax][i]) {
					thesummarytemp[x + 1] = thesummarytemp[x]; // Se desplaza uno
					x--;
				}
				thesummarytemp[x + 1] = thesummary[addressRowMax][i]; // Se copia en su hueco correspondiente
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
    int count = 0, i, j, k, x, xx;
    int32_t* nonzerovectorelements = calloc(davmatrixRows*davmatrixCols, sizeof(int32_t));
    for (i = 0; i < davmatrixRows; i++) {
        for (j = 0; j < davmatrixCols; j++) {
            if (davmatrix [i][j] != false) {
                nonzerovectorelements[count] = j * davmatrixCols + i;
                count++; //Sabemos cuantos elementos, pero alguno puede estar repetido
            }
        }
    }
    // Condicion si count = 0 con break de la función?
    
    nonzerovectorelements = realloc(nonzerovectorelements, count*sizeof(int32_t));
    
    //Now, transform the elements in pairs creating a matrix with 2 columns.
    int32_t** relatedpairs = calloc(count, sizeof(int32_t*));
    for (i = 0; i < count; i++) {
        relatedpairs[i] = calloc(2, sizeof(int32_t));
    }
    
    /*    @simd for k1 = 1:length(nonzerovectorelements)
     # The first operation calculates the column
     @inbounds relatedpairs[k1,:] = [rem(nonzerovectorelements[k1], NRows) div(nonzerovectorelements[k1],NRows)+1]
     # REM means "remainder of division", DIV integer division.
     end*/
    for (i = 0; i < count; i++) {

        if((nonzerovectorelements[i] % davmatrixRows != 0)){
        relatedpairs[i][0] = nonzerovectorelements[i] % davmatrixRows+1;
        }else{
            relatedpairs[i][0] = nonzerovectorelements[i] % davmatrixRows;
        }
        
        relatedpairs[i][1] = (nonzerovectorelements[i] / davmatrixRows)+1;
    }
    
    free(nonzerovectorelements);
    /* Good, we have created a matrix containing the related pairs. Now, we must
     # group them in larger events.
     
     # First step: A matrix initilizing the events. Perhaps too large, but it is better
     # to have spare space.
     */
    int32_t** thesummary = calloc(count, sizeof(int32_t*));
    int32_t thesummaryRows = count, thesummaryCols;
    if (count < NMaxAddressesInEvent) { // Mínimo es newCount
        thesummaryCols = count;
        for (i = 0; i < count; i++) {
            thesummary[i] = calloc(count, sizeof(int32_t));
        }
    } else { // Mínimo es NMaxAddressesInEvent
        thesummaryCols = NMaxAddressesInEvent;
        for (i = 0; i < count; i++) {
            thesummary[i] = calloc(NMaxAddressesInEvent, sizeof(int32_t));
        }
    }
    //Now, we save the first pair in the first column:
    *nTotalEvents = 0;
    int totalEvents = 0, thelargestMCU = 0;
    int32_t firstAddress, secondAddress, firstAddressRow = 0, secondAddressRow = 0;
    for (k = 0; k < count; k++) {
        //There are several cases. Let us list them in order.
        //Positions are extrated from the pairs.
        firstAddress = relatedpairs[k][0];
        secondAddress = relatedpairs[k][1];
        //TENEMOS QUE MIRAR SI firstAddress y secondAddress están en thesummary
        bool firstAddressFound = false, secondAddressFound = false;
        for (i = 0; i < thesummaryRows; i++) {
            for (j = 0; j < thesummaryCols; j++) {
                if ((thesummary[i][j] == firstAddress) && (firstAddressFound == false)) {
                    firstAddressFound = true;
                    firstAddressRow = i;
                }
                if ((thesummary[i][j] == secondAddress) && (secondAddressFound == false)) {
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
            thesummary[totalEvents][0] = relatedpairs[k][0];
            thesummary[totalEvents][1] = relatedpairs[k][1];
            totalEvents++;
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
                    thesummary[addressRowMax][i] = 0;
                }
                //Also, as two events have been merged, the number of total events is lower:
                totalEvents--;
            }
        }
    }
    
    liberaInt(relatedpairs, count);
    /*
     #Really good. Now, we will separate the cells with information from those with
     #zeros.
     # The number of rows is easy to calculate: NTotalEvents. Concerning the other
     # element:
     */
    int ind = 0, sum = 0;
    *largestMCU = 0;
    for(j = 2; j < thesummaryCols; j++) {
        for (i = 0; i < thesummaryRows; i++) {
            if(thesummary[i][j] != 0){
                if(j > ind){
                    ind = j;
                }
                sum++;
            }
        }
        if (sum == 0) {
            thelargestMCU = j;
            break;
        } else {
            thelargestMCU = ind+1;
        }
    }
    int32_t** result = calloc(totalEvents, sizeof(int32_t*));
    for (i = 0; i < totalEvents; i++) {
        result[i] = calloc(thelargestMCU, sizeof(int32_t));
    }
    for (x = 0; x < totalEvents; x++) {
        for (xx = 0; xx < thelargestMCU; xx++) {
            result[x][xx] = thesummary[x][xx];
        }
    }
    *nTotalEvents = totalEvents;
    *largestMCU = thelargestMCU;

    liberaInt(thesummary, thesummaryRows);
    return result;
}

// Devuelve el número de elementos que son diferentes de el número a buscar en una matriz
int findNotEqualMatrix(int32_t** matrix, int matrixRows, int matrixCols, int search){
    int count = 0, i, j;
    for (i = 0; i < matrixRows; i++) {
        for (j = 0; j < matrixCols; j++) {
            if (matrix[i][j] != search) {
                count++;
            }
        }
    }
    return count;
}

// Devuelve el número de elementos que son iguales de el número a buscar en una matriz
int findEqualMatrix(int32_t** matrix, int matrixRows, int matrixCols, int search){
    int count = 0, i, j;
    for (i = 0; i < matrixRows; i++) {
        for (j = 0; j < matrixCols; j++) {
            if (matrix[i][j] == search) {
                count++;
            }
        }
    }
    return count;
}

int findEqualVector(int32_t* vector, int length, int search){
    int count = 0, i;
    for (i = 0; i < length; i++) {
        if (vector[i] == search) {
            count++;
        }
    }
    return count;
}

// Devuelve el número de elementos que son diferentes de el número a buscar, en una fila o una columna de una matriz
int findNotEqual(int32_t** matrix, int matrixRows, int matrixCols, int selectedRow, int selectedCol, int search){
    int count = 0, i;
    if (selectedCol == -1) { // Búsqueda en columnas de fila seleccionada
        for (i = 0; i < matrixCols; i++) {
            if (matrix[selectedRow][i] != search) {
                count++;
            }
        }
    } else { // Búsqueda en filas de columna seleccionada
        for (i = 0; i < matrixRows; i++) {
            if (matrix[i][selectedCol] != search) {
                count++;
            }
        }
    }
    return count;
}

// Devuelve el número de elementos que son iguales de el número a buscar, en una fila o una columna de una matriz
int findEqual(int32_t** matrix, int matrixRows, int matrixCols, int selectedRow, int selectedCol, int search){
    int count = 0, i;
    if (selectedCol == -1) { // Búsqueda en columnas de fila seleccionada
        for (i = 0; i < matrixCols; i++) {
            if (matrix[selectedRow][i] == search) {
                count++;
            }
        }
    } else { // Búsqueda en filas de columna seleccionada
        for (i = 0; i < matrixRows; i++) {
            if (matrix[i][selectedCol] == search) {
                count++;
            }
        }
    }
    return count;
}

/*
 *This is a simple function to cut files and rows of redundant zeros in 1 or 2
 *dimensions matrix.
 */

int32_t** cutZerosFromArray(int32_t** matrix, int32_t matrixRows, int32_t matrixCols, int* newRows, int* newCols, bool* isMatrix, bool* isVector){
	int i, j = 0, countIni = 0, countFin = 0;
    isMatrix = false;
    isVector = false;
	*newRows = 0;
	*newCols = 0;

 /*   for (i = 0; i < matrixRows; i++) {
     for (j = 0; j < matrixCols; j++) {
     printf("%d ", matrix[i][j]);
     }
        printf("\n");
     }
     printf("\n");*/

	if (matrixRows > 1 && matrixCols > 1) {
        isMatrix = true;
		// this is a classical matrix.
		// Establecemos las esquinas de la matriz
		int firstRow = 0;
		int lastRow = matrixRows;
		int firstCol = 0;
		int lastCol = matrixCols;
		for (i = 0; i < matrixRows; i++) { // Desde la primera fila
			if (findNotEqual(matrix, matrixRows, matrixCols, i, -1, 0) != 0) {
				break;
			}
			firstRow++;
		}
		for (i = matrixRows-1; i >= 0; i--) { // Desde la última fila
			if (findNotEqual(matrix, matrixRows, matrixCols, i, -1, 0) != 0) {
				break;
			}
			lastRow--;
		}
		for (i = 0; i < matrixCols; i++) { // Desde la primera columna
			if (findNotEqual(matrix, matrixRows, matrixRows, -1, i, 0) != 0) {
				break;
			}
			firstCol++;
		}
		for (i = matrixCols-1; i >= 0; i--) { // Desde la última columna
			if (findNotEqual(matrix, matrixRows, matrixCols, -1, i, 0) != 0) {
				break;
			}
			lastCol--;
		}
		// Calculamos el nuevo tamaño de matriz y la declaramos
		*newRows = (lastRow - firstRow);
		*newCols = (lastCol - firstCol);
		int32_t** result = calloc(*newRows, sizeof(int32_t*));
		for (i = 0; i < *newRows; i++) {
			result[i] = calloc(*newCols, sizeof(int32_t));
		}
		// Copiamos la antigua en la nueva
		for (i = firstRow; i < lastRow; i++) {
			for (j = firstCol; j <= lastCol; j++) {
				result[i][j] = matrix[i][j];
			}
		}
		return result;
	}
	else if (matrixRows == 1 || matrixCols == 1){
		// Upss, this is a vector. Be careful.
        isVector = true;
		if (matrixRows == 1) { // Única fila, recorrer en horizontal
			while (matrix[countIni] == 0 && countIni < matrixCols) {
				countIni++; // Incrementamos el contador mientras existan 0s, para conocer cuántos hay delante de los valores
			}
			i = matrixCols;
			while (matrix[i] == 0 && i >= 0) {
				countFin++; // Incrementamos el contador mientras existan 0s, para conocer cuántos hay detrás de los valores
				i--;
			}
			// Una vez acotados los lados, se crea la variable y se copian los datos con las nuevas dimensiones
			*newCols = matrixCols - (countIni + countFin);
			int32_t** result = calloc(*newCols, sizeof(int32_t));
			for (i = countIni; i < (matrixCols - countFin); i++) {
				result[j][0] = matrix[0][i];
				j++;
			}
			return result;
		}
		else { // Única columna recorrer en vertical
			while (matrix[countIni] == 0 && countIni < matrixRows) {
				countIni++; // Incrementamos el contador mientras existan 0s, para conocer cuántos hay delante de los valores
			}
			i = matrixRows;
			while (matrix[i] == 0 && i >= 0) {
				countFin++; // Incrementamos el contador mientras existan 0s, para conocer cuántos hay detrás de los valores
				i--;
			}
			// Una vez acotados los lados, se crea la variable y se copian los datos con las nuevas dimensiones
			*newCols = matrixRows - (countIni + countFin);
			int32_t** result = calloc(*newCols, sizeof(int32_t));
			for (i = countIni; i < (matrixRows - countFin); i++) {
				result[j][0] = matrix[i][0];
				j++;
			}
			return result;
		}
	}
	else {
		printf("\tMatriz dimension different from 1 or 2. Exiting.");
	}
	return matrix;
}

/*
# A very specific function to reduce the size of the summaries with many
# involved rounds.
*/
int32_t*** cutZerosFromMCUsummary(int32_t*** MCUSummary, int* rows, int* cols, int dims){
	int i, j, z, maxRows = 0, maxCols = 0;
	int newRows = 0;
	int newCols = 0;
    int tmpMatrixRows = *rows, tmpMatrixCols = *cols;
    bool isMatrix;
    bool isVector;

	int32_t** dimensionsMCU = calloc(dims, sizeof(int32_t*));
	for (i = 0; i < dims; i++) {
		dimensionsMCU[i] = calloc(2, sizeof(int32_t));
	}
	for (i = 0; i < dims; i++) { // Cada matriz de cada dimensión separar para cutzeros?
        // MCU matriz para cada test
        int32_t** tmpMatrix = calloc(tmpMatrixRows, sizeof(int32_t*));
        for (j = 0; j < tmpMatrixRows; j++) {
            tmpMatrix[j] = calloc(tmpMatrixCols, sizeof(int32_t));
        }
        for (j = 0; j < tmpMatrixRows; j++) {
            for (z = 0; z < tmpMatrixCols; z++) {
                tmpMatrix[j][z] = MCUSummary[j][z][i];
            }
        }
		cutZerosFromArray(tmpMatrix, tmpMatrixRows, tmpMatrixCols, &newRows, &newCols, &isMatrix, &isVector);
		dimensionsMCU[i][0] = newRows;
		dimensionsMCU[i][1] = newCols;
		// Es necesario saber los mayores valores, para adecuar la matriz total
		if (newRows > maxRows) {
			maxRows = newRows;
		}
		if (newCols > maxCols) {
			maxCols = newCols;
		}
        liberaInt(tmpMatrix, tmpMatrixRows);
	}
	// Una vez guardadas en la matriz de valores todos los nuevos tamaños, se copia la matriz en una matriz nueva
	int32_t*** simpleMCUSummary = calloc(maxRows, sizeof(int32_t**));
	for (i = 0; i < maxRows; i++) {
		simpleMCUSummary[i] = calloc(maxCols, sizeof(int32_t*));
		for (j = 0; j < maxCols; j++) {
			simpleMCUSummary[i][j] = calloc(dims, sizeof(int32_t));
		}
	}
	//Ojo, apuntado papel
	for (i = 0; i < maxRows; i++) {
		for (j = 0; j < maxCols; j++) {
			for (z = 0; z < dims; z++) {
				simpleMCUSummary[i][j][z] = MCUSummary[i][j][z];
			}
		}
	}
    *rows = maxRows;
    *cols = maxCols;
	//liberaInt(dimensionsMCU,2);
	return simpleMCUSummary;
}

/*
 * Función que recibe un vector y una matriz, devuelve un vector con los elementos que estén en el vector inicial, 
 * pero que no estén en la matriz. Si no, devuelve array vacío
 */
int32_t* setdiffTraceRule(int32_t* vector, int vectorLenght, int32_t** matrix, int matrixRows, int matrixCols, int* newLenght){
    int i, x, y, cont = 0;
    bool find = false;
    int32_t* newVector = calloc(vectorLenght, sizeof(int32_t));
    for (i = 0; i < vectorLenght; i++) {
        for (x = 0; x < matrixRows; x++) {
            for (y = 0; y < matrixCols; y++) {
                if (vector[i] == matrix[x][y]) {
                    find = true;
                }
            }
        }
        if (find == false) {
            newVector[cont] = vector[i];
            cont++;
        }
        find = false;
    }
    
    newVector = realloc(newVector, cont*sizeof(int32_t));
    *newLenght = cont;
    return newVector;
}

//UNION de max 3 vectores, no ordena
int32_t* unionVec(int32_t* v1, int v1Elems, int32_t* v2, int v2Elems, int32_t* v3, int v3Elems, int* tam){
    int32_t* unionVec;
    int i, j, cont = 0;
    if (v3 == NULL) { //Solo unión de dos vectores
        unionVec = calloc((v1Elems+v2Elems), sizeof(int32_t));
        for (i = 0; i < v1Elems; i++) { // Se copia el primer vector eliminando los elementos repetidos
            if (i == 0) {
                unionVec[i] = v1[i];
                cont++;
            } else {
                j = 0;
                while (v1[i] != unionVec[j] && j < cont) {
                    j++;
                }
                if (j == cont) { // No se ha encontrado elemento igual
                    unionVec[cont] = v1[i];
                    cont++;
                }
            }
        }
        for (i = 0; i < v2Elems; i++) {
            j = 0;
            while (j < cont && v2[i] != unionVec[j]) {
                j++;
            }
            if (j == cont) { // No se ha encontrado elemento igual
                unionVec[cont] = v2[i];
                cont++;
            }
        }
        unionVec = realloc(unionVec, cont*sizeof(int32_t));
    }else { // Unión de 3 vectores
        unionVec = calloc((v1Elems+v2Elems+v3Elems), sizeof(int32_t));
        for (i = 0; i < v1Elems; i++) { // Se copia el primer vector eliminando los elementos repetidos
            if (i == 0) {
                unionVec[i] = v1[i];
                cont++;
            } else {
                j = 0;
                while (v1[i] != unionVec[j] && j < cont) {
                    j++;
                }
                if (j == cont) { // No se ha encontrado elemento igual
                    unionVec[cont] = v1[i];
                    cont++;
                }
            }
        }
        for (i = 0; i < v2Elems; i++) {
            j = 0;
            while (j < cont && v2[i] != unionVec[j]) {
                j++;
            }
            if (j == cont) { // No se ha encontrado elemento igual
                unionVec[cont] = v2[i];
                cont++;
            }
        }
        for (i = 0; i < v3Elems; i++) {
            j = 0;
            while (j < cont && v3[i] != unionVec[j]) {
                j++;
            }
            if (j == cont) { // No se ha encontrado elemento igual
                unionVec[cont] = v3[i];
                cont++;
            }
        }
        unionVec = realloc(unionVec, cont*sizeof(int32_t));
    }
    *tam = cont;
    return unionVec;
}

/*
# An alternative implementation of the trace rule.

### Anomal XOR values is used to avoid the repetition of elements. It is just
### XORExtractedValues[1,:]
*/
int32_t* traceRule(int32_t** xorDVtotalrepetitions,int xorDVtotalrepetitionsLenght, uint32_t** totalDVhistogram, int32_t** anomalXORvalues, int anomalXORvaluesRows, long int LN, int* traceRuleLong){
	// Anomal XOR values is used to avoid the repetition of elements. It is just
    int i, j, z, index;
    int32_t nBits = (log(LN + 1) / log(2)); // ejemplo nBits = 24
    
	// Elements with 1-trace
	printf("\n\tInvestigating Trace = 1: ");
	//Crea un array de 24 elementos 0 -> 24-1 con valores potencia de 2
	int32_t* candidatesT1 = calloc(nBits, sizeof(int32_t));
	for (i = 0; i < nBits; i++){
		candidatesT1[i] = pow(2, i);
	}

	int32_t nT1Threshold = excessiveRepetitions(xorDVtotalrepetitions, xorDVtotalrepetitionsLenght, 1, LN, "xor", randomnessThreshold);
    
    int32_t* selectedT1 = calloc(nBits, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT1Lenght = 0;
    for (i = 0; i < nBits; i++) { // Se recorre el totalDVhistogram, en las posiciones de los candidatos - 1, pq nuestro histg empieza 0
        if (totalDVhistogram[candidatesT1[i]-1][1] >= nT1Threshold) {
            selectedT1[sT1Lenght] = candidatesT1[i];
            sT1Lenght++;
        }
    }
    selectedT1 = realloc(selectedT1, sT1Lenght*sizeof(int32_t));
    free(candidatesT1);

    // Devolver array con los elementos que estén en selectedt1, pero que no estén en anomalxorvalues
    int nWinnersT1 = 0;
    int32_t* winnersT1 = setdiffTraceRule(selectedT1, sT1Lenght, anomalXORvalues, anomalXORvaluesRows, 2, &nWinnersT1);
    
    free(selectedT1);
    printf("%d candidate", nWinnersT1);
    if (nWinnersT1 != 1)
        printf("s");
    printf(" found.");
    
    /// Elements with 2-trace
    printf("\n\tInvestigating Trace = 2: ");
    int32_t candidatesT2Length = 0.5*nBits*(nBits-1);
    int32_t* candidatesT2 = calloc(candidatesT2Length, sizeof(int32_t));
    index = 0;
    for (i = 0; i < (nBits-1); i++) {
        for (j = i+1; j < nBits; j++) {
            candidatesT2[index] = (pow(2, i))+(pow(2, j)); // CandidatesT2[index]=2^k1+2^k2
            index++;
        }
    }
    
    int32_t nT2Threshold = excessiveRepetitions(xorDVtotalrepetitions, xorDVtotalrepetitionsLenght, 1, LN, "xor", randomnessThreshold);
   
    int32_t* selectedT2 = calloc(candidatesT2Length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT2Lenght = 0;
    for (i = 0; i < candidatesT2Length; i++) {
        if (totalDVhistogram[candidatesT2[i]-1][1] >= nT2Threshold) {
            selectedT2[sT2Lenght] = candidatesT2[i];
            sT2Lenght++;
        }
    }
    selectedT2 = realloc(selectedT2, sT2Lenght*sizeof(int32_t));
    free(candidatesT2);
    
    int nWinnersT2 = 0;
    int32_t* winnersT2 = setdiffTraceRule(selectedT2, sT2Lenght, anomalXORvalues, anomalXORvaluesRows, 2, &nWinnersT2);
    
    free(selectedT2);
    printf("%d candidate", nWinnersT2);
    if (nWinnersT2 != 1)
        printf("s");
    printf(" found.");
    
    /// Elements with 3-trace
    printf("\n\tInvestigating Trace = 3: ");
    int32_t candidatesT3Length = nBits*(nBits-1)*(nBits-2)/6;
    int32_t* candidatesT3 = calloc(candidatesT3Length, sizeof(int32_t));
    index = 0;
    for (i = 0; i < (nBits-2); i++) {
        for (j = i+1; j < (nBits-1); j++) {
            for (z = j+1; z < nBits; z++) {
                candidatesT3[index] = (pow(2, i))+(pow(2, j))+(pow(2, z)); // 2^k1+2^k2+2^k3
                index++;
            }
        }
    }
    int32_t* selectedT3 = calloc(candidatesT3Length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT3Lenght = 0;
    for (i = 0; i < candidatesT3Length; i++) {
        if (totalDVhistogram[candidatesT3[i]-1][1] >= nT2Threshold) {
            selectedT3[sT3Lenght] = candidatesT3[i];
            sT3Lenght++;
        }
    }
    selectedT3 = realloc(selectedT3, sT3Lenght*sizeof(int32_t));
    
    free(candidatesT3);
    
    int nWinnersT3 = 0;
    int32_t* winnersT3 = setdiffTraceRule(selectedT3, sT3Lenght, anomalXORvalues, anomalXORvaluesRows, 2, &nWinnersT3);
    printf("%d candidate", nWinnersT3);
    free(selectedT3);
    if (nWinnersT3 != 1)
        printf("s");
    printf(" found.");

    printf("\n\tEnd of search.\n");
    
    int tam = 0;
    int32_t* vectorUnion = unionVec(winnersT1, nWinnersT1, winnersT2, nWinnersT2, winnersT3, nWinnersT3, &tam);
    *traceRuleLong = tam;
    free(winnersT1);
    free(winnersT2);
    free(winnersT3);
    
    return vectorUnion;
}


int32_t*** propose_MCUs(char* op, uint32_t*** XORDVmatrix3D, int numRows, int numCols, int dim, uint32_t*** POSDVmatrix3D, int32_t** anomalXOR, int32_t anomalXORlength, int32_t** anomalPOS, int32_t anomalPOSlength, int* nAddressesInRound, int* rowss, int* colss, int* dimss){

	int i, j, z, ktest;
    *dimss = dim;
	int rows = ceil(numRows / 2 - 1);
	if (strcmp(op, "xor") == 0){
		int32_t*** XOR_summary = calloc(rows, sizeof(int32_t**));
		for (i = 0; i < rows; i++) {
			XOR_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < UnrealMCUsize; j++) {
				XOR_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}
		for (ktest = 0; ktest < dim; ktest++){
			int32_t** xordvmatrix = calloc(nAddressesInRound[ktest], sizeof(int32_t*));
			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				xordvmatrix[i] = calloc(nAddressesInRound[ktest], sizeof(int32_t));
			}
            copyOfMatrix3to2(XORDVmatrix3D, xordvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], ktest);

			bool** XORmarked_pairs = marking_addresses(xordvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], anomalXOR, anomalXORlength);
        
			//The following results only concern an operation.Just present for illustrative
			//operations.In practical use, only CMB should be used.
			int nTotalEvents = 0, largestMCUs = 0;
			int32_t** XORproposed_MCUs = agrupate_mcus(XORmarked_pairs, nAddressesInRound[ktest], nAddressesInRound[ktest], &nTotalEvents, &largestMCUs);

			//Let us save everything.
			for (i = 0; i < nTotalEvents; i++) {
				for (j = 0; j < largestMCUs; j++) {
					XOR_summary[i][j][ktest] = XORproposed_MCUs[i][j];
				}
			}
		}
        // However, as there are too many zeros, it is interesting to reshape the matrix.
        // Thus, some memory is saved.But we must be careful with the different dimensions
        // of every partial XY matrix in XYZ arrays.A simple but unelegant way is :

        int xorSummaryRows = rows;
        int xorSummaryCols = UnrealMCUsize;
        int32_t*** finalXOR_summary = cutZerosFromMCUsummary(XOR_summary, &xorSummaryRows, &xorSummaryCols, dim);
        *rowss = xorSummaryRows;
        *colss = xorSummaryCols;
       // libera3D(XOR_summary, rows, UnrealMCUsize);
        return finalXOR_summary;
	}
	else if (strcmp(op, "pos") == 0){
		int32_t*** POS_summary = calloc(rows, sizeof(int32_t**));
		for (i = 0; i < rows; i++) {
			POS_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < UnrealMCUsize; j++) {
				POS_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < dim; ktest++){
			int32_t** posdvmatrix = calloc(nAddressesInRound[ktest], sizeof(int32_t*));
			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				posdvmatrix[i] = calloc(nAddressesInRound[ktest], sizeof(int32_t));
			}
            copyOfMatrix3to2(POSDVmatrix3D, posdvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], ktest);

			bool** POSmarked_pairs = marking_addresses(posdvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], anomalPOS, anomalPOSlength);

			//The following results only concern an operation.Just present for illustrative
			//operations.In practical use, only CMB should be used.
			int nTotalEvents = 0, largestMCUs = 0;
			int32_t** POSproposed_MCUs = agrupate_mcus(POSmarked_pairs, nAddressesInRound[ktest], nAddressesInRound[ktest], &nTotalEvents, &largestMCUs);


			//Let us save everything.
			for (i = 0; i < nTotalEvents; i++) {
				for (j = 0; j < largestMCUs; j++) {
					POS_summary[i][j][ktest] = POSproposed_MCUs[i][j];
				}
			}
		}
        // However, as there are too many zeros, it is interesting to reshape the matrix.
        // Thus, some memory is saved.But we must be careful with the different dimensions
        // of every partial XY matrix in XYZ arrays.A simple but unelegant way is :
        int posSummaryRows = rows;
        int posSummaryCols = UnrealMCUsize;
        int32_t*** finalPOS_summary = cutZerosFromMCUsummary(POS_summary, &posSummaryRows, &posSummaryCols, dim);
        *rowss = posSummaryRows;
        *colss = posSummaryCols;
        libera3D(POS_summary,rows,UnrealMCUsize);
        return finalPOS_summary;
	}
	else{//op == CMB
		int32_t*** CMB_summary = calloc(rows, sizeof(int32_t**));
		for (i = 0; i < rows; i++) {
			CMB_summary[i] = calloc(UnrealMCUsize, sizeof(int32_t*));
			for (j = 0; j < UnrealMCUsize; j++) {
				CMB_summary[i][j] = calloc(dim, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < dim; ktest++){
			int32_t** xordvmatrix = calloc(nAddressesInRound[ktest], sizeof(int32_t*));
			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				xordvmatrix[i] = calloc(nAddressesInRound[ktest], sizeof(int32_t));
			}
            copyOfMatrix3to2(XORDVmatrix3D, xordvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], ktest);

			int32_t** posdvmatrix = calloc(nAddressesInRound[ktest], sizeof(int32_t*));
			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				posdvmatrix[i] = calloc(nAddressesInRound[ktest], sizeof(int32_t));

			}
            copyOfMatrix3to2(POSDVmatrix3D, posdvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], ktest);


			bool** XORmarked_pairs = marking_addresses(xordvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], anomalXOR, anomalXORlength);
			bool** POSmarked_pairs = marking_addresses(posdvmatrix, nAddressesInRound[ktest], nAddressesInRound[ktest], anomalPOS, anomalPOSlength);

			bool** XOR_POS_marked_pairs = calloc(nAddressesInRound[ktest], sizeof(bool*));
			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				XOR_POS_marked_pairs[i] = calloc(nAddressesInRound[ktest], sizeof(bool));
			}

			for (i = 0; i < nAddressesInRound[ktest]; i++) {
				for (j = 0; j < nAddressesInRound[ktest]; j++) {
					XOR_POS_marked_pairs[i][j] = XORmarked_pairs[i][j] | POSmarked_pairs[i][j];
				}
			}

			int nTotalEvents = 0, largestMCUs = 0;
			int32_t** CMBproposed_MCUs = agrupate_mcus(XOR_POS_marked_pairs, nAddressesInRound[ktest], nAddressesInRound[ktest], &nTotalEvents, &largestMCUs);

			//Let us save everything.
			for (i = 0; i < nTotalEvents; i++) {
				for (j = 0; j < largestMCUs; j++) {
					CMB_summary[i][j][ktest] = CMBproposed_MCUs[i][j];
                }
			}
		}
        // However, as there are too many zeros, it is interesting to reshape the matrix.
        // Thus, some memory is saved.But we must be careful with the different dimensions
        // of every partial XY matrix in XYZ arrays.A simple but unelegant way is :

        int cmbSummaryRows = rows;
        int cmbSummaryCols = UnrealMCUsize;
        int32_t*** finalCMB_summary = cutZerosFromMCUsummary(CMB_summary, &cmbSummaryRows, &cmbSummaryCols, dim);
        *rowss = cmbSummaryRows;
        *colss = cmbSummaryCols;
        //libera3D(CMB_summary, rows, UnrealMCUsize);
        return finalCMB_summary;

	}
    int32_t*** aux = NULL;
    return aux;
}

int32_t** condensate_summary(int32_t*** summary3D, int summary3DRows, int summary3DCols, int summary3DDims, int* nAddressesInRound, int* condensateRows, int* condensateCols){
	int i, j, k = 1, sum = 0;
	int32_t** condensateSummary = calloc(summary3DCols, sizeof(int32_t*));
	for (i = 0; i < summary3DCols; i++) {
		condensateSummary[i] = calloc(summary3DDims + 1, sizeof(int32_t));
        condensateSummary[i][0] = k;
        k++;
	}
	for (i = 1; i < summary3DDims + 1; i++) {
		// elems tendrá la longitud del nuevo vector
		int nEvents = 0;
		/*uint32_t* columnvector = (uint32_t*)malloc(sizeof(uint32_t)*(elems));
		for (k = 0; k < elems; ++k) {
			columnvector[k] = content[k][3*i-3];
		}
		uint32_t* addresses = extract_addressvector(columnvector, &elems);
        free(columnvector);
        free(addresses);*/
		//Ahora para cada test: Summary = Summary3D[:,:,ktest]
        int32_t** summary2D = calloc(summary3DRows, sizeof(int32_t*));
        for (k = 0; k < summary3DRows; k++) {
            summary2D[k] = calloc(summary3DCols, sizeof(int32_t*));
        }
        copyOfMatrix3to2(summary3D, summary2D, summary3DRows, summary3DCols, (i-1));

		for (k = summary3DCols-1; k > 0; k--) {
            nEvents = findNotEqual(summary2D, summary3DRows, summary3DCols, -1, k, 0);
            for (j = k; j < summary3DCols; j++) {
                sum += condensateSummary[j][i];
            }
            condensateSummary[k][i] = nEvents - sum;
            sum = 0;
		}
		// Apparently, the following solution is a bit faster.
        condensateSummary[0][i] = nAddressesInRound[i-1] - findNotEqualMatrix(summary2D, summary3DRows, summary3DCols, 0);
	}
    *condensateRows = summary3DCols;
    *condensateCols = summary3DDims + 1;
	return condensateSummary;
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
uint32_t** locate_mbus(int32_t** content, int32_t contentRows, int32_t contentCols, int32_t** pattern, int nRoundsInPattern, int datawidth, int* newRows, int* newCols){
	int i, j, z, nDetected = 0;
	uint32_t** summary = calloc((3 * contentRows), sizeof(uint32_t*));

	for (i = 0; i < (3 * contentRows); i++) {
		summary[i] = calloc(4, sizeof(uint32_t));
	}

	for (i = 0; i < nRoundsInPattern; i++) {
		for (j = 0; j < contentRows; j++) {
			if (content[j][3 * i] == 0xFFFFFFFF) {
				break;
			}
			else {
                int32_t corruptedBitsLength = 0;

				int32_t* corruptedBits = flipped_bits(content[j][3 * i + 1], pattern[i][1], datawidth, &corruptedBitsLength);
				if (corruptedBitsLength > 1){
                    summary[nDetected][0] = i;
                    summary[nDetected][1] = j;
                    summary[nDetected][2] = content[j][3*i];
                    summary[nDetected][3] = corruptedBitsLength;
                    nDetected++;
				}
			}
		}
	}
	/*int a,b;
    for (a = 0; a < 3 * contentRows; a++) {
        for (b = 0; b < 4; b++) {
            printf("%d ", summary[a][b]);
        }
    }*/
    return cutZerosFromArray(summary, (3*contentRows), 4, newRows, newCols, NULL, NULL);
}


// Esta función depende de la operación
int32_t** extractAnomalDVSelfConsistency(char* op, int32_t** opDVtotalrepetitions, int32_t opDVtotalrepetitionsRows, uint32_t** totalDVhistogram, int32_t histogramLenght, long int LN, int* nAddressesInRound, int nRoundsInPattern, int32_t** opdvhistogram, int32_t*** opdvmatrixbackup, int opdvmatrixbackupRows, int* XORANOMALS, int* lenght){
	int32_t opNthreshold;
	int* n_anomalous_values = malloc(sizeof(int));
    int testOPaLenght = 0;
	int32_t** testOPa;
	int32_t** oPExtracted_values = NULL;
	bool oPSelfConsistence = true;
	int n_anomalous_repetitions = 1, i, j, test, newLenght, propMCUSrows = 0, propMCUScols = 0, aux = 1;
	printf("\n\tDetermining the threshold for repetition excess:\n");
    
    if (strcmp(op, "xor") == 0) {
		printf("\n\t\tXOR operation...");
		opNthreshold = excessiveRepetitions(opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, LN, "xor", randomnessThreshold);
        // Redondeo hacia arriba posible solución?
        testOPa = find_anomalies_histogram(totalDVhistogram, histogramLenght, 1, opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, opNthreshold, &testOPaLenght, n_anomalous_values);
        *XORANOMALS = *n_anomalous_values;

		printf("\n\tXOR operation:\n");
		oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &newLenght, n_anomalous_values, 1);

		while (oPSelfConsistence) {
			printf("\t\tStep %d ", n_anomalous_repetitions);
            int nOPextrValues;
			oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &nOPextrValues, n_anomalous_values, n_anomalous_repetitions);
            
			int32_t** opPartRepsSelfCons = calloc(nOPextrValues, sizeof(int32_t*));
			for (i = 0; i < nOPextrValues; i++) {
				opPartRepsSelfCons[i] = calloc(nRoundsInPattern, sizeof(int32_t));
			}
			for (i = 0; i < nOPextrValues; i++) {
                for (j = 0; j < nRoundsInPattern; j++) {
					opPartRepsSelfCons[i][j] = opdvhistogram[(oPExtracted_values[i][0])-1][j+1];
				}
			}
            
			for (test = 0; test < nRoundsInPattern; test++) {
				printf("Test %d", test+1);
				int32_t** opdvmatrix = calloc(nAddressesInRound[test], sizeof(int32_t*));
				for (i = 0; i < nAddressesInRound[test]; i++) {
					opdvmatrix[i] = calloc(nAddressesInRound[test], sizeof(int32_t));
				}
                copyOfMatrix3to2(opdvmatrixbackup, opdvmatrix, nAddressesInRound[test], nAddressesInRound[test], test);
				bool** opMarkedPairs = marking_addresses(opdvmatrix, nAddressesInRound[test], nAddressesInRound[test], oPExtracted_values, nOPextrValues);
                
				int32_t** opProposedMcus = agrupate_mcus(opMarkedPairs, nAddressesInRound[test], nAddressesInRound[test], &propMCUSrows, &propMCUScols);
                int largestMCUSize = propMCUScols;
                
                int continuation = findEqual(opPartRepsSelfCons, nOPextrValues, nRoundsInPattern, 0, test, largestMCUSize);
                
                if (continuation != 0) {
                    oPSelfConsistence = false;
                }
                
			}
			if (oPSelfConsistence) {
				n_anomalous_repetitions++;
				printf("\n");
				if (n_anomalous_repetitions > *n_anomalous_values) {
					oPSelfConsistence = false;
					printf("\n\t\tNo more anomalous elements to check. Exiting\n");
				}
            }
            else {
                printf("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n");
                n_anomalous_repetitions--;
                oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &newLenght, n_anomalous_values, n_anomalous_repetitions);
                    
            }
        }

	} else if (strcmp(op, "pos") == 0){
		printf("\n\t\tPOS operation (It can take a long if there are too many addresses)...\n");
		opNthreshold = excessiveRepetitions(opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, LN, "pos", randomnessThreshold);
        testOPa = find_anomalies_histogram(totalDVhistogram, histogramLenght, 2, opDVtotalrepetitions, opDVtotalrepetitionsRows, 1, opNthreshold, &testOPaLenght, n_anomalous_values);

		printf("\n\tPOS operation:\n");
        oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &newLenght, n_anomalous_values, 1);
    
        while (oPSelfConsistence) {
            int nOPextrValues;
			printf("\t\tStep %d ", n_anomalous_repetitions);
            oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &nOPextrValues, n_anomalous_values, n_anomalous_repetitions);
            
            int32_t** opPartRepsSelfCons = calloc(nOPextrValues, sizeof(int32_t*));
            for (i = 0; i < nOPextrValues; i++) {
                opPartRepsSelfCons[i] = calloc(nRoundsInPattern, sizeof(int32_t));
            }
            for (i = 0; i < nOPextrValues; i++) {
                for (j = 0; j < nRoundsInPattern; j++) {
                    opPartRepsSelfCons[i][j] = opdvhistogram[(oPExtracted_values[i][0])-1][j+1];
                }
            }
            
			for (test = 0; test < nRoundsInPattern; test++) {
				printf("Test %d", test);

                int32_t** opdvmatrix = calloc(nAddressesInRound[test], sizeof(int32_t*));
                for (i = 0; i < nAddressesInRound[test]; i++) {
                    opdvmatrix[i] = calloc(nAddressesInRound[test], sizeof(int32_t));
                }
                copyOfMatrix3to2(opdvmatrixbackup, opdvmatrix, nAddressesInRound[test], nAddressesInRound[test], test);
                bool** opMarkedPairs = marking_addresses(opdvmatrix, nAddressesInRound[test], nAddressesInRound[test], oPExtracted_values, nOPextrValues);
                
                int32_t** opProposedMcus = agrupate_mcus(opMarkedPairs, nAddressesInRound[test], nAddressesInRound[test], &propMCUSrows, &propMCUScols);
                int largestMCUSize = propMCUScols;
                
                int continuation = findEqual(opPartRepsSelfCons, nOPextrValues, nRoundsInPattern, 0, test, largestMCUSize);
                if (continuation != 0) {
                    oPSelfConsistence = false;
                }
                
			}
			if (oPSelfConsistence) {
                n_anomalous_repetitions++;
                printf("\n");
                if (n_anomalous_repetitions > *XORANOMALS) {
                    oPSelfConsistence = false;
                    printf("\n\t\tNo more anomalous elements to check. Exiting\n");
                }
            }else {
                printf("\n\t\tViolation of self-consistence. Returning to previous state and exiting.\n");
                n_anomalous_repetitions--;
                oPExtracted_values = extract_some_critical_values(testOPa, testOPaLenght, &newLenght, n_anomalous_values, n_anomalous_repetitions);
            }
        }

	}
    *lenght = newLenght;
	return oPExtracted_values;
}

int32_t*** index2address(int32_t*** indexMatrix, int indexMatrixRows, int indexMatrixCols, int indexMatrixDims, int32_t** addressMatrix){
	int i, j, z;
	int32_t*** result = calloc(indexMatrixRows, sizeof(int32_t**));
	for (i = 0; i < indexMatrixRows; i++) {
		result[i] = calloc(indexMatrixCols, sizeof(int32_t*));
		for (j = 0; j < indexMatrixCols; j++) {
			result[i][j] = calloc(indexMatrixDims, sizeof(int32_t));
		}
	}
	for (z = 0; z < indexMatrixDims; z++) {
		for (i = 0; i < indexMatrixRows; i++) {
			for (j = 0; j < indexMatrixCols; j++) {
				if (indexMatrix[i][j][z] == 0) {
					break;
				}
				else {
					result[i][j][z] = addressMatrix[indexMatrix[i][j][z]][z];
				}
			}
		}
	}
	return result;
}

int32_t* criticalXORValuesFromClusters(int32_t*** PrelMCUSummary,int numRows, int numCols, int nRounds, int32_t** AddressMatrix, int32_t**XORextracted_values, int XORextractedRows, int XORextractedCols, int32_t** xorDVtotalrepetitions, int xorDvtotalrepetitionsLenght, int32_t** DVHistogram, long int LN, int* rows){

	int32_t* newCandidates = calloc(numRows*numCols*nRounds, sizeof(int32_t));
    int newCandidatesLength = 0;
	int ktest, kRow, kAddress1, kAddress2, i, j;
    int32_t index1, index2, address1, address2, candidate;
	for (ktest = 0; ktest < nRounds; ktest++){
		for (kRow = 0; kRow < numRows; kRow++){
			for (kAddress1 = 1; kAddress1 < numCols; kAddress1++){

				index1 = PrelMCUSummary[kRow][kAddress1][ktest];
				if (index1 == 0){
					break;
				}
				else{
					address1 = AddressMatrix[index1-1][ktest*3];
					for (kAddress2 = 0; kAddress2 < kAddress1; kAddress2++){
						index2 = PrelMCUSummary[kRow][kAddress2][ktest];
						address2 = AddressMatrix[index2-1][ktest*3];
						candidate = address1 ^ address2;
						if (findEqualMatrix(XORextracted_values, XORextractedRows, XORextractedCols, candidate) == 0){
                            int32_t* vCandidate = calloc(1, sizeof(int32_t));
                            vCandidate[0] = candidate;
                            newCandidates = unionVec(newCandidates, newCandidatesLength, vCandidate, 1, NULL, 0, &newCandidatesLength);
                            free(vCandidate);
						}
					}
				}		
			}
		}		
    }
    newCandidates = realloc(newCandidates, newCandidatesLength*sizeof(int32_t));

	int32_t xorNthreshold = excessiveRepetitions(xorDVtotalrepetitions, xorDvtotalrepetitionsLenght, 1, LN, "xor", randomnessThreshold);

    int32_t* purgedCandidatesXOR = calloc(newCandidatesLength, sizeof(int32_t));
    int purgedCandidatesLenght = 0;
    for (i = 0; i < newCandidatesLength; i++) {
        if (DVHistogram[newCandidates[i]-1][1] >= xorNthreshold) {
            purgedCandidatesXOR[purgedCandidatesLenght] = newCandidates[i];
            purgedCandidatesLenght++;
        }
    }
    purgedCandidatesXOR = realloc(purgedCandidatesXOR, purgedCandidatesLenght*sizeof(int32_t));
	*rows = purgedCandidatesLenght;

	return purgedCandidatesXOR;
}

void extractAnomalDVfromClusters(int32_t** content, int contentRows, int contentCols, int32_t** XORextracted_values, int XOREVRows, int XOREVCols, int32_t** POSextracted_values, int POSEVRows, int POSEVCols, int32_t** xorDVtotalrepetitions, int xorDVtotalrepetitionsLenght, int32_t** posDVtotalrepetitions, uint32_t** totalDVhistogram, int32_t*** xordvmatrixbackup, int xordvmatrixbackupRows, int xordvmatrixbackupCols, int xordvmatrixbackupDims, int32_t*** posdvmatrixbackup, int* nAddressesInRound, long int LN, int32_t** discoveredXORDVs, int* discoveredXORDVsLength, int32_t** discoveredPOSDVs, int* discoveredPOSDVsLength, int32_t*** tempXORDVvalues, int* tempXORDVvaluesLength, int32_t*** tempPOSDVvalues, int* tempPOSDVvaluesLength){

    // Cuatro últimos elementos de la función son las devoluciones, matrices por referecia?
    bool foundNewDVValues = true;
    int step = 0, i, j;
    
    *tempXORDVvalues = calloc(XOREVRows, sizeof(int32_t*));
     for (i = 0; i < XOREVRows; i++) {
         (*tempXORDVvalues)[i] = calloc(XOREVCols, sizeof(int32_t));
     }

     *tempPOSDVvalues = calloc(POSEVRows, sizeof(int32_t*));
     for (i = 0; i < POSEVRows; i++) {
         (*tempPOSDVvalues)[i] = calloc(POSEVCols, sizeof(int32_t));
     }

    //COPIA OBLIGADA, Se podría hacer función de copia de matriz simple a otra matriz simple
    for (i = 0; i < XOREVRows; i++) {
        for (j = 0; j < XOREVCols; j++) {
            (*tempXORDVvalues)[i][j] = XORextracted_values[i][j];
        }
    }
    for (i = 0; i < POSEVRows; i++) {
        for (j = 0; j < POSEVCols; j++) {
            (*tempPOSDVvalues)[i][j] = POSextracted_values[i][j];
        }
    }
    *tempXORDVvaluesLength = XOREVRows;
    *tempPOSDVvaluesLength = POSEVRows;
    
    while (foundNewDVValues) {
        step++;
        printf("\nStep: %d", step);
        printf("Creating preliminary organization...");
        int32_t*** tempCMB_summary;
        int cmbRows, cmbCols, cmbDims;
        int32_t*** tempXOR_summary;
        int xorRows, xorCols, xorDims;
        int32_t*** tempPOS_summary;
        int posRows, posCols, posDims;
        tempXOR_summary = propose_MCUs("xor", xordvmatrixbackup, xordvmatrixbackupRows, xordvmatrixbackupCols, xordvmatrixbackupDims, posdvmatrixbackup, *tempXORDVvalues, XOREVRows, *tempPOSDVvalues, POSEVRows, nAddressesInRound, &xorRows, &xorCols, &xorDims);
        tempPOS_summary = propose_MCUs("pos", xordvmatrixbackup, xordvmatrixbackupRows, xordvmatrixbackupCols, xordvmatrixbackupDims, posdvmatrixbackup, *tempXORDVvalues, XOREVRows, *tempPOSDVvalues, POSEVRows, nAddressesInRound, &posRows, &posCols, &posDims);
        tempCMB_summary = propose_MCUs("cmb", xordvmatrixbackup, xordvmatrixbackupRows, xordvmatrixbackupCols, xordvmatrixbackupDims, posdvmatrixbackup, *tempXORDVvalues, XOREVRows, *tempPOSDVvalues, POSEVRows, nAddressesInRound, &cmbRows, &cmbCols, &cmbDims);

        printf(" Ended. Searching new Elements...");

        printf(" XOR...");

        int nProposedXORDV = 0;
        int32_t* proposedXORDV = criticalXORValuesFromClusters(tempCMB_summary, cmbRows, cmbCols, cmbDims, content, *tempXORDVvalues, XOREVRows, XOREVCols, xorDVtotalrepetitions, xorDVtotalrepetitionsLenght, totalDVhistogram, LN, &nProposedXORDV);
		int32_t* proposedPOSDV = NULL;
        int nProposedPOSDV = 0;
        printf(" POS. SUB...");

        printf(" Ended. ");
        
        if (nProposedXORDV != 0) {
            printf("\n\t\tXOR operation:  %d new DV element", nProposedXORDV);
            if (nProposedXORDV != 1) {
                printf("s");
            }
            printf("found.");
            int nDiscoveredXORDVs = 0;
            (*discoveredXORDVs) = unionVec(proposedXORDV, nProposedXORDV, *discoveredXORDVs, 0, NULL, 0, &nDiscoveredXORDVs);
            *discoveredXORDVsLength = nDiscoveredXORDVs;

            int32_t** auxXORDV = calloc(XOREVRows + nProposedXORDV, sizeof(int32_t*));

            for (i = 0; i < XOREVRows + nProposedXORDV; i++) {
                auxXORDV[i] = calloc(XOREVCols, sizeof(int32_t));
            }
            for (i = 0; i < XOREVRows; i++) {
                for (j = 0; j < XOREVCols; j++) {
                    auxXORDV[i][j] = (*tempXORDVvalues)[i][j];
                }
            }
            int index = XOREVRows;
            for (int a = 0; a < nProposedXORDV; a++) {
                auxXORDV[index][0] = proposedXORDV[a];
                auxXORDV[index][1] = totalDVhistogram[proposedXORDV[a]-1][1];
                index++;
            }

            *tempXORDVvalues = auxXORDV;

            *tempXORDVvaluesLength = XOREVRows + nProposedXORDV;
            XOREVRows = XOREVRows + nProposedXORDV;

        }
        if (nProposedPOSDV != 0) {
            printf("\n\t\tPOS operation:  %d new DV element", nProposedPOSDV);
            if (nProposedPOSDV != 1) {
                printf("s");
            }
            printf("found.");
            int nDiscoveredPOSDVs = 0;
            (*discoveredPOSDVs) = unionVec(proposedPOSDV, nProposedPOSDV, (*discoveredPOSDVs), 0, NULL, 0, &nDiscoveredPOSDVs);
            *discoveredPOSDVsLength = nDiscoveredPOSDVs;

            int32_t** auxPOSDV = calloc(POSEVRows + nProposedPOSDV, sizeof(int32_t*));

            for (i = 0; i < POSEVRows + nProposedPOSDV; i++) {
                auxPOSDV[i] = calloc(POSEVCols, sizeof(int32_t));
            }
            for (i = 0; i < POSEVRows; i++) {
                for (j = 0; j < POSEVCols; j++) {
                    auxPOSDV[i][j] = (*tempPOSDVvalues)[i][j];
                }
            }
            int index = POSEVRows;
            for (int a = 0; a < nProposedPOSDV; a++) {
                auxPOSDV[index][0] = proposedPOSDV[a];
                auxPOSDV[index][1] = totalDVhistogram[proposedPOSDV[a]-1][2];
                index++;
            }
            (*tempPOSDVvalues) = auxPOSDV;
            *tempPOSDVvaluesLength = POSEVRows + nProposedPOSDV;
            POSEVRows = POSEVRows + nProposedPOSDV;
            
        }
        if ((nProposedXORDV == 0) && (nProposedPOSDV == 0)) {
            foundNewDVValues = false;
            printf("\n\n\t\tNO MORE ELEMENTS DISCOVERED. Stopping analysis.");
        }

    }
}


void criticalXORvaluesFromXORingRule(int32_t** XORextracted_values, int XORextracted_valuesRows, int32_t** XORDVtotalrepetitions, int xorDvtotalrepetitionsLength, uint32_t** totalDVhistogram, long int LN, int32_t*** newXORDVvalues, int* newXORDVvaluesLength, int32_t** candidates, int** candidatesLength){
	int i, j;
	bool candidateFind = false, prelCandidateFind = false;
	printf("\n\tCase 1: XORing known values... ");
	int32_t* prelCandidates = calloc(1, sizeof(int32_t));

	int32_t* oldXORValues = calloc(XORextracted_valuesRows, sizeof(int32_t));
    for (i = 0; i < XORextracted_valuesRows; i++) {
        oldXORValues[i] = XORextracted_values[i][0];
    }

    sort(oldXORValues, XORextracted_valuesRows);

	int k1, k2, prelCandidatesTam = 0,pos = 0,posCandidatesCase1 = 0;
	for (k1 = 0; k1 < XORextracted_valuesRows; k1++){
		int old1 = oldXORValues[k1];
		for (k2 = k1 + 1; k2 < XORextracted_valuesRows; k2++){
			candidateFind = false, prelCandidateFind = false;
			int old2 = oldXORValues[k2];
			int candidate = old1^old2;
            if ((findEqualVector(oldXORValues, XORextracted_valuesRows, candidate) == 0) && (findEqualVector(prelCandidates, prelCandidatesTam, candidate) == 0)) {
                int32_t* vCandidate = calloc(1, sizeof(int32_t));
                vCandidate[0] = candidate;
                prelCandidates = unionVec(prelCandidates, prelCandidatesTam, vCandidate, 1, NULL, 0, &prelCandidatesTam);
            }
		}

	}
	sort(prelCandidates, prelCandidatesTam);

	int nThresholdCase1 = excessiveRepetitions(XORDVtotalrepetitions, xorDvtotalrepetitionsLength, 1, LN, "xor", randomnessThreshold);

    int32_t* candidatesCase1 = calloc(prelCandidatesTam, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int candidatesCase1Length = 0;
    for (i = 0; i < prelCandidatesTam; i++) {
        if (totalDVhistogram[prelCandidates[i]-1][1] >= nThresholdCase1) {
            candidatesCase1[candidatesCase1Length] = prelCandidates[i];
            candidatesCase1Length++;
        }
    }
    candidatesCase1 = realloc(candidatesCase1, candidatesCase1Length*sizeof(int32_t));

	if (candidatesCase1Length > 0){
		printf("%d new value", candidatesCase1Length);
		if (candidatesCase1Length != 1)
			printf("s");
		printf(" obtained. \n");
	}
	else{
		printf("No new values. Going on.\n");
		}
	int newTamOldXORValues = 0;
	oldXORValues = unionVec(oldXORValues, XORextracted_valuesRows, candidatesCase1, candidatesCase1Length, NULL, 0, &newTamOldXORValues);
	sort(oldXORValues, newTamOldXORValues);

/* Case 2: XORingDV elements with confirmed XORs to obtain another one.*/

		printf("\tCASE 2: XORing DV elements with confirmed values... ");
		int32_t* DVElements = calloc(LN, sizeof(int32_t));
		int index = 0;
		for (i = 0; i < LN; i++){
			if (totalDVhistogram[i][1] != 0){
				DVElements[index] = i;
				index++;
			}
		}
    // This allows reobtaining the elements in the DV set.
		//bool* presencePrelCandidates = calloc(LN+1, sizeof(bool));

    // This +1 is placed since, sometimes, XORED can be 0.
		int32_t* selectedPrelCandidates = calloc(1, sizeof(int32_t));
        int selectedPrelCandidatesLength = 0;
		printf(" Searching... ");
		int kdv, kxor, cont = 0;
		bool xoredFind = false, kdvFind = false;
		for (kdv = 0; kdv < index; kdv++){
			for (kxor = 0; kxor < newTamOldXORValues; kxor++){
				int xored = DVElements[kdv] ^ oldXORValues[kxor];
				//presencePrelCandidates[xored + 1] = true;
                if ((findEqualVector(oldXORValues, newTamOldXORValues, xored) != 0) && (findEqualVector(oldXORValues, newTamOldXORValues, kdv) == 0)) {
                    int32_t* vKdv = calloc(1, sizeof(int32_t));
                    vKdv[0] = kdv;
                    selectedPrelCandidates = unionVec(selectedPrelCandidates, selectedPrelCandidatesLength, vKdv, 1, NULL, 0, &selectedPrelCandidatesLength);
                }
			}
		}

        // Get rid of possible ZEROS.
			//free(presencePrelCandidates);
        //releasing RAM.
        int nThresholdCase2 = excessiveRepetitions(XORDVtotalrepetitions, xorDvtotalrepetitionsLength, 1, LN, "xor", randomnessThreshold);

    int32_t* candidatesCase2 = calloc(selectedPrelCandidatesLength, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int candidatesCase2Length = 0;
    for (i = 0; i < selectedPrelCandidatesLength; i++) {
        if (totalDVhistogram[prelCandidates[i]-1][1] >= nThresholdCase2) {
            candidatesCase2[candidatesCase2Length] = selectedPrelCandidates[i];
            candidatesCase2Length++;
        }
    }
    candidatesCase2 = realloc(candidatesCase2, candidatesCase2Length*sizeof(int32_t));

    if (candidatesCase2Length > 0){
        printf("%d new value", candidatesCase2Length);
        if (candidatesCase2Length != 1)
            printf("s");
            printf(" obtained. \n");
    }else{
        printf("No new values. Going on.\n");
    }

    int tamCandidates = 0;
    candidates = unionVec(candidatesCase1, candidatesCase1Length, candidatesCase2, candidatesCase2Length, NULL, 0, &tamCandidates);
    *newXORDVvaluesLength = XORextracted_valuesRows+tamCandidates;
    (*newXORDVvalues) = calloc((XORextracted_valuesRows+tamCandidates), sizeof(int32_t*));
    for (i = 0; i < *newXORDVvaluesLength; i++) {
        (*newXORDVvalues)[i] = calloc(2, sizeof(int32_t));
    }
    index = 0;
    for (i = 0; i < XORextracted_valuesRows; i++) {
        for (j = 0; j < 2; j++) {
            (*newXORDVvalues)[i][j] = XORextracted_values[i][j];
        }
        index++;
    }
    for (i = 0; i < tamCandidates; i++) {
        (*newXORDVvalues)[index][0] = (*candidates)[i];
        (*newXORDVvalues)[index][1] = totalDVhistogram[(*candidates)[i]-1][1];
        index++;
    }

}


#endif
