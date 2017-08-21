#ifndef libs_h
#define libs_h

#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "strucs.h"

/**** CONSTANTS ****/

/* The following constant is extremely important: it determines the
* expected number values that separates randomness from causality. We set
* this value to 0.05 (1 chance out of 20) */
#define randomnessThreshold 0.05

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
 * - Modificar los libera de las matrices para adaptarlos a los structs
 * -
 */

/*
 * Function to free type int32_t 2D matrix's dynamic memory.
 */
void liberaInt(int32_t** matrix, int rows){
    int i;
    for (i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

/*
 * Function to free type uint32_t 2Dmatrix's dynamic memory.
 */
void liberaUint(uint32_t** matrix, int rows){
    int i;
    for (i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}

/*
 * Function to free type int32_t 3Dmatrix's dynamic memory.
 */
void libera3D(int32_t*** matrix, int row, int col){
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}

/*
 * Function to open a file to a concrete pointer.
 * Gets the pointer which points the file.
 * Returns false if the file is missing.
 */
FILE* openFile(char* name){
    FILE* file;
    file = fopen(name, "r");
    if (file == NULL) {
        printf("ERROR OPENING FILE.\n");
    }
    return file;
}

/*
 * This function allows getting characters from pointer's file, transforming them into fields.
 * Returns the complete field.
 */
char* getfield(FILE* file, char* taken){
	int c = 0, elems = 0;
	char* ch = malloc(sizeof(char));
	c = getc(file);
	while ((c != ',') && (c != '\n')) { // Important characters on csvs
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
 * Function which sorts a vector into ascending order.
 */
void sort(vectorInt32Struct vector){
	int min, i, j, aux;
	for (i = 0; i < vector.length - 1; i++){
		min = i;
		for (j = i + 1; j < vector.length; j++)
			if (vector.data[min] > vector.data[j])
				min = j;

		aux = vector.data[min];
		vector.data[min] = vector.data[i];
		vector.data[i] = aux;
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
	/* Comparing last position of the array with 0xFFFFFFFF == 4294967295 */
	if (columnVector[*elems - 1] == 0xFFFFFFFF){
		i = 0;
		uint32_t* newColumnVector;
		newColumnVector = (uint32_t*)malloc(sizeof(uint32_t)*(*elems));
		while (columnVector[i] != 4294967295){
			newColumnVector[i] = columnVector[i];
			i++;
		}
		*elems = i;
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
vectorInt32Struct flipped_bits(int32_t word, int32_t pattern, int wordWidth){
	int k, nFound = 0;
    vectorInt32Struct result;
    result.data = calloc(wordWidth, sizeof(int32_t));
	for (k = 0; k < wordWidth; k++) {
		result.data[k] = 1;
		result.data[k] *= (wordWidth + 1);
	}

	int32_t bitflips = word ^ pattern;
	for (k = 0; k < wordWidth; k++){
		if ((bitflips % 2) == 1) {
            result.data[nFound] = k;
            nFound++;
        }
		bitflips = bitflips >> 1;
	}

    if (nFound > 0) {
        result.data = realloc(result.data, nFound*sizeof(int32_t));
    }
    result.length = nFound;
	return result;
}

/*
 * Function to transpose a matrix.
 */
matrixUint322DStruct transposed_matrix(matrixUint322DStruct matrix){
    matrixUint322DStruct transposedMatrix;
    transposedMatrix.rows = matrix.rows;
    transposedMatrix.cols = matrix.cols;
    transposedMatrix.data = calloc(matrix.rows, sizeof(uint32_t*));
	for (int i = 0; i < transposedMatrix.rows; i++)
		transposedMatrix.data[i] = calloc(transposedMatrix.cols, sizeof(uint32_t));

	for (int i = 0; i < matrix.rows; i++){
		for (int j = 0; j < matrix.cols; j++){
			transposedMatrix.data[i][j] = matrix.data[j][i];
		}
	}
	return transposedMatrix;
}

//Funcion de relleno de fila superior y columna izquierda vector, elemens y matriz para modificarla


/*
 * Function to create the DV sets associated with the XOR and positive operations.
 * The idea is simple: The address vector is replicated to create a square matrix
 * traspose it, and operating. Every element in elements above the main
 * diagonal is an element of the DV.
 */
matrixUint322DStruct create_DVmatrix(uint32_t* addresses, int elems, char* op, uint32_t* RWcycles, int nbits4blocks, long int ln){
    int i, j;
    matrixUint322DStruct DVmatrix;
    DVmatrix.rows = elems;
    DVmatrix.cols = elems;
    DVmatrix.data = calloc(DVmatrix.rows, sizeof(uint32_t *));
    for (i = 0; i < DVmatrix.rows; i++){
        DVmatrix.data[i] = calloc(DVmatrix.cols, sizeof(uint32_t));
    }
    
    /* DVmatrixOP is a matrix to operate */
    matrixUint322DStruct DVmatrixOP;
    DVmatrixOP.rows = elems;
    DVmatrixOP.cols = elems;
    DVmatrixOP.data = calloc(DVmatrixOP.rows, sizeof(uint32_t *));
    for (i = 0; i < DVmatrixOP.rows; i++){
        DVmatrixOP.data[i] = calloc(DVmatrixOP.cols, sizeof(uint32_t));
    }
    for (i = 0; i < DVmatrixOP.rows; i++) {
        for (j = 0; j < DVmatrixOP.cols; j++) {
            DVmatrixOP.data[i][j] = addresses[i];
        }
    }
    
    /* DVmatrixTrans is the transposed matrix to operate */
    matrixUint322DStruct DVmatrixTrans = transposed_matrix(DVmatrixOP);

	if (strcmp(op, "xor") == 0){
		for (i = 0; i < DVmatrixOP.rows; i++){
            for (j = 0; j < DVmatrixOP.cols; j++){
				DVmatrix.data[i][j] = DVmatrixOP.data[i][j] ^ DVmatrixTrans.data[i][j];
			}
        }
	}
	else if (strcmp(op, "pos") == 0){
		for (i = 0; i < DVmatrixOP.rows; i++){
            for (j = 0; j < DVmatrixOP.cols; j++){
                 DVmatrix.data[i][j] = abs(DVmatrixOP.data[i][j] - DVmatrixTrans.data[i][j]);
			}
		}
	}
	else {
		printf("\n A problem creating the DV set. Operation not recognized\n");
		EXIT_FAILURE;
	}
    
    /* The elements of the lower triangular must be deleted. */
    for (i = 0; i < DVmatrix.rows; i++) {
        for (j = 0; j < DVmatrix.cols; j++) {
            if ((i == j) || (j < i)) {
                DVmatrix.data[i][j] = 0;
            }
        }
    }

	/* Now, we will incorporate the information about the R-W cycles in the same
	 * way. We create a matrix with the vector, traspose it and look for coincidence
	 * element to element. Later, the matrix with True/False elements is converted
	 * into integer and multiplied element to element with the DVmatrix. */
    for (i = 0; i < DVmatrix.rows; i++) {
        for (j = 0; j < DVmatrix.cols; j++) {
            DVmatrixOP.data[i][j] = RWcycles[i];
        }
    }
    
    DVmatrixTrans = transposed_matrix(DVmatrixOP);
    /* coincidence is a matrix to discover the elements which are in both operation matrices. */
    matrixUint322DStruct coincidence;
    coincidence.rows = elems;
    coincidence.cols = elems;
    coincidence.data = calloc(coincidence.rows, sizeof(uint32_t *));
    for (i = 0; i < elems; i++){
        coincidence.data[i] = calloc(coincidence.cols, sizeof(uint32_t));
    }
    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            if (DVmatrixOP.data[i][j] == DVmatrixTrans.data[i][j]) {
                coincidence.data[i][j] = 1;
            } else {
                coincidence.data[i][j] = 0;
            }
        }
    }
    
    /* The elements of the lower triangular must be deleted. */
    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            if ((i == j) || (j < i)) {
                coincidence.data[i][j] = 0;
            }
        }
    }

    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            DVmatrix.data[i][j] = DVmatrix.data[i][j] * coincidence.data[i][j];
        }
    }

	/* The same for the block division. We multiply the address vector by 2^bits4blocks/LN
	 * and round to the floor integer. Traspose and multiply. */
	int32_t* addressblocks = calloc(elems, sizeof(int32_t));
	for (i = 0; i < coincidence.rows; i++) {
		addressblocks[i] = floor(addresses[i] / ln);
	}
    
    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            DVmatrixOP.data[i][j] = addressblocks[i];
        }
    }
    
    DVmatrixTrans = transposed_matrix(DVmatrixOP);
    
    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            if (DVmatrixOP.data[i][j] == DVmatrixTrans.data[i][j]) {
                coincidence.data[i][j] = 1;
            } else {
                coincidence.data[i][j] = 0;
            }
        }
    }
    
    /* The elements of the lower triangular must be deleted */
    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            if ((i == j) || (j < i)) {
                DVmatrix.data[i][j] = 0;
            }
        }
    }

    for (i = 0; i < coincidence.rows; i++) {
        for (j = 0; j < coincidence.cols; j++) {
            DVmatrix.data[i][j] = DVmatrix.data[i][j] * coincidence.data[i][j];
        }
    }

    liberaUint(coincidence.data, elems);
    liberaUint(DVmatrixOP.data, elems);
    liberaUint(DVmatrixTrans.data, elems);
	return DVmatrix;
}

/*
 * This function creates the DVvector cotaining the elements of the triu in
 * the DVmatrix. elems is the vector's size.
 * Returns that vector.
 */
vectorUint32Struct create_DVvector(matrixUint322DStruct* DVmatrix){
	int ndv = 0.5*DVmatrix->rows*(DVmatrix->rows - 1);
    // Vector with the elements on the upper triangular of the selected matrix
    vectorUint32Struct DVvector;
    DVvector.data = calloc(ndv, sizeof(uint32_t));

	/* The elements of the matrix are recovered. Elements are added fixing the
	* column and going down. */
	int index = 0;
    
    for (int kcol = 1; kcol < DVmatrix->rows; kcol++) {
        for (int krow = 0; krow < kcol; krow++) {
            if (DVmatrix->data[krow][kcol] != 0) {
                DVvector.data[index] = DVmatrix->data[krow][kcol];
                index++;
            }
        }
    }
    DVvector.data = realloc(DVvector.data, index*sizeof(uint32_t));
    DVvector.length = index;
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

/*
 * This function transforms a 3D matrix on its selected dimension on a 2D matrix
 */
matrixUint322DStruct copyOfMatrix3UintTo2(matrixUint323DStruct* matrix, matrixUint322DStruct* matrixbackup, int ktest){
	int i, j;
	for (i = 0; i < matrixbackup->rows; i++){
        for (j = 0; j < matrixbackup->cols; j++) {
            matrixbackup->data[i][j] = matrix->data[i][j][ktest];
        }
	}
    return *matrixbackup;
}

/*
 * This function transforms a 3D matrix on its selected dimension on a 2D matrix
 */
matrixUint322DStruct copyOfMatrix3to2(matrixInt323DStruct* matrix, matrixUint322DStruct* matrixbackup, int ktest){
    int i, j;
    for (i = 0; i < matrixbackup->rows; i++){
        for (j = 0; j < matrixbackup->cols; j++) {
            matrixbackup->data[i][j] = matrix->data[i][j][ktest];
        }
    }
    return *matrixbackup;
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
 * in the correct order
 */
int32_t countsOfElems(matrixInt322DStruct histogram, int maxValue, int col){
	int i;
    int32_t repetitionstmp = 0;
	for (i = 0; i < histogram.rows; i++) {
		if (histogram.data[i][col] > 0 && histogram.data[i][col] < maxValue) {
			repetitionstmp++;
		}
	}
	return repetitionstmp;
}

/*
 *
 */
vectorInt32Struct countsInt(matrixInt322DStruct* matrix, int col, int maxValue){
    int i, j, count = 0;
    vectorInt32Struct result;
    result.length = maxValue + 1;
    result.data = calloc(result.length, sizeof(int32_t));
    for (i = 0; i < result.length; i++) {
        for (j = 0; j < matrix->rows; j++) {
            if (matrix->data[j][col] == i) {
                count++;
            }
        }
        result.data[i] = count;
        count = 0;
    }
    return result;
}

/*
 *
 */
vectorInt32Struct countsUint(matrixUint322DStruct* matrix, int col, int maxValue){
    int i, j, count = 0;
    vectorInt32Struct result;
    result.length = maxValue + 1;
    result.data = calloc(result.length, sizeof(int32_t));
    for (i = 0; i < result.length; i++) {
        for (j = 0; j < matrix->rows; j++) {
            if (matrix->data[j][col] == i) {
                count++;
            }
        }
        result.data[i] = count;
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
int lookForElem(vectorUint32Struct* vector, int position){
	int sum = 0, i;
	for (i = 0; i < vector->length; i++) {
		if (vector->data[position] == vector->data[i]) {
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
vectorInt32Struct create_histogram(vectorUint32Struct* vector, long int ln){
    vectorInt32Struct histogram;
    histogram.length = ln;
    histogram.data = calloc(ln, sizeof(int32_t));
	int k = 0, sum = 0;
	for (k = 0; k < vector->length; k++) {
		sum = lookForElem(vector, k);
		histogram.data[vector->data[k]] = sum;
        sum = 0;
	}
	return histogram;
}

/*
 *
 */
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

/*
 * Concatenates two matrices with the same number of columns
 */
matrixInt322DStruct vcat(matrixInt322DStruct* matrixA, matrixInt322DStruct* matrixB){
	int numRow = matrixA->rows + matrixB->rows;
    matrixInt322DStruct matrixC;
    matrixC.rows = numRow;
    matrixC.cols = matrixA->cols;
	matrixC.data = calloc(matrixC.rows, sizeof(int32_t *));
    for (int i = 0; i < matrixC.rows; i++){
		matrixC.data[i] = calloc(matrixA->cols, sizeof(int32_t));
    }
	for (int i = 0; i < matrixA->rows; i++){
		for (int j = 0; j < matrixA->cols; j++) {
			matrixC.data[i][j] = matrixA->data[i][j];
		}
	}
	for (int i = matrixA->rows; i < matrixC.rows; i++){
		for (int j = 0; j < matrixB->cols; j++){
			matrixC.data[i][j] = matrixB->data[i][j];
		}
	}
	return matrixC;
}

/* 
 * This function just implements the theoretical expresions to determine expected
 * repetitions. Variable names are meaningful.
 */
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
 * This function allows calculating the lowest number of repetitions from which
 * the expected number is below threshold.
 * RepInHistVector is a column vector containing the statistics.It is assumed
 * that the first element correspond to 0 repetitions.
 */
int32_t excessiveRepetitions(matrixInt322DStruct* repInHistVector, int repInHistVectorCol, long int LN, char* operation, double threshold){
	int32_t nthreshold = 0, ndv = 0, k, i;
    vectorInt32Struct occurrenceIndex;
    occurrenceIndex.length = repInHistVector->rows;
	occurrenceIndex.data = calloc(occurrenceIndex.length, sizeof(int32_t));
	for (i = 0; i < repInHistVector->rows; i++)
		occurrenceIndex.data[i] = i;

	for (i = 0; i < occurrenceIndex.length; i++)
		ndv += ((int32_t)repInHistVector->data[i][repInHistVectorCol] * occurrenceIndex.data[i]);

	k = 0;
	if (strcmp(operation, "pos") == 0){
        if ((excessiveRepetitions(repInHistVector, repInHistVectorCol, LN, "xor", threshold) - 1) > 0){
			k = excessiveRepetitions(repInHistVector, repInHistVectorCol, LN, "xor", threshold) - 1;
        } else {
            k = 0;
        }
    } else {
        k = 0;
    }

    for (i = k; i < occurrenceIndex.length; i++) {
        if (expectedRepetitions(i, LN, ndv, operation) < threshold) {
            nthreshold = i;
            break;
        }
    }
    free(occurrenceIndex.data);
	return nthreshold;
}

/*
 * This function lists the element anomalously repeated in the histogram.
 * occurrences is a column vector [V(0)...V(K)] where V(K) indicates the number
 * of elements that are repeated K times.
 * nthreshold indicates the limit value from which repetitions should not
 * be expected in only-SBU experiments.
 */
matrixInt322DStruct find_anomalies_histogram(matrixUint322DStruct *histogram, int histogramCol, matrixInt322DStruct *occurrences, int occurrencesCol, int32_t nthreshold, int* n_anomalous_values){ //n_anomalous_values por referencia porque necesitamos el valor
	/* The total number of elements anomalously repeated. We use +1 since the first
	* element of the vector is related to 0 repetitions.
	*/
	int i, index = 0, k, k2;
	int32_t nAnomalies = 0;

	for (i = nthreshold; i < occurrences->rows; i++) {
		nAnomalies += occurrences->data[i][occurrencesCol];
	}
	/*The number of repetitions of the most repeated element. Coincides roughly with the vector length. */
	int32_t highestAnomaly = occurrences->rows - 1;
    matrixInt322DStruct anomalies;
    anomalies.rows = nAnomalies;
    anomalies.cols = 2;
    anomalies.data = calloc(anomalies.rows, sizeof(int32_t*));
	for (i = 0; i < anomalies.rows; i++) {
		anomalies.data[i] = calloc(anomalies.cols, sizeof(int32_t));
	}
	/* the following variable will be used to indicate the number of elements different
	* from 0, above nthreshold, in the occurrences vector. */
	*n_anomalous_values = 0;
    anomalies.rows = nAnomalies;
    for (k = highestAnomaly; k > nthreshold-1; k--) {
		/* I've prefered to implement this solution instead of using the native function
		* "find()" since I believe this solution is faster as it includes breaks once
		* the number of anomalous occurrences are achieved. */
		if (occurrences->data[k][occurrencesCol] != 0){
			int n_occurrence_value = 0;
			int occurrence_value = k;
			*n_anomalous_values += 1;
			for (k2 = 0; k2 < histogram->rows; k2++) {
				if (histogram->data[k2][histogramCol] == occurrence_value){
					anomalies.data[index][0] = k2 + 1;
					anomalies.data[index][1] = occurrence_value;
					index++;
					n_occurrence_value++;
					if (n_occurrence_value == occurrences->data[k][occurrencesCol]) {
						break;
					}
				}
			}
		}
	}
	return anomalies;
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
matrixInt322DStruct extract_some_critical_values(matrixInt322DStruct *anomalous_repetitions, int* n_anomalous_repetitions, int n){
	int observed_values = 1, nrow = 1, k;
	if (*n_anomalous_repetitions <= n) {
        return *anomalous_repetitions;
	} else {
		for (k = 1; k < anomalous_repetitions->rows; k++){
			if (anomalous_repetitions->data[k][1] != anomalous_repetitions->data[k - 1][1]) {
				if (observed_values == n) {
					nrow = k;
				}
				observed_values++;
			}
		}
        matrixInt322DStruct extracted_values;
        extracted_values.rows = nrow;
        extracted_values.cols = 2;
		extracted_values.data = calloc(extracted_values.rows, sizeof(int32_t*));
		for (k = 0; k < extracted_values.rows; k++) {
			extracted_values.data[k] = calloc(extracted_values.cols, sizeof(int32_t));
		}
		for (k = 0; k < extracted_values.rows; k++) {
			extracted_values.data[k][0] = anomalous_repetitions->data[k][0];
			extracted_values.data[k][1] = anomalous_repetitions->data[k][1];
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
matrixBool2DStruct marking_addresses(matrixInt322DStruct* davmatrix, matrixInt322DStruct* critical_values){
	int i, x, y;
	int32_t value;
    matrixBool2DStruct result;
    result.rows = davmatrix->rows;
    result.cols = davmatrix->cols;
    result.data = calloc(result.rows, sizeof(bool*));
	for (i = 0; i < result.rows; i++) {
		result.data[i] = calloc(result.cols, sizeof(bool));
	}

	for (i = 0; i < critical_values->rows; i++) {
		value = critical_values->data[i][0];
		for (x = 0; x < result.rows; x++) { /* High cost because its necessary to go through the entire matrix for each element */
			for (y = 0; y < result.cols; y++) {
				if (davmatrix->data[x][y] == value) {
					result.data[x][y] = true;
				}
			}
		}
	}
	return result;
}

/*
 * Esta función, crea un vector temporal de tamaño thesummaryCols, para hacer la unión entre las filas addressRowMin y
 * addressRowMax, colocando los elementos en orden mientras se insertan.
 * This function makes a temp vector with the same dimensions as thesummary, to do the union between the rows
 * addressRowMin and addressRowMax, placing the elements in order during the insertion
 */
vectorInt32Struct unionAgrupate_mcus(matrixInt322DStruct* thesummary, int32_t addressRowMin, int32_t addressRowMax){
	int i, j, x;
    vectorInt32Struct thesummarytemp;
    thesummarytemp.length = thesummary->cols;
	thesummarytemp.data = calloc(thesummarytemp.length, sizeof(int32_t));
	// We need a vector to save the two rows' elements for the union
	j = 0; // thesummarytemp index
	for (i = 0; i < thesummary->cols; i++) {
		if (thesummary->data[addressRowMin][i] != 0) {
			if (thesummarytemp.data[0] == 0) { // Empty
				thesummarytemp.data[j] = thesummary->data[addressRowMin][i];
				j++;
			}
			else if (thesummarytemp.data[j - 1] < thesummary->data[addressRowMin][i]){ // Last value lowest
				thesummarytemp.data[j] = thesummary->data[addressRowMin][i];
				j++;
			}
			else if (thesummarytemp.data[j - 1] > thesummary->data[addressRowMin][i]){ // Last value highest, need a gap
				x = j - 1;
				while (thesummarytemp.data[x] > thesummary->data[addressRowMin][i]) {
					thesummarytemp.data[x + 1] = thesummarytemp.data[x]; // One scrolls
					x--;
				}
				thesummarytemp.data[x + 1] = thesummary->data[addressRowMin][i]; // The value is copied in the correspondent gap
			}
		}
	}
    for (i = 0; i < thesummary->cols; i++) {
		if (thesummary->data[addressRowMax][i] != 0) {
			if (thesummarytemp.data[0] == 0) { // Empty
				thesummarytemp.data[j] = thesummary->data[addressRowMax][i];
				j++;
			}
			else if (thesummarytemp.data[j - 1] < thesummary->data[addressRowMax][i]){ // Last value lowest
				thesummarytemp.data[j] = thesummary->data[addressRowMax][i];
				j++;
			}
			else if (thesummarytemp.data[j - 1] > thesummary->data[addressRowMax][i]){ // Last value highest, need a gap
				x = j - 1;
				while (thesummarytemp.data[x] > thesummary->data[addressRowMax][i]) {
					thesummarytemp.data[x + 1] = thesummarytemp.data[x]; // One scrolls
					x--;
				}
				thesummarytemp.data[x + 1] = thesummary->data[addressRowMax][i]; // The value is copied in the correspondent gap
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

matrixInt322DStruct agrupate_mcus(matrixBool2DStruct* davmatrix){
    /* Let us locate the elements */
    /* Find the indexes with nonzero values */
    int count = 0, i, j, k, x, xx;
    vectorInt32Struct nonzerovectorelements;
    nonzerovectorelements.length = davmatrix->rows * davmatrix->cols;
    nonzerovectorelements.data = calloc(davmatrix->rows*davmatrix->cols, sizeof(int32_t));
    for (i = 0; i < davmatrix->rows; i++) {
        for (j = 0; j < davmatrix->cols; j++) {
            if (davmatrix->data[i][j] != false) {
                nonzerovectorelements.data[count] = j * davmatrix->cols + i;
                count++; /* Number of elements known but they can be repeated */
            }
        }
    }
    nonzerovectorelements.length = count;
    nonzerovectorelements.data = realloc(nonzerovectorelements.data, nonzerovectorelements.length * sizeof(int32_t));
    
    /* Now, transform the elements in pairs creating a matrix with 2 columns. */
    matrixInt322DStruct relatedpairs;
    relatedpairs.rows = nonzerovectorelements.length;
    relatedpairs.cols = 2;
    relatedpairs.data = calloc(relatedpairs.rows, sizeof(int32_t*));
    for (i = 0; i < relatedpairs.rows; i++) {
        relatedpairs.data[i] = calloc(relatedpairs.cols, sizeof(int32_t));
    }
    
    /*    @simd for k1 = 1:length(nonzerovectorelements)
     # The first operation calculates the column
     @inbounds relatedpairs[k1,:] = [rem(nonzerovectorelements[k1], NRows) div(nonzerovectorelements[k1],NRows)+1]
     # REM means "remainder of division", DIV integer division.
     end*/
    for (i = 0; i < relatedpairs.rows; i++) {
        if((nonzerovectorelements.data[i] % davmatrix->rows != 0)){
        relatedpairs.data[i][0] = nonzerovectorelements.data[i] % davmatrix->rows+1;
        }else{
            relatedpairs.data[i][0] = nonzerovectorelements.data[i] % davmatrix->rows;
        }
        relatedpairs.data[i][1] = (nonzerovectorelements.data[i] / davmatrix->rows)+1;
    }
    
    free(nonzerovectorelements.data);
    /* Good, we have created a matrix containing the related pairs. Now, we must
     # group them in larger events.
     
     # First step: A matrix initilizing the events. Perhaps too large, but it is better
     # to have spare space.
     */
    matrixInt322DStruct thesummary;
    thesummary.rows = count;
    thesummary.data = calloc(count, sizeof(int32_t*));
    if (count < NMaxAddressesInEvent) { /* If the minimun is newCount */
        thesummary.cols = count;
        for (i = 0; i < count; i++) {
            thesummary.data[i] = calloc(thesummary.cols, sizeof(int32_t));
        }
    } else { /* If the minimun is NMaxAddressesInEvent */
        thesummary.cols = NMaxAddressesInEvent;
        for (i = 0; i < count; i++) {
            thesummary.data[i] = calloc(thesummary.cols, sizeof(int32_t));
        }
    }
    
    /* Now, we save the first pair in the first column: */
    int totalEvents = 0, thelargestMCU = 0;
    int32_t firstAddress, secondAddress, firstAddressRow = 0, secondAddressRow = 0;
    for (k = 0; k < count; k++) {
        //There are several cases. Let us list them in order.
        //Positions are extrated from the pairs.
        firstAddress = relatedpairs.data[k][0];
        secondAddress = relatedpairs.data[k][1];
        //TENEMOS QUE MIRAR SI firstAddress y secondAddress están en thesummary
        bool firstAddressFound = false, secondAddressFound = false;
        for (i = 0; i < thesummary.rows; i++) {
            for (j = 0; j < thesummary.cols; j++) {
                if ((thesummary.data[i][j] == firstAddress) && (firstAddressFound == false)) {
                    firstAddressFound = true;
                    firstAddressRow = i;
                }
                if ((thesummary.data[i][j] == secondAddress) && (secondAddressFound == false)) {
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
            thesummary.data[totalEvents][0] = relatedpairs.data[k][0];
            thesummary.data[totalEvents][1] = relatedpairs.data[k][1];
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
            while ((thesummary.data[firstAddressRow][x] != 0) && (x < thesummary.cols)) {
                x++;
            }
            //Fixed problem. We do know exactly the position to put the new value.
            if (x < thesummary.cols) {
                thesummary.data[firstAddressRow][x] = secondAddress;
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
            while ((thesummary.data[secondAddressRow][x] != 0) && (x < thesummary.cols)) {
                x++;
            }
            //Fixed problem. We do know exactly the position to put the new value.
            if (x < thesummary.cols) {
                thesummary.data[secondAddressRow][x] = firstAddress;
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
                vectorInt32Struct thesummarytemp;
                thesummarytemp = unionAgrupate_mcus(&thesummary, addressRowMin, addressRowMax);

                /*
                 * But this is a column vector where the elements are not repeated but
                 * not disposed in increasing order. Let us solve this:
                 */
                /*
                 * Vec vectorizes the elements of the matrix.
                 * Time to put the elements in the row.
                 */
                for (i = 0; i < thesummary.cols; i++) {
                    thesummary.data[addressRowMin][i] = thesummarytemp.data[i]; // Se copia el vector de la unión en la matriz
                    thesummary.data[addressRowMax][i] = 0;
                }
                //Also, as two events have been merged, the number of total events is lower:
                totalEvents--;
            }
        }
    }
    
    liberaInt(relatedpairs.data, count);
    /*
     #Really good. Now, we will separate the cells with information from those with
     #zeros.
     # The number of rows is easy to calculate: NTotalEvents. Concerning the other
     # element:
     */
    int ind = 0, sum = 0;
    for(j = 2; j < thesummary.cols; j++) {
        for (i = 0; i < thesummary.rows; i++) {
            if(thesummary.data[i][j] != 0){
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
    matrixInt322DStruct result;
    result.rows = totalEvents;
    result.cols = thelargestMCU;
    result.data = calloc(result.rows, sizeof(int32_t*));
    for (i = 0; i < result.rows; i++) {
        result.data[i] = calloc(result.cols, sizeof(int32_t));
    }
    for (x = 0; x < result.rows; x++) {
        for (xx = 0; xx < result.cols; xx++) {
            result.data[x][xx] = thesummary.data[x][xx];
        }
    }

    liberaInt(thesummary.data, thesummary.rows);
    return result;
}

/*
 * This funtion returns the number of elements which are different to the number searched
 */
int findNotEqualMatrix(matrixInt322DStruct* matrix, int search){
    int count = 0, i, j;
    for (i = 0; i < matrix->rows; i++) {
        for (j = 0; j < matrix->cols; j++) {
            if (matrix->data[i][j] != search) {
                count++;
            }
        }
    }
    return count;
}

/*
 * This function returns the number of elements which are equal to the number searched
 */
int findEqualMatrix(matrixInt322DStruct* matrix, int search){
    int count = 0, i, j;
    for (i = 0; i < matrix->rows; i++) {
        for (j = 0; j < matrix->cols; j++) {
            if (matrix->data[i][j] == search) {
                count++;
            }
        }
    }
    return count;
}

/*
 * Returns the number of elements which are equal to the number searched
 */
int findEqualVector(vectorInt32Struct* vector, int search){
    int count = 0, i;
    for (i = 0; i < vector->length; i++) {
        if (vector->data[i] == search) {
            count++;
        }
    }
    return count;
}

/*
 * This function returns the number of the elements which are diferent to the number searched on a selected
 * row or a selected column
 */
int findNotEqual(matrixInt322DStruct* matrix, int selectedRow, int selectedCol, int search){
    int count = 0, i;
    if (selectedCol == -1) { /* Searching on columns of the selected row */
        for (i = 0; i < matrix->cols; i++) {
            if (matrix->data[selectedRow][i] != search) {
                count++;
            }
        }
    } else { /* Searching on rows of the selected column */
        for (i = 0; i < matrix->rows; i++) {
            if (matrix->data[i][selectedCol] != search) {
                count++;
            }
        }
    }
    return count;
}

/*
 * This function returns the number of the elements which are equal to the number searched on a selected
 * row or a selected column
 */
int findEqual(matrixInt322DStruct* matrix, int selectedRow, int selectedCol, int search){
    int count = 0, i;
    if (selectedCol == -1) { /* Searching on columns of the selected row */
        for (i = 0; i < matrix->cols; i++) {
            if (matrix->data[selectedRow][i] == search) {
                count++;
            }
        }
    } else { /* Searching on rows of the selected column */
        for (i = 0; i < matrix->rows; i++) {
            if (matrix->data[i][selectedCol] == search) {
                count++;
            }
        }
    }
    return count;
}

/*
 * This function allows to cut files and rows of redundant zeros in 1 or 2
 * dimensions matrix.
 */
matrixInt322DStruct cutZerosFromArray(matrixInt322DStruct *matrix, bool* isMatrix, bool* isVector){
	int i, j = 0, countIni = 0, countFin = 0, newRows = 0, newCols = 0;
    
    if (isMatrix != NULL && isVector != NULL) {
        *isMatrix = false;
        *isVector = false;
    }

	if (matrix->rows > 1 && matrix->cols > 1) { /* This is a classical matrix. */
        if (isMatrix != NULL) {
            *isMatrix = true;
        }
		/* The corners of the matrix are established */
		int firstRow = 0;
		int lastRow = matrix->rows;
		int firstCol = 0;
		int lastCol = matrix->cols;
		for (i = 0; i < matrix->rows; i++) { // From the first row
			if (findNotEqual(matrix, i, -1, 0) != 0) {
				break;
			}
			firstRow++;
        }
		for (i = matrix->rows-1; i >= 0; i--) { // From the last row
			if (findNotEqual(matrix, i, -1, 0) != 0) {
				break;
			}
			lastRow--;
		}
		for (i = 0; i < matrix->cols; i++) { // From the first column
			if (findNotEqual(matrix, -1, i, 0) != 0) {
				break;
			}
			firstCol++;
		}
		for (i = matrix->cols-1; i >= 0; i--) { // From the last column
			if (findNotEqual(matrix, -1, i, 0) != 0) {
				break;
			}
			lastCol--;
		}
		/* The new size of the matrix is calculated */
		newRows = abs(lastRow - firstRow);
		newCols = abs(lastCol - firstCol);
        matrixInt322DStruct result;
        result.rows = newRows;
        result.cols = newCols;
        result.data = calloc(result.rows, sizeof(int32_t*));
		for (i = 0; i < result.rows; i++) {
			result.data[i] = calloc(result.cols, sizeof(int32_t));
		}
		/* We copy the old matrix in the new one */
		for (i = firstRow; i < lastRow; i++) {
			for (j = firstCol; j < lastCol; j++) {
				result.data[i][j] = matrix->data[i][j];
			}
		}
        matrix->rows = result.rows;
        matrix->cols = result.cols;
		return result;
	}
	else if (matrix->rows == 1 || matrix->rows == 1){ /* This is a vector. Be careful. */
        if (isVector != NULL) {
            *isVector = true;
        }
        if (matrix->rows == 1) { /* Single row, roam it horizontally */
			while (matrix->data[countIni] == 0 && countIni < matrix->cols) {
				countIni++; /* The counter increments while 0s exist, to know how many are before de valid values */
			}
			i = matrix->cols;
			while (matrix->data[i] == 0 && i >= 0) {
				countFin++; /* The counter increments while 0s exist, to know how many are after de valid values */
				i--;
			}
            /* With limits on the sides we copy the data */
			newCols = matrix->cols - (countIni + countFin);
            matrixInt322DStruct result;
            result.rows = newCols;
            result.cols = (matrix->cols - countFin);
            result.data = calloc(result.rows, sizeof(int32_t));
            for (i = countIni; i < result.cols; i++) {
				result.data[j][0] = matrix->data[0][i];
				j++;
			}
			return result;
		}
		else { /* Single column, roam it vertically */
			while (matrix->data[countIni] == 0 && countIni < matrix->rows) {
				countIni++; /* The counter increments while 0s exist, to know how many are before de valid values */
			}
			i = matrix->rows;
			while (matrix->data[i] == 0 && i >= 0) {
				countFin++; /* The counter increments while 0s exist, to know how many are after de valid values */
				i--;
			}
            /* With limits on the sides we copy the data */
			newCols = matrix->rows - (countIni + countFin);
            matrixInt322DStruct result;
            result.rows = newCols;
            result.cols = (matrix->cols - countFin);
            result.data = calloc(result.rows, sizeof(int32_t));
			for (i = countIni; i < (matrix->rows - countFin); i++) {
				result.data[j][0] = matrix->data[i][0];
				j++;
			}
            matrix->rows = result.rows;
            matrix->cols = result.cols;
			return result;
		}
	}
	else {
		printf("\tMatriz dimension different from 1 or 2. Exiting.");
	}
	return *matrix;
}

/*
 * A very specific function to reduce the size of the summaries with many
 * involved rounds.
 */
matrixInt323DStruct cutZerosFromMCUsummary(matrixInt323DStruct* MCUSummary){
	int i, j, z, maxRows = 0, maxCols = 0;
    bool isMatrix;
    bool isVector;

    matrixInt322DStruct dimensionsMCU;
    dimensionsMCU.rows = 2;
    dimensionsMCU.cols = MCUSummary->dims;
	dimensionsMCU.data = calloc(dimensionsMCU.rows, sizeof(int32_t*));
	for (i = 0; i < dimensionsMCU.rows; i++) {
		dimensionsMCU.data[i] = calloc(dimensionsMCU.cols, sizeof(int32_t));
	}
	for (i = 0; i < dimensionsMCU.cols; i++) {
        /* We need to separate each matrix on the test */
        matrixInt322DStruct tmpMatrix;
        tmpMatrix.rows = MCUSummary->rows;
        tmpMatrix.cols = MCUSummary->cols;
        tmpMatrix.data = calloc(tmpMatrix.rows, sizeof(int32_t*));
        for (j = 0; j < tmpMatrix.rows; j++) {
            tmpMatrix.data[j] = calloc(tmpMatrix.cols, sizeof(int32_t));
        }
        for (j = 0; j < tmpMatrix.rows; j++) {
            for (z = 0; z < tmpMatrix.cols; z++) {
                tmpMatrix.data[j][z] = MCUSummary->data[j][z][i];
            }
        }
		cutZerosFromArray(&tmpMatrix, &isMatrix, &isVector);
		dimensionsMCU.data[0][i] = tmpMatrix.rows;
		dimensionsMCU.data[1][i] = tmpMatrix.cols;
        /* Its necessary to know the maximun values to adapt the total matrix */
		if (tmpMatrix.rows > maxRows) {
			maxRows = tmpMatrix.rows;
		}
		if (tmpMatrix.cols > maxCols) {
			maxCols = tmpMatrix.cols;
		}
        //liberaInt(tmpMatrix.data, tmpMatrixRows);
	}
	// Una vez guardadas en la matriz de valores todos los nuevos tamaños, se copia la matriz en una matriz nueva
    matrixInt323DStruct simpleMCUSummary;
    simpleMCUSummary.rows = maxRows;
    simpleMCUSummary.cols = maxCols;
    simpleMCUSummary.dims = MCUSummary->dims;
    simpleMCUSummary.data = calloc(simpleMCUSummary.rows, sizeof(int32_t**));
	for (i = 0; i < simpleMCUSummary.rows; i++) {
		simpleMCUSummary.data[i] = calloc(simpleMCUSummary.cols, sizeof(int32_t*));
		for (j = 0; j < simpleMCUSummary.cols; j++) {
			simpleMCUSummary.data[i][j] = calloc(simpleMCUSummary.dims, sizeof(int32_t));
		}
	}

    for (i = 0; i < simpleMCUSummary.rows; i++) {
		for (j = 0; j < simpleMCUSummary.cols; j++) {
			for (z = 0; z < simpleMCUSummary.dims; z++) {
				simpleMCUSummary.data[i][j][z] = MCUSummary->data[i][j][z];
			}
		}
	}
    simpleMCUSummary.rows = maxRows;
    simpleMCUSummary.cols = maxCols;
	//liberaInt(dimensionsMCU,2);
	return simpleMCUSummary;
}

/*
 * Función que recibe un vector y una matriz, devuelve un vector con los elementos que estén en el vector inicial, 
 * pero que no estén en la matriz. Si no, devuelve array vacío
 */
vectorInt32Struct setdiffTraceRule(vectorInt32Struct *vector, matrixInt322DStruct *matrix){
    int i, x, y, cont = 0;
    bool find = false;
    vectorInt32Struct newVector;
    newVector.length = vector->length;
    newVector.data = calloc(newVector.length, sizeof(int32_t));
    for (i = 0; i < newVector.length; i++) {
        for (x = 0; x < matrix->rows; x++) {
            for (y = 0; y < matrix->cols; y++) {
                if (vector->data[i] == matrix->data[x][y]) {
                    find = true;
                }
            }
        }
        if (find == false) {
            newVector.data[cont] = vector->data[i];
            cont++;
        }
        find = false;
    }
    
    newVector.data = realloc(newVector.data, cont*sizeof(int32_t));
    newVector.length = cont;
    return newVector;
}

/*
 * This function joins two or three vectors without ordering their elements
 */
vectorInt32Struct unionVec(vectorInt32Struct *v1, vectorInt32Struct *v2, vectorInt32Struct *v3){
    vectorInt32Struct unionVec;
    int i, j, cont = 0;
    if (v3 == NULL) { // Two vectors union
        unionVec.length = v1->length + v2->length;
        unionVec.data = calloc(unionVec.length, sizeof(int32_t));
        for (i = 0; i < v1->length; i++) { // The fisrt vector is copied without the repeated elements
            if (i == 0) {
                unionVec.data[i] = v1->data
                [i];
                cont++;
            } else {
                j = 0;
                while (v1->data[i] != unionVec.data[j] && j < cont) {
                    j++;
                }
                if (j == cont) { //  There is not equal element
                    unionVec.data[cont] = v1->data[i];
                    cont++;
                }
            }
        }
        for (i = 0; i < v2->length; i++) {
            j = 0;
            while (j < cont && v2->data[i] != unionVec.data[j]) {
                j++;
            }
            if (j == cont) { // There is not equal element
                unionVec.data[cont] = v2->data[i];
                cont++;
            }
        }
        unionVec.data = realloc(unionVec.data, cont*sizeof(int32_t));
    }else { // Three vectors union
        unionVec.length = v1->length + v2->length + v3->length;
        unionVec.data = calloc(unionVec.length, sizeof(int32_t));
        for (i = 0; i < v1->length; i++) { // Se copia el primer vector eliminando los elementos repetidos
            if (i == 0) {
                unionVec.data[i] = v1->data[i];
                cont++;
            } else {
                j = 0;
                while (v1->data[i] != unionVec.data[j] && j < cont) {
                    j++;
                }
                if (j == cont) { // No se ha encontrado elemento igual
                    unionVec.data[cont] = v1->data[i];
                    cont++;
                }
            }
        }
        for (i = 0; i < v2->length; i++) {
            j = 0;
            while (j < cont && v2->data[i] != unionVec.data[j]) {
                j++;
            }
            if (j == cont) { // No se ha encontrado elemento igual
                unionVec.data[cont] = v2->data[i];
                cont++;
            }
        }
        for (i = 0; i < v3->length; i++) {
            j = 0;
            while (j < cont && v3->data[i] != unionVec.data[j]) {
                j++;
            }
            if (j == cont) { // No se ha encontrado elemento igual
                unionVec.data[cont] = v3->data[i];
                cont++;
            }
        }
        unionVec.data = realloc(unionVec.data, cont*sizeof(int32_t));
    }
    unionVec.length = cont;
    return unionVec;
}

/*
# An alternative implementation of the trace rule.

### Anomal XOR values is used to avoid the repetition of elements. It is just
### XORExtractedValues[1,:]
*/
vectorInt32Struct traceRule(matrixInt322DStruct *xorDVtotalrepetitions, matrixUint322DStruct *totalDVhistogram, matrixInt322DStruct *anomalXORvalues, long int LN){
    int i, j, z, index;
    int32_t nBits = (log(LN + 1) / log(2)); // ejemplo nBits = 24
    
	// Elements with 1-trace
	printf("\n\tInvestigating Trace = 1: ");
	//Crea un array de 24 elementos 0 -> 24-1 con valores potencia de 2
    vectorInt32Struct candidatesT1;
    candidatesT1.length = nBits;
	candidatesT1.data = calloc(candidatesT1.length, sizeof(int32_t));
	for (i = 0; i < candidatesT1.length; i++){
		candidatesT1.data[i] = pow(2, i);
	}

	int32_t nT1Threshold = excessiveRepetitions(xorDVtotalrepetitions, 1, LN, "xor", randomnessThreshold);
    
    vectorInt32Struct selectedT1;
    selectedT1.length = nBits;
    selectedT1.data = calloc(selectedT1.length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT1Lenght = 0;
    for (i = 0; i < nBits; i++) { // Se recorre el totalDVhistogram, en las posiciones de los candidatos - 1, pq nuestro histg empieza 0
        if (totalDVhistogram->data[candidatesT1.data[i]-1][1] >= nT1Threshold) {
            selectedT1.data[sT1Lenght] = candidatesT1.data[i];
            sT1Lenght++;
        }
    }
    selectedT1.length = sT1Lenght;
    selectedT1.data = realloc(selectedT1.data, selectedT1.length*sizeof(int32_t));
    free(candidatesT1.data);

    // Devolver array con los elementos que estén en selectedt1, pero que no estén en anomalxorvalues
    vectorInt32Struct winnersT1;
    winnersT1 = setdiffTraceRule(&selectedT1, anomalXORvalues);
    
    free(selectedT1.data);
    printf("%d candidate", winnersT1.length);
    if (winnersT1.length != 1)
        printf("s");
    printf(" found.");
    
    /// Elements with 2-trace
    printf("\n\tInvestigating Trace = 2: ");
    vectorInt32Struct candidatesT2;
    candidatesT2.length = 0.5*nBits*(nBits-1);
    candidatesT2.data = calloc(candidatesT2.length, sizeof(int32_t));
    index = 0;
    for (i = 0; i < (nBits-1); i++) {
        for (j = i+1; j < nBits; j++) {
            candidatesT2.data[index] = (pow(2, i))+(pow(2, j)); // CandidatesT2[index]=2^k1+2^k2
            index++;
        }
    }
    
    int32_t nT2Threshold = excessiveRepetitions(xorDVtotalrepetitions, 1, LN, "xor", randomnessThreshold);
   
    vectorInt32Struct selectedT2;
    selectedT2.length = candidatesT2.length;
    selectedT2.data = calloc(candidatesT2.length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT2Lenght = 0;
    for (i = 0; i < candidatesT2.length; i++) {
        if (totalDVhistogram->data[candidatesT2.data[i]-1][1] >= nT2Threshold) {
            selectedT2.data[sT2Lenght] = candidatesT2.data[i];
            sT2Lenght++;
        }
    }
    selectedT2.length = sT2Lenght;
    selectedT2.data = realloc(selectedT2.data, sT2Lenght*sizeof(int32_t));
    free(candidatesT2.data);
    
    vectorInt32Struct winnersT2;
    winnersT2 = setdiffTraceRule(&selectedT2, anomalXORvalues);
    
    free(selectedT2.data);
    printf("%d candidate", winnersT2.length);
    if (winnersT2.length != 1)
        printf("s");
    printf(" found.");
    
    /// Elements with 3-trace
    printf("\n\tInvestigating Trace = 3: ");
    vectorInt32Struct candidatesT3;
    candidatesT3.length = nBits*(nBits-1)*(nBits-2)/6;
    candidatesT3.data = calloc(candidatesT3.length, sizeof(int32_t));
    index = 0;
    for (i = 0; i < (nBits-2); i++) {
        for (j = i+1; j < (nBits-1); j++) {
            for (z = j+1; z < nBits; z++) {
                candidatesT3.data[index] = (pow(2, i))+(pow(2, j))+(pow(2, z)); // 2^k1+2^k2+2^k3
                index++;
            }
        }
    }
    vectorInt32Struct selectedT3;
    selectedT3.data = calloc(candidatesT3.length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    int sT3Lenght = 0;
    for (i = 0; i < candidatesT3.length; i++) {
        if (totalDVhistogram->data[candidatesT3.data[i]-1][1] >= nT2Threshold) {
            selectedT3.data[sT3Lenght] = candidatesT3.data[i];
            sT3Lenght++;
        }
    }
    selectedT3.length = sT3Lenght;
    selectedT3.data = realloc(selectedT3.data, selectedT3.length*sizeof(int32_t));
    
    free(candidatesT3.data);
    
    vectorInt32Struct winnersT3;
    winnersT3 = setdiffTraceRule(&selectedT3, anomalXORvalues);
    printf("%d candidate", winnersT3.length);
    free(selectedT3.data);
    if (winnersT3.length != 1)
        printf("s");
    printf(" found.");
    printf("\n\tEnd of search.\n");
    
    vectorInt32Struct vectorUnion = unionVec(&winnersT1, &winnersT2, &winnersT3);
    free(winnersT1.data);
    free(winnersT2.data);
    free(winnersT3.data);
    
    return vectorUnion;
}


matrixInt323DStruct propose_MCUs(char* op, matrixUint323DStruct *XORDVmatrix3D, matrixUint323DStruct *POSDVmatrix3D, matrixInt322DStruct *anomalXOR, matrixInt322DStruct *anomalPOS, vectorIntStruct* nAddressesInRound){
	int i, j, ktest;
	int rows = ceil(XORDVmatrix3D->rows / 2 - 1);
	if (strcmp(op, "xor") == 0){ /* OP == XOR */
        matrixInt323DStruct XOR_summary;
        XOR_summary.rows = rows;
        XOR_summary.cols = UnrealMCUsize;
        XOR_summary.dims = XORDVmatrix3D->dims;
        XOR_summary.data = calloc(XOR_summary.rows, sizeof(int32_t**));
		for (i = 0; i < XOR_summary.rows; i++) {
			XOR_summary.data[i] = calloc(XOR_summary.cols, sizeof(int32_t*));
			for (j = 0; j < XOR_summary.cols; j++) {
				XOR_summary.data[i][j] = calloc(XOR_summary.dims, sizeof(int32_t));
			}
		}
		for (ktest = 0; ktest < XOR_summary.dims; ktest++){
            matrixUint322DStruct xordvmatrix;
            xordvmatrix.rows = nAddressesInRound->data[ktest];
            xordvmatrix.cols = nAddressesInRound->data[ktest];
            xordvmatrix.data = calloc(xordvmatrix.rows, sizeof(int32_t*));
			for (i = 0; i < xordvmatrix.rows; i++) {
				xordvmatrix.data[i] = calloc(xordvmatrix.cols, sizeof(int32_t));
			}
            copyOfMatrix3to2(XORDVmatrix3D, &xordvmatrix, ktest);
            
			matrixBool2DStruct XORmarked_pairs = marking_addresses(&xordvmatrix, anomalXOR);
        
			/* The following results only concern an operation.Just present for illustrative
			 * operations.In practical use, only CMB should be used. */
			matrixInt322DStruct XORproposed_MCUs = agrupate_mcus(&XORmarked_pairs);
			for (i = 0; i < XORproposed_MCUs.rows; i++) {
				for (j = 0; j < XORproposed_MCUs.cols; j++) {
					XOR_summary.data[i][j][ktest] = XORproposed_MCUs.data[i][j];
				}
			}
		}
        
        /* However, as there are too many zeros, it is interesting to reshape the matrix.
         * Thus, some memory is saved.But we must be careful with the different dimensions
         * of every partial XY matrix in XYZ arrays.A simple but unelegant way is : */
        matrixInt323DStruct finalXOR_summary = cutZerosFromMCUsummary(&XOR_summary);
        libera3D(XOR_summary.data, XOR_summary.rows, XOR_summary.cols);
        return finalXOR_summary;
	}
	else if (strcmp(op, "pos") == 0){ /* OP == POS */
        matrixInt323DStruct POS_summary;
        POS_summary.rows = rows;
        POS_summary.cols = UnrealMCUsize;
        POS_summary.dims = XORDVmatrix3D->dims;
        POS_summary.data = calloc(POS_summary.rows, sizeof(int32_t**));
		for (i = 0; i < POS_summary.rows; i++) {
			POS_summary.data[i] = calloc(POS_summary.cols, sizeof(int32_t*));
			for (j = 0; j < POS_summary.cols; j++) {
				POS_summary.data[i][j] = calloc(POS_summary.dims, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < XORDVmatrix3D->dims; ktest++){
            matrixInt322DStruct posdvmatrix;
            posdvmatrix.rows = nAddressesInRound->data[ktest];
            posdvmatrix.cols = nAddressesInRound->data[ktest];
            posdvmatrix.data = calloc(posdvmatrix.rows, sizeof(int32_t*));
			for (i = 0; i < posdvmatrix.rows; i++) {
				posdvmatrix.data[i] = calloc(posdvmatrix.cols, sizeof(int32_t));
			}
            copyOfMatrix3to2(POSDVmatrix3D, &posdvmatrix, ktest);
            
            matrixBool2DStruct POSmarked_pairs;
            POSmarked_pairs = marking_addresses(&posdvmatrix, anomalPOS);

			/* The following results only concern an operation.Just present for illustrative
			 * operations.In practical use, only CMB should be used. */
			matrixInt322DStruct POSproposed_MCUs = agrupate_mcus(&POSmarked_pairs);
			for (i = 0; i < POSproposed_MCUs.rows; i++) {
				for (j = 0; j < POSproposed_MCUs.rows; j++) {
					POS_summary.data[i][j][ktest] = POSproposed_MCUs.data[i][j];
				}
			}
		}
        /* However, as there are too many zeros, it is interesting to reshape the matrix.
         * Thus, some memory is saved.But we must be careful with the different dimensions
         * of every partial XY matrix in XYZ arrays.A simple but unelegant way is : */
        matrixInt323DStruct finalPOS_summary = cutZerosFromMCUsummary(&POS_summary);
        libera3D(POS_summary.data, POS_summary.rows, POS_summary.cols);
        return finalPOS_summary;
	}
	else{ /* OP == CMB */
        matrixInt323DStruct CMB_summary;
        CMB_summary.rows = rows;
        CMB_summary.cols = UnrealMCUsize;
        CMB_summary.dims = XORDVmatrix3D->dims;
        CMB_summary.data = calloc(CMB_summary.rows, sizeof(int32_t**));
		for (i = 0; i < CMB_summary.rows; i++) {
			CMB_summary.data[i] = calloc(CMB_summary.cols, sizeof(int32_t*));
			for (j = 0; j < CMB_summary.cols; j++) {
				CMB_summary.data[i][j] = calloc(CMB_summary.dims, sizeof(int32_t));
			}
		}

		for (ktest = 0; ktest < XORDVmatrix3D->dims; ktest++){
            matrixInt322DStruct xordvmatrix;
            xordvmatrix.rows = nAddressesInRound->data[ktest];
            xordvmatrix.cols = nAddressesInRound->data[ktest];
            xordvmatrix.data = calloc(xordvmatrix.rows, sizeof(int32_t*));
			for (i = 0; i < xordvmatrix.rows; i++) {
				xordvmatrix.data[i] = calloc(xordvmatrix.cols, sizeof(int32_t));
			}
            copyOfMatrix3to2(XORDVmatrix3D, &xordvmatrix, ktest);
            
            matrixInt322DStruct posdvmatrix;
            posdvmatrix.rows = nAddressesInRound->data[ktest];
            posdvmatrix.cols = nAddressesInRound->data[ktest];
            posdvmatrix.data = calloc(posdvmatrix.rows, sizeof(int32_t*));
			for (i = 0; i < posdvmatrix.rows; i++) {
				posdvmatrix.data[i] = calloc(posdvmatrix.cols, sizeof(int32_t));

			}
            copyOfMatrix3to2(POSDVmatrix3D, &posdvmatrix, ktest);
            
            matrixBool2DStruct XORmarked_pairs = marking_addresses(&xordvmatrix, anomalXOR);
			matrixBool2DStruct POSmarked_pairs = marking_addresses(&posdvmatrix, anomalPOS);

            matrixBool2DStruct XOR_POS_marked_pairs;
            XOR_POS_marked_pairs.rows = nAddressesInRound->data[ktest];
            XOR_POS_marked_pairs.cols = nAddressesInRound->data[ktest];
            XOR_POS_marked_pairs.data = calloc(nAddressesInRound->data[ktest], sizeof(bool*));
			for (i = 0; i < nAddressesInRound->data[ktest]; i++) {
				XOR_POS_marked_pairs.data[i] = calloc(nAddressesInRound->data[ktest], sizeof(bool));
			}

			for (i = 0; i < nAddressesInRound->data[ktest]; i++) {
				for (j = 0; j < nAddressesInRound->data[ktest]; j++) {
					XOR_POS_marked_pairs.data[i][j] = XORmarked_pairs.data[i][j] | POSmarked_pairs.data[i][j];
				}
			}

			matrixInt322DStruct CMBproposed_MCUs = agrupate_mcus(&XOR_POS_marked_pairs);
			for (i = 0; i < CMBproposed_MCUs.rows; i++) {
				for (j = 0; j < CMBproposed_MCUs.cols; j++) {
					CMB_summary.data[i][j][ktest] = CMBproposed_MCUs.data[i][j];
                }
			}
		}
        /* However, as there are too many zeros, it is interesting to reshape the matrix.
         * Thus, some memory is saved.But we must be careful with the different dimensions
         * of every partial XY matrix in XYZ arrays.A simple but unelegant way is : */
        matrixInt323DStruct finalCMB_summary;
        finalCMB_summary = cutZerosFromMCUsummary(&CMB_summary);
        libera3D(CMB_summary.data, CMB_summary.rows, CMB_summary.cols);
        return finalCMB_summary;
	}
}

matrixInt322DStruct condensate_summary(matrixInt323DStruct* summary3D, vectorIntStruct* nAddressesInRound){
	int i, j, k = 1, sum = 0;
    matrixInt322DStruct condensateSummary;
    condensateSummary.rows = summary3D->cols;
    condensateSummary.cols = summary3D->dims;
    condensateSummary.data = calloc(condensateSummary.rows, sizeof(int32_t*));
	for (i = 0; i < condensateSummary.rows; i++) {
		condensateSummary.data[i] = calloc(condensateSummary.cols + 1, sizeof(int32_t));
        condensateSummary.data[i][0] = k;
        k++;
	}
	for (i = 1; i < summary3D->dims + 1; i++) {
		// elems tendrá la longitud del nuevo vector
		int nEvents = 0;
        matrixInt322DStruct summary2D;
        summary2D.rows = summary3D->rows;
        summary2D.cols = summary3D->cols;
        summary2D.data = calloc(summary2D.rows, sizeof(int32_t*));
        for (k = 0; k < summary2D.rows; k++) {
            summary2D.data[k] = calloc(summary2D.cols, sizeof(int32_t*));
        }
        copyOfMatrix3to2(summary3D, &summary2D, (i-1));

		for (k = summary3D->cols-1; k > 0; k--) {
            nEvents = findNotEqual(&summary2D, -1, k, 0);
            for (j = k; j < summary3D->cols; j++) {
                sum += condensateSummary.data[j][i];
            }
            condensateSummary.data[k][i] = nEvents - sum;
            sum = 0;
		}
		// Apparently, the following solution is a bit faster.
        condensateSummary.data[0][i] = nAddressesInRound->data[i-1] - findNotEqualMatrix(&summary2D, 0);
	}
    condensateSummary.rows = summary3D->cols;
    condensateSummary.cols = summary3D->dims + 1;
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
matrixInt322DStruct locate_mbus(matrixInt322DStruct* content, int** pattern, int nRoundsInPattern, int datawidth){
	int i, j, nDetected = 0;
    matrixInt322DStruct summary;
    summary.rows = 3 * content->rows;
    summary.cols = 4;
    summary.data = calloc(summary.rows, sizeof(uint32_t*));
	for (i = 0; i < summary.rows; i++) {
		summary.data[i] = calloc(summary.cols, sizeof(uint32_t));
	}

	for (i = 0; i < nRoundsInPattern; i++) {
		for (j = 0; j < content->rows; j++) {
			if (content->data[j][3 * i] == 0xFFFFFFFF) {
				break;
			}
			else {
				vectorInt32Struct corruptedBits = flipped_bits(content->data[j][3 * i + 1], pattern[i][1], datawidth);
				if (corruptedBits.length > 1){
                    summary.data[nDetected][0] = i;
                    summary.data[nDetected][1] = j;
                    summary.data[nDetected][2] = content->data[j][3*i];
                    summary.data[nDetected][3] = corruptedBits.length;
                    nDetected++;
				}
			}
		}
	}

    return cutZerosFromArray(&summary, NULL, NULL);
}

/*
 *
 */
matrixInt322DStruct extractAnomalDVSelfConsistency(char* op, matrixInt322DStruct* opDVtotalrepetitions, matrixUint322DStruct* totalDVhistogram, vectorIntStruct* nAddressesInRound, matrixInt322DStruct* opdvhistogram, matrixInt323DStruct* opdvmatrixbackup, int* XORANOMALS){
	int32_t opNthreshold;
	int* n_anomalous_values = malloc(sizeof(int));
    matrixInt322DStruct testOPa;
    testOPa.rows = 0;
    testOPa.cols = 0;
	matrixInt322DStruct oPExtracted_values;
	bool oPSelfConsistence = true;
	int n_anomalous_repetitions = 1, i, j, test;
	printf("\n\tDetermining the threshold for repetition excess:\n");
    
    if (strcmp(op, "xor") == 0) {
		printf("\n\t\tXOR operation...");
		opNthreshold = excessiveRepetitions(opDVtotalrepetitions, 1, totalDVhistogram->rows, "xor", randomnessThreshold);
        // Redondeo hacia arriba posible solución?
        testOPa = find_anomalies_histogram(totalDVhistogram, 1, opDVtotalrepetitions, 1, opNthreshold, n_anomalous_values);
        *XORANOMALS = *n_anomalous_values;

		printf("\n\tXOR operation:\n");
		oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, 1);

		while (oPSelfConsistence) {
			printf("\t\tStep %d ", n_anomalous_repetitions);
			oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, n_anomalous_repetitions);
            
            matrixInt322DStruct opPartRepsSelfCons;
            opPartRepsSelfCons.rows = oPExtracted_values.rows;
            opPartRepsSelfCons.cols = nAddressesInRound->length;
            opPartRepsSelfCons.data = calloc(opPartRepsSelfCons.rows, sizeof(int32_t*));
			for (i = 0; i < opPartRepsSelfCons.rows; i++) {
				opPartRepsSelfCons.data[i] = calloc(opPartRepsSelfCons.cols, sizeof(int32_t));
			}
			for (i = 0; i < opPartRepsSelfCons.rows; i++) {
                for (j = 0; j < opPartRepsSelfCons.cols; j++) {
					opPartRepsSelfCons.data[i][j] = opdvhistogram->data[(oPExtracted_values.data[i][0])-1][j+1];
				}
			}
            
			for (test = 0; test < nAddressesInRound->length; test++) {
				printf("Test %d", test+1);
                matrixUint322DStruct opdvmatrix;
                opdvmatrix.rows = nAddressesInRound->data[test];
                opdvmatrix.cols = nAddressesInRound->data[test];
                opdvmatrix.data = calloc(opdvmatrix.rows, sizeof(int32_t*));
				for (i = 0; i < opdvmatrix.rows; i++) {
					opdvmatrix.data[i] = calloc(opdvmatrix.cols, sizeof(int32_t));
				}
                copyOfMatrix3to2(opdvmatrixbackup, &opdvmatrix, test);
                matrixBool2DStruct opMarkedPairs;
                opMarkedPairs = marking_addresses(&opdvmatrix, &oPExtracted_values);
                
                matrixInt322DStruct opProposedMcus;
                opProposedMcus = agrupate_mcus(&opMarkedPairs);
                int largestMCUSize = opProposedMcus.cols;
                
                int continuation = findEqual(&opPartRepsSelfCons, 0, test, largestMCUSize);
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
                oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, n_anomalous_repetitions);
                    
            }
        }

	} else if (strcmp(op, "pos") == 0){
		printf("\n\t\tPOS operation (It can take a long if there are too many addresses)...\n");
		opNthreshold = excessiveRepetitions(opDVtotalrepetitions, 1, totalDVhistogram->rows, "pos", randomnessThreshold);
        testOPa = find_anomalies_histogram(totalDVhistogram, 2, opDVtotalrepetitions, 1, opNthreshold, n_anomalous_values);

		printf("\n\tPOS operation:\n");
        oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, 1);
    
        while (oPSelfConsistence) {
			printf("\t\tStep %d ", n_anomalous_repetitions);
            oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, n_anomalous_repetitions);
            
            matrixInt322DStruct opPartRepsSelfCons;
            opPartRepsSelfCons.rows = oPExtracted_values.rows;
            opPartRepsSelfCons.cols = nAddressesInRound->length;
            opPartRepsSelfCons.data = calloc(opPartRepsSelfCons.rows, sizeof(int32_t*));
            for (i = 0; i < opPartRepsSelfCons.rows; i++) {
                opPartRepsSelfCons.data[i] = calloc(opPartRepsSelfCons.cols, sizeof(int32_t));
            }
            for (i = 0; i < opPartRepsSelfCons.rows; i++) {
                for (j = 0; j < opPartRepsSelfCons.cols; j++) {
                    opPartRepsSelfCons.data[i][j] = opdvhistogram->data[(oPExtracted_values.data[i][0])-1][j+1];
                }
            }
            
			for (test = 0; test < nAddressesInRound->length; test++) {
				printf("Test %d", test);

                matrixUint322DStruct opdvmatrix;
                opdvmatrix.rows = nAddressesInRound->data[test];
                opdvmatrix.cols = nAddressesInRound->data[test];
                opdvmatrix.data = calloc(opdvmatrix.rows, sizeof(int32_t*));
                for (i = 0; i < opdvmatrix.rows; i++) {
                    opdvmatrix.data[i] = calloc(opdvmatrix.cols, sizeof(int32_t));
                }
                copyOfMatrix3to2(opdvmatrixbackup, &opdvmatrix, test);
                matrixBool2DStruct opMarkedPairs = marking_addresses(&opdvmatrix, &oPExtracted_values);
                
                matrixInt322DStruct opProposedMcus = agrupate_mcus(&opMarkedPairs);
                int largestMCUSize = opProposedMcus.cols;
                
                int continuation = findEqual(&opPartRepsSelfCons, 0, test, largestMCUSize);
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
                oPExtracted_values = extract_some_critical_values(&testOPa, n_anomalous_values, n_anomalous_repetitions);
            }
        }

	}
	return oPExtracted_values;
}

matrixInt323DStruct index2address(matrixInt323DStruct* indexMatrix, matrixInt322DStruct* addressMatrix){
	int i, j, z;
    matrixInt323DStruct result;
    result.rows = indexMatrix->rows;
    result.cols = indexMatrix->cols;
    result.dims = indexMatrix->dims;
    result.data = calloc(result.rows, sizeof(int32_t**));
	for (i = 0; i < result.rows; i++) {
		result.data[i] = calloc(result.cols, sizeof(int32_t*));
		for (j = 0; j < result.cols; j++) {
			result.data[i][j] = calloc(result.dims, sizeof(int32_t));
		}
	}
	for (z = 0; z < indexMatrix->dims; z++) {
		for (i = 0; i < indexMatrix->rows; i++) {
			for (j = 0; j < indexMatrix->cols; j++) {
				if (indexMatrix->data[i][j][z] == 0) {
					break;
				}
				else {
					result.data[i][j][z] = addressMatrix->data[indexMatrix->data[i][j][z]][z];
				}
			}
		}
	}
	return result;
}

vectorInt32Struct criticalXORValuesFromClusters(matrixInt323DStruct* PrelMCUSummary, matrixInt322DStruct* AddressMatrix, matrixInt322DStruct* XORextracted_values, matrixInt322DStruct* xorDVtotalrepetitions, matrixInt322DStruct* DVHistogram){
    vectorInt32Struct newCandidates;
    newCandidates.length = 0;
    newCandidates.data = calloc(PrelMCUSummary->rows * PrelMCUSummary->cols * PrelMCUSummary->dims, sizeof(int32_t));
	int ktest, kRow, kAddress1, kAddress2, i;
    int32_t index1, index2, address1, address2, candidate;
	for (ktest = 0; ktest < PrelMCUSummary->dims; ktest++){
		for (kRow = 0; kRow < PrelMCUSummary->rows; kRow++){
			for (kAddress1 = 1; kAddress1 < PrelMCUSummary->cols; kAddress1++){
				index1 = PrelMCUSummary->data[kRow][kAddress1][ktest];
				if (index1 == 0){
					break;
				}
				else{
					address1 = AddressMatrix->data[index1-1][ktest*3];
					for (kAddress2 = 0; kAddress2 < kAddress1; kAddress2++){
						index2 = PrelMCUSummary->data[kRow][kAddress2][ktest];
						address2 = AddressMatrix->data[index2-1][ktest*3];
						candidate = address1 ^ address2;
						if (findEqualMatrix(XORextracted_values, candidate) == 0){
                            vectorInt32Struct vCandidate;
                            vCandidate.length = 1;
                            vCandidate.data = calloc(1, sizeof(int32_t));
                            vCandidate.data[0] = candidate;
                            newCandidates = unionVec(&newCandidates, &vCandidate, NULL);
                            free(vCandidate.data);
						}
					}
				}		
			}
		}		
    }
    newCandidates.data = realloc(newCandidates.data, newCandidates.length * sizeof(int32_t));
    
    printf("\n");
    for (i = 0; i < newCandidates.length; i++) {
        printf("%x ", newCandidates.data[i]);
    }

	int32_t xorNthreshold = excessiveRepetitions(xorDVtotalrepetitions, 1, DVHistogram->rows, "xor", randomnessThreshold);

    vectorInt32Struct purgedCandidatesXOR;
    purgedCandidatesXOR.length = newCandidates.length;
    purgedCandidatesXOR.data = calloc(purgedCandidatesXOR.length, sizeof(int32_t));
    int purgedCandidatesLenght = 0;
    for (i = 0; i < purgedCandidatesXOR.length; i++) {
        if (DVHistogram->data[newCandidates.data[i]-1][1] >= xorNthreshold) {
            purgedCandidatesXOR.data[purgedCandidatesLenght] = newCandidates.data[i];
            purgedCandidatesLenght++;
        }
    }
    purgedCandidatesXOR.length = purgedCandidatesLenght;
    purgedCandidatesXOR.data = realloc(purgedCandidatesXOR.data, purgedCandidatesXOR.length*sizeof(int32_t));
    
	return purgedCandidatesXOR;
}

void extractAnomalDVfromClusters(matrixInt322DStruct* content, matrixInt322DStruct* XORextracted_values, matrixInt322DStruct* POSextracted_values, matrixInt322DStruct* xorDVtotalrepetitions, matrixInt322DStruct* posDVtotalrepetitions, matrixUint322DStruct* totalDVhistogram, matrixInt323DStruct* xordvmatrixbackup, matrixInt323DStruct* posdvmatrixbackup, vectorIntStruct* nAddressesInRound, long int LN, vectorInt32Struct* discoveredXORDVs, vectorInt32Struct* discoveredPOSDVs, matrixInt322DStruct* tempXORDVvalues, matrixInt322DStruct* tempPOSDVvalues){

    bool foundNewDVValues = true;
    int step = 0, i, j;
    
    tempXORDVvalues->rows = XORextracted_values->rows;
    tempXORDVvalues->cols = XORextracted_values->cols;
    tempXORDVvalues->data = calloc(tempXORDVvalues->rows, sizeof(int32_t*));
    for (i = 0; i < tempXORDVvalues->rows; i++) {
        tempXORDVvalues->data[i] = calloc(tempXORDVvalues->cols, sizeof(int32_t));
    }

    tempPOSDVvalues->rows = POSextracted_values->rows;
    tempPOSDVvalues->cols = POSextracted_values->cols;
    tempPOSDVvalues->data = calloc(tempPOSDVvalues->rows, sizeof(int32_t*));
    for (i = 0; i < tempPOSDVvalues->rows; i++) {
        tempPOSDVvalues->data[i] = calloc(tempPOSDVvalues->cols, sizeof(int32_t));
    }

    //COPIA OBLIGADA, Se podría hacer función de copia de matriz simple a otra matriz simple
    for (i = 0; i < tempXORDVvalues->rows; i++) {
        for (j = 0; j < tempXORDVvalues->cols; j++) {
            tempXORDVvalues->data[i][j] = XORextracted_values->data[i][j];
        }
    }
    for (i = 0; i < tempPOSDVvalues->rows; i++) {
        for (j = 0; j < tempPOSDVvalues->cols; j++) {
            tempPOSDVvalues->data[i][j] = POSextracted_values->data[i][j];
        }
    }
    
    while (foundNewDVValues) {
        step++;
        printf("\nStep: %d", step);
        printf("Creating preliminary organization...");
        matrixInt323DStruct tempCMB_summary;
        matrixInt323DStruct tempXOR_summary;
        matrixInt323DStruct tempPOS_summary;
        
        tempXOR_summary = propose_MCUs("xor", xordvmatrixbackup, posdvmatrixbackup, tempXORDVvalues, tempPOSDVvalues, nAddressesInRound);
        tempPOS_summary = propose_MCUs("pos", xordvmatrixbackup, posdvmatrixbackup, tempXORDVvalues, tempPOSDVvalues, nAddressesInRound);
        tempCMB_summary = propose_MCUs("cmb", xordvmatrixbackup, posdvmatrixbackup, tempXORDVvalues, tempPOSDVvalues, nAddressesInRound);
        
       
        printf("\n");
        for (int i = 0; i < tempXOR_summary.rows; i++) {
            for (int j = 0; j < tempXOR_summary.cols; j++) {
                printf("%d ", tempXOR_summary.data[i][j][0]);
            }
            printf("\n");
        }
        
        printf("\n");
        for (int i = 0; i < tempXOR_summary.rows; i++) {
            for (int j = 0; j < tempXOR_summary.cols; j++) {
                printf("%d ", tempXOR_summary.data[i][j][1]);
            }
            printf("\n");
        }
        
        printf("\n");
        for (int i = 0; i < tempXOR_summary.rows; i++) {
            for (int j = 0; j < tempXOR_summary.cols; j++) {
                printf("%d ", tempXOR_summary.data[i][j][2]);
            }
            printf("\n");
        }
        
        
        
        printf(" Ended. Searching new Elements...");

        printf(" XOR...");

        int nProposedXORDV = 0;
        vectorInt32Struct proposedXORDV;
        proposedXORDV = criticalXORValuesFromClusters(&tempCMB_summary, content, tempXORDVvalues, xorDVtotalrepetitions, totalDVhistogram);
		
            for (int j = 0; j < proposedXORDV.length; j++) {
                printf("%d ", proposedXORDV.data[j]);
            }
            printf("\n");
    
        
        vectorInt32Struct proposedPOSDV;
        nProposedXORDV = proposedXORDV.length;
        int nProposedPOSDV = 0;
        printf(" POS. SUB...");

        printf(" Ended. ");
        
        if (nProposedXORDV != 0) {
            printf("\n\t\tXOR operation:  %d new DV element", nProposedXORDV);
            if (nProposedXORDV != 1) {
                printf("s");
            }
            printf("found.");
            *discoveredXORDVs = unionVec(&proposedXORDV, discoveredXORDVs, NULL);

            matrixInt322DStruct auxXORDV;
            auxXORDV.rows = XORextracted_values->rows + nProposedXORDV;
            auxXORDV.cols = XORextracted_values->cols;
            auxXORDV.data = calloc(auxXORDV.rows, sizeof(int32_t*));
            for (i = 0; i < auxXORDV.rows; i++) {
                auxXORDV.data[i] = calloc(auxXORDV.cols, sizeof(int32_t));
            }
            for (i = 0; i < XORextracted_values->rows; i++) {
                for (j = 0; j < XORextracted_values->cols; j++) {
                    auxXORDV.data[i][j] = tempXORDVvalues->data[i][j];
                }
            }
            int index = XORextracted_values->rows;
            for (int a = 0; a < nProposedXORDV; a++) {
                auxXORDV.data[index][0] = proposedXORDV.data[a];
                auxXORDV.data[index][1] = totalDVhistogram->data[proposedXORDV.data[a]-1][1];
                index++;
            }

            tempXORDVvalues->data = auxXORDV.data;

            tempXORDVvalues->rows = XORextracted_values->rows + nProposedXORDV;
            XORextracted_values->rows = XORextracted_values->rows + nProposedXORDV;

        }
        if (nProposedPOSDV != 0) {
            printf("\n\t\tPOS operation:  %d new DV element", nProposedPOSDV);
            if (nProposedPOSDV != 1) {
                printf("s");
            }
            printf("found.");
            *discoveredPOSDVs = unionVec(&proposedPOSDV, discoveredPOSDVs, NULL);

            matrixInt322DStruct auxPOSDV;
            auxPOSDV.rows = POSextracted_values->rows + nProposedPOSDV;
            auxPOSDV.cols = POSextracted_values->cols;
            auxPOSDV.data = calloc(auxPOSDV.rows, sizeof(int32_t*));
            for (i = 0; i < auxPOSDV.rows; i++) {
                auxPOSDV.data[i] = calloc(auxPOSDV.cols, sizeof(int32_t));
            }
            for (i = 0; i < auxPOSDV.rows; i++) {
                for (j = 0; j < auxPOSDV.cols; j++) {
                    auxPOSDV.data[i][j] = tempPOSDVvalues->data[i][j];
                }
            }
            int index = POSextracted_values->rows;
            for (int a = 0; a < nProposedPOSDV; a++) {
                auxPOSDV.data[index][0] = proposedPOSDV.data[a];
                auxPOSDV.data[index][1] = totalDVhistogram->data[proposedPOSDV.data[a]-1][2];
                index++;
            }

            tempPOSDVvalues->data = auxPOSDV.data;
            tempXORDVvalues->rows = POSextracted_values->rows + nProposedPOSDV;
            POSextracted_values->rows = POSextracted_values->rows + nProposedPOSDV;
            
        }
        if ((nProposedXORDV == 0) && (nProposedPOSDV == 0)) {
            foundNewDVValues = false;
            printf("\n\n\t\tNO MORE ELEMENTS DISCOVERED. Stopping analysis.");
        }

    }
}


matrixInt322DStruct criticalXORvaluesFromXORingRule(matrixInt322DStruct* XORextracted_values, matrixInt322DStruct* XORDVtotalrepetitions, matrixUint322DStruct* totalDVhistogram, vectorInt32Struct* candidates){
	int i, j;
	bool candidateFind = false, prelCandidateFind = false;
	printf("\n\tCase 1: XORing known values... ");
    vectorInt32Struct prelCandidates;
    prelCandidates.data = calloc(1, sizeof(int32_t));
    prelCandidates.length = 0;

    vectorInt32Struct oldXORValues;
    oldXORValues.length = XORextracted_values->rows;
    oldXORValues.data = calloc(oldXORValues.length, sizeof(int32_t));
    for (i = 0; i < oldXORValues.length; i++) {
        oldXORValues.data[i] = XORextracted_values->data[i][0];
    }

    sort(oldXORValues);

	int k1, k2;
	for (k1 = 0; k1 < XORextracted_values->rows; k1++){
		int old1 = oldXORValues.data[k1];
		for (k2 = k1 + 1; k2 < XORextracted_values->rows; k2++){
			candidateFind = false, prelCandidateFind = false;
			int old2 = oldXORValues.data[k2];
			int candidate = old1^old2;
            if ((findEqualVector(&oldXORValues, candidate) == 0) && (findEqualVector(&prelCandidates, candidate) == 0)) {
                vectorInt32Struct vCandidate;
                vCandidate.length = 1;
                vCandidate.data = calloc(1, sizeof(int32_t));
                vCandidate.data[0] = candidate;
                prelCandidates = unionVec(&prelCandidates, &vCandidate, NULL);
            }
		}
	}
	sort(prelCandidates);

	int nThresholdCase1 = excessiveRepetitions(XORDVtotalrepetitions, 1, totalDVhistogram->rows, "xor", randomnessThreshold);

    vectorInt32Struct candidatesCase1;
    candidatesCase1.length = 0;
    candidatesCase1.data = calloc(prelCandidates.length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    for (i = 0; i < prelCandidates.length; i++) {
        if (totalDVhistogram->data[prelCandidates.data[i]-1][1] >= nThresholdCase1) {
            candidatesCase1.data[candidatesCase1.length] = prelCandidates.data[i];
            candidatesCase1.length++;
        }
    }
    candidatesCase1.data = realloc(candidatesCase1.data, candidatesCase1.length*sizeof(int32_t));

	if (candidatesCase1.length > 0){
		printf("%d new value", candidatesCase1.length);
		if (candidatesCase1.length != 1)
			printf("s");
		printf(" obtained. \n");
	}
	else{
		printf("No new values. Going on.\n");
		}

    oldXORValues = unionVec(&oldXORValues, &candidatesCase1, NULL);
	sort(oldXORValues);

    /* Case 2: XORingDV elements with confirmed XORs to obtain another one.*/

    printf("\tCASE 2: XORing DV elements with confirmed values... ");
    vectorInt32Struct DVElements;
    DVElements.length = totalDVhistogram->rows;
    DVElements.data = calloc(DVElements.length, sizeof(int32_t));
		int index = 0;
		for (i = 0; i < DVElements.length; i++){
			if (totalDVhistogram->data[i][1] != 0){
				DVElements.data[index] = i;
				index++;
			}
		}
    // This allows reobtaining the elements in the DV set.

    // This +1 is placed since, sometimes, XORED can be 0.
    vectorInt32Struct selectedPrelCandidates;
    selectedPrelCandidates.length = 0;
    selectedPrelCandidates.data = calloc(1, sizeof(int32_t));
    printf(" Searching... ");
    int kdv, kxor;
	bool xoredFind = false, kdvFind = false;
	for (kdv = 0; kdv < index; kdv++){
		for (kxor = 0; kxor < oldXORValues.length; kxor++){
			int xored = DVElements.data[kdv] ^ oldXORValues.data[kxor];
			//presencePrelCandidates[xored + 1] = true;
               if ((findEqualVector(&oldXORValues, xored) != 0) && (findEqualVector(&oldXORValues, kdv) == 0)) {
                   vectorInt32Struct vKdv;
                   vKdv.length = 1;
                   vKdv.data = calloc(1, sizeof(int32_t));
                   vKdv.data[0] = kdv;
                   selectedPrelCandidates = unionVec(&selectedPrelCandidates, &vKdv, NULL);
                }
			}
		}

    int nThresholdCase2 = excessiveRepetitions(XORDVtotalrepetitions, 1, totalDVhistogram->rows, "xor", randomnessThreshold);

    vectorInt32Struct candidatesCase2;
    candidatesCase2.length = 0;
    candidatesCase2.data = calloc(selectedPrelCandidates.length, sizeof(int32_t)); //Luego se redimensiona al número de elemntos real
    for (i = 0; i < selectedPrelCandidates.length; i++) {
        if (totalDVhistogram->data[prelCandidates.data[i]-1][1] >= nThresholdCase2) {
            candidatesCase2.data[candidatesCase2.length] = selectedPrelCandidates.data[i];
            candidatesCase2.length++;
        }
    }
    candidatesCase2.data = realloc(candidatesCase2.data, candidatesCase2.length*sizeof(int32_t));

    if (candidatesCase2.length > 0){
        printf("%d new value", candidatesCase2.length);
        if (candidatesCase2.length != 1)
            printf("s");
            printf(" obtained. \n");
    }else{
        printf("No new values. Going on.\n");
    }

    matrixInt322DStruct newXORDVvalues;
    *candidates = unionVec(&candidatesCase1, &candidatesCase2, NULL);
    newXORDVvalues.rows = XORextracted_values->rows + candidates->length;
    newXORDVvalues.cols = 2;
    newXORDVvalues.data = calloc(XORextracted_values->rows, sizeof(int32_t*));
    for (i = 0; i < newXORDVvalues.rows; i++) {
        newXORDVvalues.data[i] = calloc(2, sizeof(int32_t));
    }
    index = 0;
    for (i = 0; i < XORextracted_values->rows; i++) {
        for (j = 0; j < 2; j++) {
            newXORDVvalues.data[i][j] = XORextracted_values->data[i][j];
        }
        index++;
    }
    for (i = 0; i < candidates->length; i++) {
        newXORDVvalues.data[index][0] = candidates->data[i];
        newXORDVvalues.data[index][1] = totalDVhistogram->data[candidates->data[i]-1][1];
        index++;
    }
    return newXORDVvalues;
}


#endif
