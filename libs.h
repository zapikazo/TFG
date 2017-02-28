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

typedef char tAddress[11]; //Has to be 9, 8+(\0)

char* getfield(FILE* file, char* taken){
    int c = 0, elems = 0;
    char* ch = malloc(sizeof(char));
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
	uint32_t* newColumnVector;
    newColumnVector = (uint32_t*)malloc(sizeof(uint32_t)*(*elems));
    int i;
    // Comparing last position of the array with 0xFFFFFFFF == 4294967295
    if(columnVector[*elems-1] == 0xFFFFFFFF){
        i = 0;
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
//NO SE COMO FUNCIONA ESTO
int* flipped_bits(tAddress word, tAddress pattern, int wordWidth){
    int* bitflips;
    bitflips = (int*)calloc(wordWidth, sizeof(int));
    int* result;
    result = calloc((wordWidth+1), sizeof(wordWidth)); //Array of zeros, != doc
    int nFound, i;
    nFound = 0;
    
    //We have to create the bitflips' array with xoring, skipping 0x and \0?
 /*   for(i=2; i<wordWidth-1; ++i)
        bitflips[i-2] = (char)(*word[i] ^ *pattern[i]);
    for(i=0; i<wordWidth-3; ++i)
        printf("%02X ", bitflips[i]);
    printf("\n");*/
    for(i=0; i<wordWidth-1; ++i){
        if(word[i] == pattern[i]){
            bitflips[i] = 0;
        } else {
            bitflips[i] = 1;
        }
    }
    for(i=0; i<wordWidth-1; ++i)
        printf("%02X ", bitflips[i]);
    printf("\n");
    
    
    for (i = 0; wordWidth-1; ++i) {
        if(bitflips[i] % 2 == 1){
            nFound++;
            result[nFound] = i;
        }
    }
    
    //Falta return
    return result;
}

//Funcion de relleno de fila superior y columna izquierda vector, elemens y matriz para modificarla


/*
 * Function to create the DV sets associated with the XOR and positive operations.
 * The idea is simple: The address vector is replicated to create a square matrix
 * traspose it, and operating. Every element in elements above the main
 * diagonal is an element of the DV. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 */
uint32_t** create_DVmatrix(uint32_t* addresses, int elems, char* op, uint32_t* RWcycles, int nbits4blocks, long int ln){
    int32_t **DVmatrix = (int32_t **)malloc((elems+2)*sizeof(int32_t *));
    int i, j, cont = 0, ntriu = 0;
    for (i = 0; i < elems+1; i++){
        DVmatrix[i] = (int32_t *)malloc((elems+2)*sizeof(int32_t));
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
uint32_t** copyOfMatrix(uint32_t** matrix, int elems){
    uint32_t** matrixbackup = (uint32_t **)malloc((elems+1)*sizeof(uint32_t *));
    int i, j;
    for(i = 0; i < elems+1; i++){
        for(j = 0; j < elems+1; j++){
            matrixbackup[i][j] = matrix[i][j];
        }
    }
    return matrixbackup;
}

/* Copia de vectores */
void copyOfVector(uint32_t* vector, uint32_t** vectorbackup, int ndv, int ktest){
    int i;
    for(i = 0; i < ndv+1; i++){
        vectorbackup[i][ktest-1] = vector[i];
    }
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
int32_t* create_histogram(uint32_t* vector, int ndv, int ln){
    int32_t* histogram = calloc(ln, sizeof(int32_t));
    int k, sum = 0;
    for (k = 0; k < ndv; k++) {
        sum = lookForElem(vector, ndv, k);
        histogram[k] = sum;
    }
    return histogram;
}

long long binomial(int n,int k)
{
    long long ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return ans;
}

/* This function just implements the theoretical expresions to determine expected
 * repetitions.  Variable names are meaningful. */
double ExpectedRepetitions(int m, int ndv, int LN, char* operation){
  double result = 0,logresulttmp = 0;

  if (strcmp(operation, "pos") == 0){
    //Difference of values
    if(abs(ndv/LN) < 0.1){
      logresulttmp =((1-m)*log(LN)+m*log(2)-log(m+1)+log(binomial((int)ndv,m)));
      result = exp(logresulttmp);
      //result = result*(1+(3*(m)-2*Int64(ndv))/(LN+1))*(1+2*(ndv-m)/(LN-2)/(m+2));
    }else{
      //pk=2/LN/(LN+1)*(LN+1-collect(1:1:LN));
      //result = sum(binomial(BigInt(ndv), m)*(pk.^m).*(1-pk).^(ndv-m));
    }
  }
  else if (strcmp(operation, "xor") == 0){
    // XOR
    //logresulttmp = Float64(log(binomial(BigInt(ndv),m))+(ndv-m)*log(LN-1)-(ndv-1)*log(LN));
    result = exp(logresulttmp);
  }else{
    //Badly specified function.
    printf("\n\tAn unknown operation to calculate expected repetitions. Returning a strange value\n");
    result = -1;
  }
  return result;
}






#endif
