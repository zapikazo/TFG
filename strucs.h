//
//  strucs.h
//
//  Created on 5/7/17.
//  Copyright © 2017. All rights reserved.
//

#ifndef strucs_h
#define strucs_h

typedef struct {
    int32_t** pattern;
    int patternLength;
    int nBitsAddress;
    int nWordWidth;
    int nBits4Blocks;
    int nRoundsInPattern;
} programInfo;

typedef struct {
    int32_t** data;
    int rows;
    int cols;
} matrixInt322DStruct;

typedef struct {
    int32_t*** data;
    int rows;
    int cols;
    int dims;
} matrixInt323DStruct;

typedef struct {
    uint32_t** data;
    int rows;
    int cols;
} matrixUint322DStruct;

typedef struct {
    uint32_t*** data;
    int rows;
    int cols;
    int dims;
} matrixUint323DStruct;

typedef struct {
    int* data;
    int length;
} vectorIntStruct;

typedef struct {
    int32_t* data;
    int length;
} vectorInt32Struct;

typedef struct {
    uint32_t* data;
    int length;
} vectorUint32Struct;

typedef struct {
    bool** data;
    int rows;
    int cols;
} matrixBool2DStruct;

typedef struct {
    vectorInt32Struct XORextracted_values;
    vectorInt32Struct POSextracted_values;
} previouslyKnownValues;

#endif /* strucs_h */
