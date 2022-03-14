#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <stack>

#include "batch_accumulators.hpp"

// #define ZERO(X) g.copy(X, g.zeroAffine())
// #define COPY(D,S) g.copy(D,S)
#define IS_ZERO(X) g.isZero(X)

#define ZERO(X) X=zero
#define COPY(D,S) D=S
// #define IS_ZERO(X) memcmp(&X,&zero, sizeof(zero)) == 0

template <typename Curve>
BatchAccumulators<Curve>::BatchAccumulators (Curve &_g)
    :g(_g)
{
    valuesCount = 0;
    valuesSize = 0;
    accumulators = NULL;
    accumulatorsCount = 0;
    leftValues = rightValues = resultValues = NULL;
    accumulatorIds = NULL;
    currentLoop = 0;
    g.copy(zero, g.zeroAffine());
    clearStats();
}

template <typename Curve>
int64_t BatchAccumulators<Curve>::defineAccumulators ( int64_t count )
{
    int64_t _reference = accumulatorsCount;
    accumulatorsCount += count;
    return _reference;
}

template <typename Curve>
void BatchAccumulators<Curve>::setup (int64_t _initialValues, int64_t _deltaValues)
{
    freeValues();
    initialValues = _initialValues;
    deltaValues = _deltaValues;

    valuesSize = _initialValues;
    accumulators = (Accumulator *)malloc(accumulatorsCount * sizeof(accumulators[0]));
    clear();

    leftValues = (typename Curve::PointAffine *)malloc(valuesSize * sizeof(leftValues[0]));
    rightValues = (typename Curve::PointAffine *)malloc(valuesSize * sizeof(rightValues[0]));
    resultValues = (typename Curve::PointAffine *)malloc(valuesSize * sizeof(resultValues[0]));
    accumulatorIds = (int64_t *)calloc(valuesSize, sizeof(accumulatorIds[0]));
}

template <typename Curve>
void BatchAccumulators<Curve>::clear (void)
{
    // #pragma omp parallel for
    for (int64_t index = 0; index < accumulatorsCount; ++index) {
        Accumulator &accumulator = accumulators[index];
        accumulator.ready = true;
        accumulator.lastLoop = 0;
        ZERO(accumulator.value);
        ZERO(accumulator.singleValue);
    }    

    valuesCount = 0;
}

template <typename Curve>
bool BatchAccumulators<Curve>::isZero (int64_t accumulatorId)
{
    // assert(accumulatorId >=0 && accumulatorId < accumulatorsCount);
    return IS_ZERO(accumulators[accumulatorId].value);
}


template <typename Curve>
void BatchAccumulators<Curve>::freeValues (void)
{
    valuesCount = 0;
    valuesSize = 0;

    if (accumulators) {
        free(accumulators);
        accumulators = NULL;
    }

    if (leftValues) {
        free(leftValues);
        leftValues = NULL;
    }
    if (rightValues) {
        free(rightValues);
        rightValues = NULL;
    }
    if (resultValues) {
        free(resultValues);
        resultValues = NULL;
    }

    if (accumulatorIds) {
        free(accumulatorIds);
        accumulatorIds = NULL;
    }
}

template <typename Curve>
void BatchAccumulators<Curve>::dbl ( int64_t accumulatorId )
{
    addPoint(accumulatorId, getValue(accumulatorId), false);    
}

template <typename Curve>
void BatchAccumulators<Curve>::add ( int64_t accumulatorId, int64_t valueAccumulatorId )
{
    // assert(valueAccumulatorId >= 0 && valueAccumulatorId < accumulatorsCount);
    // assert(accumulators[valueAccumulatorId].ready);
    internalAdd(accumulatorId, accumulators[valueAccumulatorId].value, false);    
}

template <typename Curve>
void BatchAccumulators<Curve>::addPoint ( int64_t accumulatorId, const typename Curve::PointAffine &value )
{
    internalAdd(accumulatorId, value, false);
}

template <typename Curve>
int64_t BatchAccumulators<Curve>::incValuesCount ( void )
{
    // printf("incValuesCount(valuesCount:%ld)\n", valuesCount);
    ++valuesCount;    
    if (valuesCount >= valuesSize) {
        printf("RESIZE!!!\n");
        resize();
    }

    if (valuesCount > maxValues) {
        maxValues = valuesCount;
    }
    return valuesCount - 1;
}

template <typename Curve>
void BatchAccumulators<Curve>::internalAdd ( int64_t accumulatorId, const typename Curve::PointAffine &value, bool internal )
{
    // assert(accumulatorId >= 0 && accumulatorId < accumulatorsCount);
    if (IS_ZERO(value)) {
        return;
    }

    Accumulator &accumulator = accumulators[accumulatorId];
   
    // TEST
    // g.add(accumulator.value, accumulator.value, value);

    #ifdef __DEBUG__
    std::cout << "\x1b[1;33m/** BEFORE ** VALUE " << g.toString(value) << " acc[" << accumulatorId << "] ready:" << accumulator.ready << " value:" << g.toString(accumulator.value) 
              << " singleValue:" << g.toString(accumulator.singleValue) 
              << " internal:" << internal
              << "\x1b[0m\n";
    #endif
    if (!internal) {
        if (accumulator.ready) {
            if (!IS_ZERO(accumulator.value)) {
                int64_t index = incValuesCount();
                accumulator.ready = false;
            // assert(accumulator.index >= 0 && accumulator.index < valuesSize);
                COPY(leftValues[index], accumulator.value);
                COPY(rightValues[index], value);
                COPY(accumulator.value, zero);
                accumulatorIds[index] = accumulatorId;
            }
            else {
                COPY(accumulator.value, value);
            }
            #ifdef __DEBUG__
            std::cout << "\x1b[33m\\** AFTER1 ** acc[" << accumulatorId << "] ready:" << accumulator.ready << " value:" << g.toString(accumulator.value) 
              << " singleValue:" << g.toString(accumulator.singleValue) 
              << " internal:" << internal
              << "\x1b[0m\n";
            #endif
            return;
        }
    }
    accumulator.ready = false;
    
    // If is there a single value, and on right to complete the pair.
    if (IS_ZERO(accumulator.singleValue)) {
        COPY(accumulator.singleValue, value);
        #ifdef __DEBUG__
        std::cout << "\x1b[33m\\** AFTER2 ** acc[" << accumulatorId << "] ready:" << accumulator.ready << " value:" << g.toString(accumulator.value) 
            << " singleValue:" << g.toString(accumulator.singleValue) 
            << " internal:" << internal
            << "\x1b[0m\n";
        #endif
        return;
    }

    int64_t index = incValuesCount();
    COPY(leftValues[index], accumulator.singleValue);
    COPY(rightValues[index], value);
    ZERO(accumulator.singleValue);
    accumulatorIds[index] = accumulatorId;
    #ifdef __DEBUG__
    std::cout << "\x1b[33m\\** AFTER3 ** acc[" << accumulatorId << "] ready:" << accumulator.ready << " value:" << g.toString(accumulator.value) 
        << " singleValue:" << g.toString(accumulator.singleValue) 
        << " internal:" << internal
        << "\x1b[0m\n";
    #endif
}

template <typename Curve>
const typename Curve::PointAffine &BatchAccumulators<Curve>::getValue ( int64_t accumulatorId )
{
    // assert(accumulatorId >= 0 && accumulatorId < accumulatorsCount);
    // assert(accumulators[accumulatorId].ready);

    return accumulators[accumulatorId].value;
}


template <typename Curve>
void BatchAccumulators<Curve>::resize ( void )
{
    valuesSize += deltaValues;
    leftValues = (typename Curve::PointAffine *)realloc(leftValues, valuesSize * sizeof(leftValues[0]));
    rightValues = (typename Curve::PointAffine *)realloc(rightValues, valuesSize * sizeof(rightValues[0]));
    resultValues = (typename Curve::PointAffine *)realloc(resultValues, valuesSize * sizeof(resultValues[0]));
    accumulatorIds = (int64_t *)realloc(accumulatorIds, valuesSize * sizeof(accumulatorIds[0]));
}


template <typename Curve>
bool BatchAccumulators<Curve>::calculateOnlyOneLoop ( void )
{
    if (!valuesCount) {
        return true;
    }

    ++currentLoop;
    multiAdd();

    int64_t oldValuesCount = valuesCount;
    valuesCount = 0;

    for (int64_t index = 0; index < oldValuesCount; ++index ) {
        int64_t accumulatorId = accumulatorIds[index];  
        Accumulator &accumulator = accumulators[accumulatorId];

        #ifdef __DEBUG__
        std::cout << "\x1b[34m** INDEX ** " << index << " RESULT:" << g.toString(resultValues[index]) << " acc[" << accumulatorId << "] ready:" << accumulator.ready << " value:" << g.toString(accumulator.value) 
                << " singleValue:" << g.toString(accumulator.singleValue)
                << " lastLoop:" << accumulator.lastLoop << " currentLoop:" << currentLoop
                << "\x1b[0m\n";
        #endif

        if (accumulator.lastLoop != currentLoop) {
            accumulator.lastLoop = currentLoop;
            if (!IS_ZERO(accumulator.singleValue)) {
                internalAdd(accumulatorId, resultValues[index], false);
                continue;
            }
            accumulator.ready = true;
            COPY(accumulator.value, resultValues[index]);
            continue;
        }
        else if (accumulator.ready) {
            COPY(accumulator.singleValue, accumulator.value);
            accumulator.ready = false;
        }
        internalAdd(accumulatorId, resultValues[index], true);
    }
    return (valuesCount == 0);
}

template <typename Curve>
void BatchAccumulators<Curve>::multiAdd ( void )
{
    #ifdef BATCH_ACCUMULATORS_STATS
    ++multiAddOperations;
    totalMultiAddValues += valuesCount;
    if (maxMultiAddValues < valuesCount) {
        maxMultiAddValues = valuesCount;
    }
    if (minMultiAddValues < 0 || minMultiAddValues > valuesCount) {
        minMultiAddValues = valuesCount;
    }
    #endif

    int valuesBlock = valuesCount > 32 ? valuesCount >> 3 : valuesCount;
    int nBlocks = valuesCount / valuesBlock;
    int valuesLastBlock = valuesCount - (nBlocks - 1) * valuesBlock;
    
//    dumpStats();
//    printf("== (%s:%d) == cntToAffine: %d\n", __FILE__, __LINE__,g.cntToAffine);

    // g.multiAdd(resultValues, leftValues, rightValues, valuesCount);
    
    #pragma omp parallel for
    for (int block = 0; block < nBlocks; ++ block) {
        int count = (block == (nBlocks - 1)) ? valuesLastBlock : valuesBlock;
        int offset = block * valuesBlock;
        g.multiAdd(resultValues + offset, leftValues + offset, rightValues + offset, count);
    }
}

template <typename Curve>
void BatchAccumulators<Curve>::dumpStats ( void )
{
    printf("maxValues:%'ld\n",maxValues);
    printf("multiAddOperations:%'ld\n",multiAddOperations);
    printf("maxMultiAddValues:%'ld\n",maxMultiAddValues);
    printf("minMultiAddValues:%'ld\n",minMultiAddValues);
    printf("totalMultiAddValues:%'ld\n",totalMultiAddValues);
    printf("avgMultiAddValues:%'.02f\n",(double)totalMultiAddValues/(double)multiAddOperations);
}

template <typename Curve>
void BatchAccumulators<Curve>::clearStats ( void )
{
    maxValues=0;
    multiAddOperations=0;
    maxMultiAddValues=0;
    minMultiAddValues=-1;
    totalMultiAddValues=0;
}