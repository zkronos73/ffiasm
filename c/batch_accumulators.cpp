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
BatchAccumulators<Curve>::~BatchAccumulators ( void )
{
    freeValues();
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
    initialValues = _initialValues < 16 ? 16:_initialValues;
    deltaValues = _deltaValues < 16 ? 16:_deltaValues;

    valuesSize = _initialValues;
//    printf("[%p] SETUP valuesSize:%ld\n", this, valuesSize);
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
int64_t BatchAccumulators<Curve>::incValuesCount ( void )
{
    // printf("incValuesCount(valuesCount:%ld)\n", valuesCount);
    int64_t result;
    // #pragma omp critical
    {
        result = valuesCount;
        ++valuesCount;    
        if (valuesCount >= valuesSize) {
            ++stats.resizes;
            resize();
        }
    }
    if (valuesCount > stats.maxValues) {
        stats.maxValues = valuesCount;
    }
    return result;
}

template <typename Curve>
void BatchAccumulators<Curve>::add ( int64_t accumulatorId, int64_t valueAccumulatorId )
{
    // assert(valueAccumulatorId >= 0 && valueAccumulatorId < accumulatorsCount);
    // assert(accumulators[valueAccumulatorId].ready);
    addPoint(accumulatorId, accumulators[valueAccumulatorId].value);
}

template <typename Curve>
void BatchAccumulators<Curve>::addPoint ( int64_t accumulatorId, const typename Curve::PointAffine &value )
{
    // printf("[%p] adding value (%ld, %s @%p)\n", this, accumulatorId, g.toString(value).c_str(), &value);
    if (!nonInternalAddBlock(accumulatorId, value)) {
        internalAddBlock(accumulatorId, value);
    }
}


template <typename Curve>
bool BatchAccumulators<Curve>::nonInternalAddBlock ( int64_t accumulatorId, const typename Curve::PointAffine &value )
{
    // assert(accumulatorId >= 0 && accumulatorId < accumulatorsCount);
    if (IS_ZERO(value)) {
        return true;
    }

    // Accumulator &accumulator = accumulators[accumulatorId];
   
    // TEST
    // g.add(accumulator.value, accumulator.value, value);

    if (!accumulators[accumulatorId].ready) {
        return false;
    }

    if (!IS_ZERO(accumulators[accumulatorId].value)) {
        int64_t index = incValuesCount();
        accumulators[accumulatorId].ready = false;
    // assert(accumulator.index >= 0 && accumulator.index < valuesSize);
        COPY(leftValues[index], accumulators[accumulatorId].value);
        COPY(rightValues[index], value);
        COPY(accumulators[accumulatorId].value, zero);
        accumulatorIds[index] = accumulatorId;
    }
    else {
        COPY(accumulators[accumulatorId].value, value);
    }
    return true;
}

template <typename Curve>
void BatchAccumulators<Curve>::internalAddBlock ( int64_t accumulatorId, const typename Curve::PointAffine &value )
{
    Accumulator &accumulator = accumulators[accumulatorId];
   
    accumulator.ready = false;
    
    // If is there a single value, and on right to complete the pair.
    if (IS_ZERO(accumulator.singleValue)) {
        COPY(accumulator.singleValue, value);
        return;
    }

    int64_t index = incValuesCount();
    COPY(leftValues[index], accumulator.singleValue);
    COPY(rightValues[index], value);
    ZERO(accumulator.singleValue);
    accumulatorIds[index] = accumulatorId;
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

        if (accumulator.lastLoop != currentLoop) {
            accumulator.lastLoop = currentLoop;
            if (!IS_ZERO(accumulator.singleValue)) {
                internalAddBlock(accumulatorId, resultValues[index]);
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
        internalAddBlock(accumulatorId, resultValues[index]);
    }
    return (valuesCount == 0);
}

template <typename Curve>
void BatchAccumulators<Curve>::multiAdd ( void )
{
    const int maxBlockSize = 16*1024;
    const int minBlockSize = 64;

    int valuesBlock = valuesCount;
    if (valuesCount > minBlockSize) {
//         valuesBlock = valuesCount / omp_get_max_threads();
        valuesBlock = valuesCount / 32;
        if (valuesBlock < minBlockSize) valuesBlock = minBlockSize;
        if (valuesBlock > maxBlockSize) valuesBlock = maxBlockSize;
    }

    int nBlocks = valuesCount / valuesBlock;
    int valuesLastBlock = valuesCount - (nBlocks - 1) * valuesBlock;

    #ifdef BATCH_ACCUMULATORS_STATS
    ++stats.multiAddOperations;
    stats.totalMultiAddValues += valuesCount;
    if (stats.maxMultiAddValues < valuesCount) {
        stats.maxMultiAddValues = valuesCount;
    }
    if (stats.minMultiAddValues < 0 || stats.minMultiAddValues > valuesCount) {
        stats.minMultiAddValues = valuesCount;
    }
    if (valuesCount == 1) {
        ++stats.singleMultiAdds;
    }
    if (valuesCount <= 8) {
        ++stats.shortMultiAdds;
    }
    if (stats.maxBlock < valuesBlock) {
        stats.maxBlock = valuesBlock;
    }
    if (stats.minBlock < 0 || stats.minBlock > valuesBlock) {
        stats.minBlock = valuesBlock;
    }
    #endif
    
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
    printf("maxValues:%'ld\n",stats.maxValues);
    printf("multiAddOperations:%'ld\n",stats.multiAddOperations);
    printf("maxMultiAddValues:%'ld\n",stats.maxMultiAddValues);
    printf("minMultiAddValues:%'ld\n",stats.minMultiAddValues);
    printf("totalMultiAddValues:%'ld\n",stats.totalMultiAddValues);
    printf("avgMultiAddValues:%'.02f\n",(double)stats.totalMultiAddValues/(double)stats.multiAddOperations);
    printf("resizes:%'d\n",stats.resizes);
}

template <typename Curve>
void BatchAccumulators<Curve>::clearStats ( void )
{
    memset(&stats, 0, sizeof(stats));
    stats.minMultiAddValues = -1;
    stats.minBlock = -1;
}

template <typename Curve>
void BatchAccumulators<Curve>::prepareToJoin ( uint64_t offset )
{
    // printf("offset:%ld valuesCount:%ld\n", offset, valuesCount);
    for (int64_t index = 0; index < valuesCount; ++index) {
        accumulatorIds[index] += offset;
    }
}

template <typename Curve>
void BatchAccumulators<Curve>::join ( BatchAccumulators *aux, uint64_t offset )
{
    // must call prepareToJoin
    if (currentLoop < aux->currentLoop) {
        currentLoop = aux->currentLoop;
    }
    memcpy(accumulators + offset, aux->accumulators, aux->accumulatorsCount * sizeof(accumulators[0]));
    if (aux->valuesCount) {
        int64_t newValuesCount = valuesCount + aux->valuesCount;
        while (newValuesCount >= valuesSize) {
            resize();
        }
        memcpy(leftValues + valuesCount, aux->leftValues, aux->valuesCount * sizeof(leftValues[0]));
        memcpy(rightValues + valuesCount, aux->rightValues, aux->valuesCount * sizeof(rightValues[0]));
        memcpy(resultValues + valuesCount, aux->resultValues, aux->valuesCount * sizeof(resultValues[0]));

        // prepareToJoin
        memcpy(accumulatorIds + valuesCount, aux->accumulatorIds, aux->valuesCount * sizeof(accumulatorIds[0]));
        printf("valuesCount:%'ld (+%'ld) valuesSize:%'ld (%'ld) ", valuesCount, aux->valuesCount, valuesSize, aux->valuesSize);
        valuesCount = newValuesCount;        
        printf(" => valuesCount:%'ld valuesSize:%'ld \n", valuesCount, valuesSize);
    }
}
