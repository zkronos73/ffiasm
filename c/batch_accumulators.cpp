#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <stack>

#include "batch_accumulators.hpp"

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
    for (int64_t index = 0; index < accumulatorsCount; ++index) {
        Accumulator &accumulator = accumulators[index];
        accumulator.index = -1;
        accumulator.ready = true;
        accumulator.count = 0;
        g.copy(accumulator.value, g.zero());
    }    

    valuesCount = 0;
}

template <typename Curve>
bool BatchAccumulators<Curve>::isZero (int64_t accumulatorId)
{
    // assert(accumulatorId >=0 && accumulatorId < accumulatorsCount);
    return g.isZero(accumulators[accumulatorId].value);
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
    auto value = getValue(accumulatorId);
    add(accumulatorId, value, false);    
}

template <typename Curve>
void BatchAccumulators<Curve>::add ( int64_t accumulatorId, int64_t valueAccumulatorId )
{
    // assert(valueAccumulatorId >= 0 && valueAccumulatorId < accumulatorsCount);
    // assert(accumulators[valueAccumulatorId].ready);
    internalAdd(accumulatorId, accumulators[valueAccumulatorId].value, false);    
}

template <typename Curve>
void BatchAccumulators<Curve>::add ( int64_t accumulatorId, typename Curve::PointAffine &value )
{
    internalAdd(accumulatorId, value, false);
}

template <typename Curve>
void BatchAccumulators<Curve>::internalAdd ( int64_t accumulatorId, typename Curve::PointAffine &value, bool internal )
{
    // assert(accumulatorId >= 0 && accumulatorId < accumulatorsCount);
    if (g.isZero(value)) {
        return;
    }

    Accumulator &accumulator = accumulators[accumulatorId];
   
    // TEST
    // g.add(accumulator.value, accumulator.value, value);

    accumulator.ready = false;

    if (!internal) {
        if ((valuesCount + 1) >= valuesSize) {
            resize();
        }
        if (!accumulator.count) {
            accumulator.index = valuesCount++;
            // assert(accumulator.index >= 0 && accumulator.index < valuesSize);
            g.copy(leftValues[accumulator.index], accumulator.value);
            accumulatorIds[accumulator.index] = accumulatorId;
            ++accumulator.count;
        }
    }
    
    // If is there a single value, and on right to complete the pair.
    if (accumulator.index >= 0) {
        // assert(accumulator.index >= 0 && accumulator.index < valuesSize);
        g.copy(rightValues[accumulator.index], value);
        accumulator.index = -1;
        return;
    }

    if (!internal) {
        ++accumulator.count;
    }

    accumulator.index = valuesCount++;
    if (valuesCount > maxValues) {
        maxValues = valuesCount;
    }
    // assert(accumulator.index >= 0 && accumulator.index < valuesSize);
    g.copy(leftValues[accumulator.index], value);
    accumulatorIds[accumulator.index] = accumulatorId;
}


template <typename Curve>
void BatchAccumulators<Curve>::prepareSingleValues ( void )
{
    #pragma omp parallel for
    for (int64_t index = 0; index < accumulatorsCount; ++index) {
        if (accumulators[index].index < 0) continue;
        // assert(accumulators[index].index >= 0 && accumulators[index].index < valuesSize);
        g.copy(rightValues[accumulators[index].index], g.zero());
    }
}

template <typename Curve>
typename Curve::PointAffine BatchAccumulators<Curve>::getValue ( int64_t accumulatorId )
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
    prepareSingleValues();

    multiAdd();

    for (int64_t index = 0; index < valuesCount; ++index ) {
        int64_t accumulatorId = accumulatorIds[index];
        Accumulator &accumulator = accumulators[accumulatorId];

        if (accumulator.count > 1) continue;
        
        g.copy(accumulator.value, resultValues[index]);
        accumulator.ready = true;
    }

    prepareAccumulatorsToNextRound();
    
    int64_t oldValuesCount = valuesCount;
    valuesCount = 0;

    for (int64_t index = 0; index < oldValuesCount; ++index ) {
        int64_t accumulatorId = accumulatorIds[index];
        Accumulator &accumulator = accumulators[accumulatorId];

        if (!accumulator.count) {            
            continue;            
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

    #pragma omp parallel for
    for (int block = 0; block < nBlocks; ++ block) {
        int count = (block == (nBlocks - 1)) ? valuesLastBlock : valuesBlock;
        int offset = block * valuesBlock;
        g.multiAdd(resultValues + offset, leftValues + offset, rightValues + offset, count);
    }
}

template <typename Curve>
void BatchAccumulators<Curve>::prepareAccumulatorsToNextRound ( void )
{
    #pragma omp parallel for
    for (int64_t index = 0; index < accumulatorsCount; ++index ) {
        Accumulator &accumulator = accumulators[index];
        accumulator.index = -1;
        if (accumulator.count > 1) {            
            // divide by 2, but if it's odd add one, as up division            
            accumulator.count = (accumulator.count >> 1) + (accumulator.count & 0x01);
        }
        else {
            accumulator.count = 0;
        }
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