#include <omp.h>
#include <memory.h>
#include "misc.hpp"

//    auto _O1 = O1; 
//    auto _O2 = O2; 
//    printf("\n## %d: %s\n## %s\n +%s\n = %s\n", __LINE__, #R "=" #O1 "+" #O2,g.toString(_O1).c_str(), g.toString(_O2).c_str(), g.toString(R).c_str()); 

#define __G_ADD__TREE__(R,O1,O2) \
{ \
    g.add(R, O1, O2); \
} 

template <typename Curve>
uint32_t ParallelMultiexpBa<Curve>::getChunk(uint32_t scalarIdx, uint32_t chunkIdx) 
{
    uint32_t bitStart = chunkIdx*bitsPerChunk;
    uint32_t byteStart = bitStart/8;
    uint32_t efectiveBitsPerChunk = bitsPerChunk;
    if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
    if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;
    uint32_t shift = bitStart - byteStart*8;
    uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + byteStart);
    v = v >> shift;
    v = v & ( (1 << efectiveBitsPerChunk) - 1);
    return uint32_t(v);
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::processChunks ( void ) 
{
    uint32_t chunkValue;

    for(uint32_t i=0; i<n; i++) {
        if (g.isZero(bases[i])) continue;
        for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
            if (chunkValue = getChunk(i, idChunk)) {
                batchAcc.add(idChunk * accsPerChunk + chunkValue, bases[i]);
            }
        }
    }
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::reduce ( void ) 
{
    uint32_t p2ref, ndiv2, ndiv2abs, chunkOffset;
    uint32_t nBits = bitsPerChunk;

    while (nBits > 0) {
        ndiv2 = 1 << (nBits-1);
        for (uint32_t idChunk = 0; idChunk < nChunks; idChunk++) {
            chunkOffset = idChunk * accsPerChunk;
            ndiv2abs = chunkOffset + ndiv2;

            for (uint32_t i = 1; i< ndiv2; i++) {
                uint32_t index = chunkOffset + i;

                if (batchAcc.isZero(index + ndiv2)) continue;
                batchAcc.add(index, index + ndiv2);
                batchAcc.add(ndiv2abs, index + ndiv2);
            }
        }
        batchAcc.calculateOnlyOneLoop();        
        --nBits;
    }
    batchAcc.calculate();        

    for (uint32_t idChunk = 0; idChunk < nChunks; idChunk++) {
        nBits = bitsPerChunk;
        while (nBits > 0) {
            ndiv2 = 1 << (nBits-1);
            auto value = batchAcc.getValue(idChunk * accsPerChunk + ndiv2);
            printf("==> BA-CHUNK-NDIV2 A = 2**%d * %s\n", nBits-1, g.toString(value).c_str());
            --nBits;
        }
    }

    batchAcc.calculate();
  
    nBits = bitsPerChunk;
    while (nBits > 0) {
        ndiv2 = 1 << (nBits-1);
        for (int i = 0; i < (nBits-1); ++i) {
            for (uint32_t idChunk = 0; idChunk < nChunks; idChunk++) {
                ndiv2abs = idChunk * accsPerChunk + ndiv2;
                batchAcc.add(ndiv2abs, ndiv2abs);            
            }
            batchAcc.calculate();
        }
        for (uint32_t idChunk = 0; idChunk < nChunks; idChunk++) {
            ndiv2abs = idChunk * accsPerChunk + ndiv2;
            if (nBits == 1) {
                chunkOffset = idChunk * accsPerChunk;
                batchAcc.add(chunkResultRef + idChunk, chunkOffset + 1);
                continue;
            }
           
            if (idChunk == 0) {
                auto value = batchAcc.getValue(ndiv2abs);
                auto value2 = batchAcc.getValue(chunkResultRef + idChunk);
                printf("ndiv2abs:%d\n", ndiv2abs);
                printf("==> BA-CHUNK-ADD-RESULT[%d] = %s + %s\n", nBits, g.toString(value2).c_str(), g.toString(value).c_str());    
            }
            batchAcc.add(chunkResultRef + idChunk, ndiv2abs);
        }
        --nBits;
    }
    batchAcc.calculate();

    for (uint32_t idChunk = 0; idChunk < nChunks; idChunk++) {
        nBits = bitsPerChunk;
        while (nBits > 0) {
            ndiv2 = 1 << (nBits-1);
            auto value = batchAcc.getValue(idChunk * accsPerChunk + ndiv2);
            printf("==> BA-CHUNK-NDIV2 B = 2**%d * %s\n", nBits-1, g.toString(value).c_str());
            --nBits;
        }
    }

}

template <typename Curve>
void ParallelMultiexpBa<Curve>::multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) 
{
//     nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

//    ThreadLimit threadLimit (nThreads);

    if (n==0) {
        g.copy(r, g.zero());
        return;
    }
    if (n==1) {
        g.mulByScalar(r, bases[0], scalars, scalarSize);
        return;
    }

    bitsPerChunk = log2(n / PME2_PACK_FACTOR);
    if (bitsPerChunk > PME2_MAX_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MAX_CHUNK_SIZE_BITS;
    if (bitsPerChunk < PME2_MIN_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MIN_CHUNK_SIZE_BITS;
    nChunks = ((scalarSize*8 - 1 ) / bitsPerChunk)+1;
    accsPerChunk = 1 << bitsPerChunk;  // In the chunks last bit is always zero.

    batchAcc.defineAccumulators(nChunks * accsPerChunk);
    
    chunkResultRef = batchAcc.defineAccumulators(nChunks);
    resultRef= batchAcc.defineAccumulators(1);
    batchAcc.defineAccumulators(1000);

    batchAcc.setup(100000, 100000);
    
    processChunks();
    batchAcc.calculate();
    reduce();
    batchAcc.calculate();

    typename Curve::PointAffine value;
    value = batchAcc.getValue(chunkResultRef + nChunks - 1);
    g.copy(r, value);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r, r);
        value = batchAcc.getValue(chunkResultRef + j);
        g.add(r, r, value);
    }
    batchAcc.printStats();
}
