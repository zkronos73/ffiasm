#include <omp.h>
#include <memory.h>
#include "misc.hpp"

// #define __LEGACY__
#define __TRACE__
#ifdef __TRACE__
#define G_COPY(D, S) gCopy(__LINE__, D, S)
#define G_DBL(D, S) gDbl(__LINE__, D, S)
#define G_ADD(D, O1, O2) gAdd(__LINE__, D, O1, O2)
#define G_IS_ZERO(S) gIsZero(__LINE__, S)
#define G_MUL_BY_SCALAR(D, O, SD, SS) gMulByScalar(__LINE__, D, O, SD, SS)
#else
#define G_COPY(D, S) g.copy(D, S)
#define G_DBL(D, S) g.dbl(D, S)
#define G_ADD(D, O1, O2) g.add(D, O1, O2)
#define G_IS_ZERO(S) g.isZero(S)
#define G_MUL_BY_SCALAR(D, O, SD, SS) g.mulByScalar(D, O, SD, SS)
#endif
/*
template <typename Curve>
void ParallelMultiexp<Curve>::initAccs() {
    #pragma omp parallel for
    for (int i=0; i<nThreads; i++) {
        memset((void *)&(accs[i*accsPerChunk]), 0, accsPerChunk*sizeof(PaddedPoint));
    }
}
*/

template <typename Curve>
void ParallelMultiexp<Curve>::initAccs() {
    for (int i=0; i < accsPerChunk; i++) {
        operations.copy(accsReference + i, 0);
        // G_COPY(accs[i].p, g.zero());
    }
}

template <typename Curve>
uint32_t ParallelMultiexp<Curve>::getChunk(uint32_t scalarIdx, uint32_t chunkIdx) {
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
void ParallelMultiexp<Curve>::processChunk(uint32_t idChunk) {
    for(uint32_t i=0; i<n; i++) {
        if (g.isZero(bases[i])) continue;
        uint32_t chunkValue = getChunk(i, idChunk);
        if (chunkValue) {
            // G_ADD(accs[idThread*accsPerChunk+chunkValue].p, accs[idThread*accsPerChunk+chunkValue].p, bases[i]);
            operations.add(accsReference + chunkValue, accsReference + chunkValue, basesReference+i, __FILE__, __LINE__);
        }
    }
}


template <typename Curve>
void ParallelMultiexp<Curve>::reduce(int64_t chunkIndex, uint32_t nBits) 
{
    operations.copy(chunkResultsReference + chunkIndex, 0 );
    uint32_t ndiv2;
    
    while (true) {
        if (nBits == 1) {
            operations.add(chunkResultsReference + chunkIndex, chunkResultsReference + chunkIndex, accsReference + 1, __FILE__, __LINE__);
            operations.copy(accsReference + 1, 0 );
            return;
        }
        ndiv2 = 1 << (nBits - 1);

        for (uint32_t i = 1; i<ndiv2; i++) {
//            if (!G_IS_ZERO(accs[ndiv2 + i].p)) {
//            maybe, check if zero reference    
            operations.add(accsReference + i, accsReference + i, accsReference + ndiv2 + i, __FILE__, __LINE__);
            operations.add(accsReference + ndiv2, accsReference + ndiv2, accsReference + ndiv2 + i, __FILE__, __LINE__);
            operations.copy(accsReference + ndiv2 + i, 0 );
//            }
        }
        for (u_int32_t i=0; i<(nBits - 1); i++) {
            operations.dbl(accsReference + ndiv2, accsReference + ndiv2, __FILE__, __LINE__);
        }        
        operations.add(chunkResultsReference + chunkIndex, chunkResultsReference + chunkIndex, accsReference + ndiv2, __FILE__, __LINE__);
        operations.copy(accsReference + ndiv2, 0 );
        --nBits;
    }
}


template <typename Curve>
void ParallelMultiexp<Curve>::multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n ) {
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

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

    basesReference = operations.define("bases", n);
    accsReference = operations.define("accs", accsPerChunk);
    chunkResultsReference = operations.define("chunkResults", nChunks);

    for (int index = 0; index < n; ++index) {
        operations.set(basesReference + index, _bases[index]);
//        printf("bases[%d] @%ld %s\n", index, basesReference + index, g.toString(_bases[index]).c_str());
    }

    chunkResults = new typename Curve::Point[nChunks];
    accs = new PaddedPoint[accsPerChunk];

    printf("accsReference:%ld chunkResultsReference:%ld\n", accsReference,chunkResultsReference);

    std::cout << "InitTrees " << "\n"; 
    initAccs();
/*
    for (int index = 0; index < n; ++index) {
        typename Curve::PointAffine p1 = operations.resolve(accsReference + index);
        printf("PRE accs[%d].p = %s\n", index, g.toString(p1).c_str());
    }
*/
    for (uint32_t i=0; i<nChunks; i++) {
        processChunk(i);
    /*    for (int index = 0; index < n; ++index) {
            typename Curve::PointAffine p1 = operations.resolve(accsReference + index);
            printf("accs[%d].p = %s\n", index, g.toString(p1).c_str());
        }
        return;*/
        reduce(i, bitsPerChunk);
    }    
    delete[] accs;
    cache = false;

    BatchOperations::Reference resultReference = operations.define("result", 1);
    operations.copy(resultReference, chunkResultsReference + nChunks - 1);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) operations.dbl(resultReference, resultReference, __FILE__, __LINE__);
        operations.add(resultReference, resultReference, chunkResultsReference + j, __FILE__, __LINE__);
    }

    delete[] chunkResults; 
    printf("=== add:%ld dbl:%ld\n", operations.addCounter, operations.dblCounter);
    typename Curve::PointAffine result = operations.resolve(resultReference, 1000000000);
    printf("=== CALLS add:%ld dbl:%ld\n", operations.addCallCounter, operations.dblCallCounter);
    g.copy(r, result);    
}
