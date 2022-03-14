#include <omp.h>
#include <memory.h>
#include <time.h>
#include "misc.hpp"

inline int64_t getRealTimeClockUs ( void );
/*
{
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    return tm.tv_sec * 1000000 + tm.tv_nsec / 1000;
}
*/
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
void ParallelMultiexpBa<Curve>::processChunks ( BatchAcc &ba, uint32_t idChunk ) 
{
    uint32_t chunkValue;
    int64_t offset;
    
    // #pragma omp parallel for
    for(uint32_t i=0; i<n; i++) {
        if (g.isZero(bases[i])) continue;
        if ((chunkValue = getChunk(i, idChunk))) {
            ba.addPoint(chunkValue, bases[i]);
        }
    }
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::reduce ( BatchAcc &ba, uint32_t idChunk ) 
{
    uint32_t p2ref, ndiv2, ndiv2abs, chunkOffset;
    uint32_t nBits = bitsPerChunk;

//    batchAcc.dumpStats();
//    printf("== (%s:%d) == cntToAffine: %d\n", __FILE__, __LINE__,g.cntToAffine);
    
    while (nBits > 0) {
        ndiv2 = 1 << (nBits-1);

        for (uint32_t i = 1; i< ndiv2; i++) {

            if (ba.isZero(i + ndiv2)) continue;
            ba.add(i, i + ndiv2);
            ba.add(ndiv2, i + ndiv2);
        }
        ba.calculateOnlyOneLoop();        
        --nBits;
    }
    ba.calculate();        
//    batchAcc.dumpStats();
//    printf("== (%s:%d) == cntToAffine: %d\n", __FILE__, __LINE__,g.cntToAffine);

    nBits = bitsPerChunk;
    while (nBits > 0) {
        ndiv2 = 1 << (nBits-1);
        for (int i = 0; i < (nBits-1); ++i) {
            ba.add(ndiv2, ndiv2);            
            ba.calculate();
        }
        if (nBits == 1) {
            ba.add(chunkResultRef, 1);
            break;
        }
           
        ba.add(chunkResultRef, ndiv2);
        --nBits;
    }
    ba.calculate();
//    batchAcc.dumpStats();
//    printf("== (%s:%d) == cntToAffine: %d\n", __FILE__, __LINE__,g.cntToAffine);
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::multiexp(typename Curve::Point &r, const typename Curve::PointAffine *_bases, const uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) 
{
//     nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
// batchAcc.dumpStats();
//    printf("== (%s:%d) == cntToAffine: %d\n", __FILE__, __LINE__,g.cntToAffine);
    uint64_t t1, t2, t3, t4, t5, t6, t7;
    // t1 = getRealTimeClockUs();
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

    // ThreadLimit threadLimit(1);

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


    // t2 = getRealTimeClockUs();

    typename Curve::PointAffine chunkResults[nChunks];

    #pragma omp parallel for
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        BatchAccumulators<Curve> ba(g);
        ba.defineAccumulators(accsPerChunk);
        chunkResultRef = ba.defineAccumulators(1);
        ba.setup(n >> 1, n >> 4);

        processChunks(ba, idChunk);
        // t3 = getRealTimeClockUs();
        ba.calculate();
        // t4 = getRealTimeClockUs();
        reduce(ba, idChunk);
        // t5 = getRealTimeClockUs();
        ba.calculate();
        // t6 = getRealTimeClockUs();
        g.copy(chunkResults[idChunk], ba.getValue(chunkResultRef));
        // ba.dumpStats();
    }
    
    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r, r);
        g.add(r, r, chunkResults[j]);
    }
    // t7 = getRealTimeClockUs();
    // printf("setup:%'ld processChunks:%'ld calculate:%'ld reduce:%'ld calculate:%'ld final:%'ld\n",
    //    t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6);
    // batchAcc.dumpStats();
}
