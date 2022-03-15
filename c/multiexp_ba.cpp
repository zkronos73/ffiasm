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
void ParallelMultiexpBa<Curve>::prepareGetChunk ( void ) 
{
    chunkInfo = new ChunkInfo[nChunks];
    
    for (int idChunk = 0; idChunk < nChunks; ++idChunk) {
        uint32_t bitStart = idChunk * bitsPerChunk;
        uint32_t byteStart = bitStart/8;
        uint32_t efectiveBitsPerChunk = bitsPerChunk;
        if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
        if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;
        chunkInfo[idChunk].byteStart = byteStart;
        chunkInfo[idChunk].mask = ((1 << efectiveBitsPerChunk) - 1); 
        chunkInfo[idChunk].shift = bitStart - byteStart*8;
    }
}

template <typename Curve>
uint32_t ParallelMultiexpBa<Curve>::fastGetChunk(uint32_t scalarIdx, uint32_t idChunk) 
{
    uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + chunkInfo[idChunk].byteStart);
    v = v >> chunkInfo[idChunk].shift;
    v = v & chunkInfo[idChunk].mask;
    return uint32_t(v);
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
void ParallelMultiexpBa<Curve>::processChunks ( BatchAcc &ba, uint32_t idChunk ) 
{
    uint32_t chunkValue;
//    int64_t offset;
    
//    uint32_t chunkBlocks = accsPerChunk >> 14;

//    #pragma omp parallel for
//    for (uint32_t j=0; j < chunkBlocks; ++j) {
        for (uint32_t i=0; i<n; i++) {
            if (g.isZero(bases[i])) continue;
            chunkValue = fastGetChunk(i, idChunk);
            if (!chunkValue) continue;
            // if (!chunkValue || (chunkValue >> 14) != j) continue;
            ba.addPoint(chunkValue, bases[i]);
        }
//    }
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
    bitsPerChunk = 14;

    nChunks = ((scalarSize*8 - 1 ) / bitsPerChunk)+1;
    accsPerChunk = 1 << bitsPerChunk;  // In the chunks last bit is always zero.

    printf("chunks:%d accsPerChunk:%'ld\n", nChunks, accsPerChunk);

    // t2 = getRealTimeClockUs();

    prepareGetChunk();
    typename Curve::PointAffine chunkResults[nChunks];

#ifdef __FULL_STATS__
    uint64_t t0,t1,t2;
    uint64_t t[5][nChunks];
    BatchAccumulatorsStats stats[nChunks];
#endif
/*
    int32_t sum1 = 0;
    int32_t sum2 = 0;
    t0 = getRealTimeClockUs();
    prepareGetChunk();
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        for (uint32_t i = 0; i < n; ++i ) {
            sum1 += fastGetChunk(i, idChunk);
        }
    }
    t1 = getRealTimeClockUs();
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        for (uint32_t i = 0; i < n; ++i ) {
            sum2 += getChunk(i, idChunk);
        }
    }
    t2 = getRealTimeClockUs();
    printf("sum1:%d sum2:%d\n", sum1, sum2);
    printf("T1:%'ld T2:%'ld\n", t1-t0, t2-t1);
    exit(EXIT_FAILURE);
*/
#ifdef __FULL_STATS__
    t0 = getRealTimeClockUs();
#endif
    #pragma omp parallel for
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        BatchAccumulators<Curve> ba(g);
        ba.defineAccumulators(accsPerChunk);
        chunkResultRef = ba.defineAccumulators(1);
        ba.setup(n >> 1, n >> 4);

#ifdef __FULL_STATS__    
        t[0][idChunk] = getRealTimeClockUs();
#endif
        processChunks(ba, idChunk);
#ifdef __FULL_STATS__    
        t[1][idChunk] = getRealTimeClockUs();
#endif
        ba.calculate();
#ifdef __FULL_STATS__    
        t[2][idChunk] = getRealTimeClockUs();
#endif
        reduce(ba, idChunk);
#ifdef __FULL_STATS__    
        t[3][idChunk] = getRealTimeClockUs();
#endif
        g.copy(chunkResults[idChunk], ba.getValue(chunkResultRef));
#ifdef __FULL_STATS__    
        t[4][idChunk] = getRealTimeClockUs();
        stats[idChunk] = ba.stats;
#endif
    }
#ifdef __FULL_STATS__    
    t1 = getRealTimeClockUs();
#endif
    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r, r);
        g.add(r, r, chunkResults[j]);
    }

#ifdef __FULL_STATS__    
    t2 = getRealTimeClockUs();

    double factor = (t2 - t0) > 2000000 ? 1000000.0: 1000.0;

    if (factor == 1000.0) {
        printf("NOTE: time in milliseconds\n");
    }
    if (factor == 1000000.0) {
        printf("NOTE: time in seconds\n");        
    }
    printf("n:%'d threads-max:%d chunks:%d accsPerChunk:%'ld\n", n, omp_get_max_threads(), nChunks, accsPerChunk);
    printf("iCh|total   |process |calculat|reduce  |copy    |start   |wait    |MaxValues   |multC|singC|shortC|minB|maxBlocks|minA|maxAdds     |totalAdds   |resizes\n");
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        printf("#%2d|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f|%'12ld|%5ld|%5ld|%6ld|%4ld|%9ld|%4ld|%'12ld|%'12ld|%7d\n", idChunk, 
            ((double)(t[4][idChunk]-t[0][idChunk]))/factor,
            ((double)(t[1][idChunk]-t[0][idChunk]))/factor,
            ((double)(t[2][idChunk]-t[1][idChunk]))/factor,
            ((double)(t[3][idChunk]-t[2][idChunk]))/factor,
            ((double)(t[4][idChunk]-t[3][idChunk]))/factor,
            ((double)(t[0][idChunk]-t0))/factor,
            ((double)(t1-t[4][idChunk]))/factor,
            stats[idChunk].maxValues,
            stats[idChunk].multiAddOperations,
            stats[idChunk].singleMultiAdds,
            stats[idChunk].shortMultiAdds,
            stats[idChunk].minBlock,
            stats[idChunk].maxBlock,
            stats[idChunk].minMultiAddValues,
            stats[idChunk].maxMultiAddValues,
            stats[idChunk].totalMultiAddValues,
            stats[idChunk].resizes         
            );
    }
    printf("final: %9.04f\n", ((double)(t2-t1))/factor);
#endif
    // t7 = getRealTimeClockUs();
    // printf("setup:%'ld processChunks:%'ld calculate:%'ld reduce:%'ld calculate:%'ld final:%'ld\n",
    //    t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6);
    // batchAcc.dumpStats();
}
