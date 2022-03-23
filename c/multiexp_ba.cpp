#include <omp.h>
#include <memory.h>
#include <time.h>
#include "misc.hpp"

#ifdef __FULL_STATS__
#define __TIME_MARK(X) X=getRealTimeClockUs()
#else
#define __TIME_MARK(X)
#endif

inline int64_t getRealTimeClockUs ( void );

template <typename Curve>
void ParallelMultiexpBa<Curve>::prepareGetChunk ( void ) 
{
    freeChunkInfo();
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
ParallelMultiexpBa<Curve>::ParallelMultiexpBa ( Curve &_g )
    : g(_g), chunkInfo(NULL)
{
}

template <typename Curve>
ParallelMultiexpBa<Curve>::~ParallelMultiexpBa ( void )
{
    freeChunkInfo();
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::freeChunkInfo ( void )
{
    if (chunkInfo) {
        free(chunkInfo);
        chunkInfo = NULL;
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

#ifdef __DIVIDE_PROCESS_CHUNK__
template <typename Curve>
void ParallelMultiexpBa<Curve>::processChunks ( BatchAcc &ba, uint32_t idChunk ) 
{
    uint32_t chunkBlockBits = 2;
    uint32_t chunkBlocks = chunkBlockBits < 1 ? 1 : (1 << chunkBlockBits);
    uint32_t accsPerChunkBlock = accsPerChunk / chunkBlocks;
    BatchAccumulators<Curve> *cbba[chunkBlocks < 1 ? 1: chunkBlocks]; // chunkBlocks];

//    #pragma omp parallel for
    for (uint32_t j=0; j < chunkBlocks; ++j) {
        uint32_t shift = bitsPerChunk - chunkBlockBits;
        uint32_t accumulatorSize = ba.getValuesSize() / (chunkBlocks > 4 ? chunkBlocks/2 : 1);
        uint32_t chunkMask = (1 << shift) - 1;
        BatchAccumulators<Curve> *_ba;
        if (j) {
            cbba[j] = new BatchAccumulators<Curve>(g);        
            _ba = cbba[j];
            _ba->defineAccumulators(accsPerChunkBlock);
            _ba->setup(accumulatorSize, accumulatorSize / 4);
        }
        else {
            _ba = &ba;
        }

        uint32_t chunkValue;

        // divide in two parts to access shifted, but not improve performance.
        // TODO: analyze it's possible improve performace.

        uint32_t initial = (n / chunkBlocks) * j;
        for (uint32_t i=initial; i<n; i++) {

            chunkValue = fastGetChunk(i, idChunk);
            if (!chunkValue || (chunkValue >> shift) != j) continue;
            if (g.isZero(bases[i])) continue;   
            _ba->addPoint(chunkValue & chunkMask, bases[i]);

        }
        for (uint32_t i=0; i<initial; i++) {

            chunkValue = fastGetChunk(i, idChunk);
            if (!chunkValue || (chunkValue >> shift) != j) continue;
            if (g.isZero(bases[i])) continue;   

            _ba->addPoint(chunkValue & chunkMask, bases[i]);
        }

        _ba->calculate();
        if (j) {
            _ba->prepareToJoin(accsPerChunkBlock * j);
        }
    }
    
    for (uint32_t j=1; j < chunkBlocks; ++j) {
        ba.join(cbba[j], accsPerChunkBlock * j);
        delete cbba[j];
    }
}
#else
template <typename Curve>
void ParallelMultiexpBa<Curve>::processChunks ( BatchAcc &ba, uint32_t idChunk ) 
{
    for (uint32_t i=0; i<n; i++) {

        uint32_t chunkValue = fastGetChunk(i, idChunk);
        if (!chunkValue) continue;
        if (g.isZero(bases[i])) continue;
        ba.addPoint(chunkValue, bases[i]);
    }
}
#endif

template <typename Curve>
void ParallelMultiexpBa<Curve>::reduce ( BatchAcc &ba, uint32_t idChunk ) 
{
    uint32_t nBits = bitsPerChunk;

    while (nBits > 0) {
        uint32_t ndiv2 = 1 << (nBits-1);

        for (uint32_t i = 1; i< ndiv2; i++) {

            if (ba.isZero(i + ndiv2)) continue;
            ba.add(i, i + ndiv2);
            ba.add(ndiv2, i + ndiv2);
        }
        ba.calculateOnlyOneLoop();        
        --nBits;
    }
    ba.calculate();        

    nBits = bitsPerChunk;
    while (nBits > 0) {
        uint32_t ndiv2 = 1 << (nBits-1);
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
}

template <typename Curve>
void ParallelMultiexpBa<Curve>::multiexp(typename Curve::Point &r, const typename Curve::PointAffine *_bases, const uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) 
{
    // int nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

    // ThreadLimit threadLimit(nThreads);

    if (n==0) {
        g.copy(r, g.zero());
        return;
    }
    if (n==1) {
        g.mulByScalar(r, bases[0], scalars, scalarSize);
        return;
    }

    bitsPerChunk = log2(n / PME2_PACK_BA_FACTOR);
    if (bitsPerChunk > PME2_MAX_CHUNK_BA_SIZE_BITS) bitsPerChunk = PME2_MAX_CHUNK_BA_SIZE_BITS;
    if (bitsPerChunk < PME2_MIN_CHUNK_BA_SIZE_BITS) bitsPerChunk = PME2_MIN_CHUNK_BA_SIZE_BITS;

    nChunks = ((scalarSize*8 - 1 ) / bitsPerChunk)+1;
    accsPerChunk = 1 << bitsPerChunk;  // In the chunks last bit is always zero.

    printf("chunks:%d accsPerChunk:%'ld\n", nChunks, accsPerChunk);

    prepareGetChunk();
    typename Curve::PointAffine chunkResults[nChunks];

#ifdef __FULL_STATS__
    uint64_t t0,t1,t2;
    uint64_t t[5][nChunks];
    BatchAccumulatorsStats stats[2][nChunks];
    #warning FULL STATS ACTIVE !!
#endif

    __TIME_MARK(t0);
    // TODO: change
    // not nice, but for performance 
    chunkResultRef = accsPerChunk;

    int64_t size = n < 2048 ? n : n / 2;

    #pragma omp parallel for
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        BatchAccumulators<Curve> ba(g);
        ba.defineAccumulators(accsPerChunk);

        // inside for performance
        int32_t chunkResultRef = ba.defineAccumulators(1);
        ba.setup(size, size/2);

        __TIME_MARK(t[0][idChunk]);
        processChunks(ba, idChunk);

        __TIME_MARK(t[1][idChunk]);
        ba.calculate();

#ifdef __FULL_STATS__    
        stats[0][idChunk] = ba.stats;
        ba.clearStats();
#endif

        __TIME_MARK(t[2][idChunk]);
        reduce(ba, idChunk);
        
        __TIME_MARK(t[3][idChunk]);

        g.copy(chunkResults[idChunk], ba.getValue(chunkResultRef));

#ifdef __FULL_STATS__    
        stats[1][idChunk] = ba.stats;
        ba.clearStats();
#endif
        __TIME_MARK(t[4][idChunk]);
    }
    
    __TIME_MARK(t1);
    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r, r);
        g.add(r, r, chunkResults[j]);
    }

    __TIME_MARK(t2);
#ifdef __FULL_STATS__    
    double factor = (t2 - t0) > 2000000 ? 1000000.0: 1000.0;

    if (factor == 1000.0) {
        printf("NOTE: time in milliseconds\n");
    }
    if (factor == 1000000.0) {
        printf("NOTE: time in seconds\n");        
    }
    printf("n:%'d threads-max:%d chunks:%d accsPerChunk:%'ld\n", n, omp_get_max_threads(), nChunks, accsPerChunk);
    printf("iCh|total   |process |calculat|reduce  |copy    |start   |wait    \n");
    printf("---|MaxValues   |multC|singC|shortC|minB|maxBlocks|minA|maxAdds     |totalAdds   |resizes\n");
    for (uint32_t idChunk = 0; idChunk < nChunks; ++idChunk) {
        printf("#%2d|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f|%8.03f\n", idChunk, 
            ((double)(t[4][idChunk]-t[0][idChunk]))/factor,
            ((double)(t[1][idChunk]-t[0][idChunk]))/factor,
            ((double)(t[2][idChunk]-t[1][idChunk]))/factor,
            ((double)(t[3][idChunk]-t[2][idChunk]))/factor,
            ((double)(t[4][idChunk]-t[3][idChunk]))/factor,
            ((double)(t[0][idChunk]-t0))/factor,
            ((double)(t1-t[4][idChunk]))/factor);
        for (int j=0; j<2; ++j) {
            printf("---|%'12ld|%5ld|%5ld|%6ld|%4ld|%9ld|%4ld|%'12ld|%'12ld|%7d\n",
                stats[j][idChunk].maxValues,
                stats[j][idChunk].multiAddOperations,
                stats[j][idChunk].singleMultiAdds,
                stats[j][idChunk].shortMultiAdds,
                stats[j][idChunk].minBlock,
                stats[j][idChunk].maxBlock,
                stats[j][idChunk].minMultiAddValues,
                stats[j][idChunk].maxMultiAddValues,
                stats[j][idChunk].totalMultiAddValues,
                stats[j][idChunk].resizes         
                );
        }
    }
    printf("final: %9.04f\n", ((double)(t2-t1))/factor);
#endif
}
