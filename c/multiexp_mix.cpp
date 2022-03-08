#include <omp.h>
#include <memory.h>
#include "misc.hpp"

//    auto _O1 = O1; 
//    auto _O2 = O2; 
//    printf("\n## %d: %s\n## %s\n +%s\n = %s\n", __LINE__, #R "=" #O1 "+" #O2,g.toString(_O1).c_str(), g.toString(_O2).c_str(), g.toString(R).c_str()); 

#define __G_ADD__MIX__(R,O1,O2) \
{ \
    g.add(R, O1, O2); \
} 

template <typename Curve>
void ParallelMultiexpMix<Curve>::initAccs() {
    for (int i=0; i<accsPerChunk; i++) {
        g.copy(accs[i].p, g.zero());
    }
}

template <typename Curve>
uint32_t ParallelMultiexpMix<Curve>::getChunk(uint32_t scalarIdx, uint32_t chunkIdx) {
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
void ParallelMultiexpMix<Curve>::processChunk(uint32_t idChunk) {
    for(uint32_t i=0; i<n; i++) {
        if (g.isZero(bases[i])) continue;
        uint32_t chunkValue = getChunk(i, idChunk);
        if (chunkValue) {
            __G_ADD__MIX__(accs[chunkValue].p, accs[chunkValue].p, bases[i]);
        }
    }
}

template <typename Curve>
void ParallelMultiexpMix<Curve>::reduce(typename Curve::Point &res, uint32_t nBits) 
{
    uint32_t ndiv2;

    g.copy(res, g.zero());
    while (true) {
        if (nBits==1) {
            g.add(res, res, accs[1].p);
            g.copy(accs[1].p, g.zero());
            return;
        }
        ndiv2 = 1 << (nBits-1);

        for (uint32_t i = 1; i<ndiv2; i++) {
            if (!g.isZero(accs[ndiv2 + i].p)) {
                __G_ADD__MIX__(accs[i].p, accs[i].p, accs[ndiv2 + i].p);
                __G_ADD__MIX__(accs[ndiv2].p, accs[ndiv2].p, accs[ndiv2 + i].p);
                g.copy(accs[ndiv2 + i].p, g.zero());
            }
        }
        for (u_int32_t i=0; i<nBits-1; i++) g.dbl(accs[ndiv2].p, accs[ndiv2].p);
        __G_ADD__MIX__(res, res, accs[ndiv2].p);
        g.copy(accs[ndiv2].p, g.zero());
        --nBits;
    }    
}

template <typename Curve>
void ParallelMultiexpMix<Curve>::multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) {
    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

    ThreadLimit threadLimit (nThreads);

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

    typename Curve::Point *chunkResults = new typename Curve::Point[nChunks];
    accs = new PaddedPoint[accsPerChunk];
    initAccs();

    printf("nChunks:%d\n", nChunks);
    for (uint32_t i=0; i<nChunks; i++) {
        processChunk(i);
/*        for (int index = 0; index < n; ++index) {
            printf("accs[%d].p = %s\n", index, g.toString(accs[index].p).c_str());
        }
        return;*/
        reduce(chunkResults[i], bitsPerChunk);
    }

    delete[] accs;

    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r,r);
        __G_ADD__MIX__(r, r, chunkResults[j]);
    }

    delete[] chunkResults; 
}
