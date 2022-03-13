
#include <omp.h>
#include <memory.h>
#include "misc.hpp"

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
            g.add(accs[chunkValue].p, accs[chunkValue].p, bases[i]);
        }
    }
}

template <typename Curve>
void ParallelMultiexpMix<Curve>::reduce(typename Curve::Point &res, uint32_t nBits) 
{
    uint32_t ndiv2;

    g.copy(res, g.zero());
    while (true) {
        ndiv2 = 1 << (nBits-1);
        if (nBits==1) {
            g.add(res, res, accs[1].p);
//            printf("==> MIX-CHUNK-NDIV2 A = 2**%d * %s\n", nBits-1, g.toString(accs[ndiv2].p).c_str());
            //printf("==> MIX-CHUNK-NDIV2 B = 2**%d * %s\n", nBits-1, g.toString(accs[ndiv2].p).c_str());
            g.copy(accs[1].p, g.zero());
            return;
        }

        for (uint32_t i = 1; i<ndiv2; i++) {
            if (!g.isZero(accs[ndiv2 + i].p)) {
                g.add(accs[i].p, accs[i].p, accs[ndiv2 + i].p);
                g.add(accs[ndiv2].p, accs[ndiv2].p, accs[ndiv2 + i].p);
                g.copy(accs[ndiv2 + i].p, g.zero());
            }
        }
        for (u_int32_t i=0; i<nBits-1; i++) g.dbl(accs[ndiv2].p, accs[ndiv2].p);
//         printf("==> MIX-CHUNK-NDIV2 B = 2**%d * %s\n", nBits-1, g.toString(accs[ndiv2].p).c_str());
//        printf("==> MIX-CHUNK-ADD-RESULT[%d] = %s + %s\n", nBits, g.toString(res).c_str(), g.toString(accs[ndiv2].p).c_str());    

        g.add(res, res, accs[ndiv2].p);
        g.copy(accs[ndiv2].p, g.zero());
        --nBits;
    }    
}

template <typename Curve>
void ParallelMultiexpMix<Curve>::multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) {
//    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
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

    typename Curve::Point *chunkResults = new typename Curve::Point[nChunks];
    accs = new PaddedPoint[accsPerChunk];
    initAccs();

    printf("nChunks:%d\n", nChunks);
    for (uint32_t i=0; i<nChunks; i++) {
        processChunk(i);
/*        for (uint32_t j=1; j<accsPerChunk; j++) {
            printf("==> MIX-CHUNK[%d-%d] = %s\n", i, j, g.toString(accs[j].p).c_str());
        }*/
/*        for (int index = 0; index < n; ++index) {
            printf("accs[%d].p = %s\n", index, g.toString(accs[index].p).c_str());
        }
        return;*/
        reduce(chunkResults[i], bitsPerChunk);
//        memset(accs, 0, sizeof(accs[0]) * accsPerChunk);
    }
    // g.printCounters();
    // return;
/*
    for (uint32_t i=0; i<nChunks; i++) {
        typename Curve::Point p = chunkResults[i];
        printf("==> MIX-CHUNK-RESULT[%d] = %s\n", i, g.toString(p).c_str());
        // break;
    }    
*/
    delete[] accs;

    g.copy(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r,r);
        g.add(r, r, chunkResults[j]);
    }

    delete[] chunkResults; 
}
