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
    cache = false;
    #pragma omp parallel for
    for (int i=0; i<nThreads*accsPerChunk; i++) {
        G_COPY(accs[i].p, g.zero());
    }
    cache = true;
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
    #pragma omp parallel for
    for(uint32_t i=0; i<n; i++) {
        if (G_IS_ZERO(bases[i])) continue;
        int idThread = omp_get_thread_num();
        uint32_t chunkValue = getChunk(i, idChunk);
        if (chunkValue) {
            G_ADD(accs[idThread*accsPerChunk+chunkValue].p, accs[idThread*accsPerChunk+chunkValue].p, bases[i]);
        }
    }
}

template <typename Curve>
void ParallelMultiexp<Curve>::packThreads() {
    #pragma omp parallel for
    for(uint32_t i=0; i<accsPerChunk; i++) {
        for(uint32_t j=1; j<nThreads; j++) {
            if (!G_IS_ZERO(accs[j*accsPerChunk + i].p)) {
                G_ADD(accs[i].p, accs[i].p, accs[j*accsPerChunk + i].p);
                G_COPY(accs[j*accsPerChunk + i].p, g.zero());
            }
        }
    }
}

#ifndef __LEGACY__

template <typename Curve>
void ParallelMultiexp<Curve>::reduce(typename Curve::Point &res, uint32_t nBits) 
{
    G_COPY(res, g.zero());
    uint32_t ndiv2;
    
    while (true) {
        if (nBits == 1) {
            G_ADD(res, res, accs[1].p);
            G_COPY(accs[1].p, g.zero());
            return;
        }
        ndiv2 = 1 << (nBits - 1);

        for (uint32_t i = 1; i<ndiv2; i++) {
            if (!G_IS_ZERO(accs[ndiv2 + i].p)) {
                G_ADD(accs[i].p, accs[i].p, accs[ndiv2 + i].p);
                G_ADD(accs[ndiv2].p, accs[ndiv2].p, accs[ndiv2 + i].p);
                G_COPY(accs[ndiv2 + i].p, g.zero());
            }
        }
        for (u_int32_t i=0; i<(nBits - 1); i++) {
            G_DBL(accs[ndiv2].p, accs[ndiv2].p);
        }
        G_ADD(res, res, accs[ndiv2].p);
        G_COPY(accs[ndiv2].p, g.zero());
        --nBits;
    }
}

#else 

template <typename Curve>
void ParallelMultiexp<Curve>::reduce(typename Curve::Point &res, uint32_t nBits) 
{
    printf("@V1: %s\n", g.toString(accs[1].p).c_str());
    if (nBits==1) {
        G_COPY(res, accs[1].p);
        printf("@V[1]: %s\n", g.toString(accs[1].p).c_str());
        // G_COPY(accs[1].p, g.zero());
        return;
    }
    uint32_t ndiv2 = 1 << (nBits-1);


    // PaddedPoint *sall = new PaddedPoint[nThreads];
    // memset(sall, 0, sizeof(PaddedPoint)*nThreads);

//    #pragma omp parallel for
    for (uint32_t i = 1; i<ndiv2; i++) {
//        int idThread = omp_get_thread_num();
        if (!G_IS_ZERO(accs[ndiv2 + i].p)) {
            G_ADD(accs[i].p, accs[i].p, accs[ndiv2 + i].p);
            G_ADD(accs[ndiv2].p, accs[ndiv2].p, accs[ndiv2 + i].p);
            // g.add(sall[idThread].p, sall[idThread].p, accs[ndiv2 + i].p);
            G_COPY(accs[ndiv2 + i].p, g.zero());
        }
    }
    /*
    for (u_int32_t i=0; i<nThreads; i++) {
        G_ADD(accs[ndiv2].p, accs[ndiv2].p, sall[i].p);
    }
    */
    for (u_int32_t i=0; i<nBits-1; i++) G_DBL(accs[ndiv2].p, accs[ndiv2].p);

    typename Curve::Point p1;
    reduce(p1, nBits-1);

//    for (u_int32_t i=0; i<nBits-1; i++) G_DBL(accs[ndiv2].p, accs[ndiv2].p);
    G_ADD(res, p1, accs[ndiv2].p);
    printf("@V[%d]: %s\n", ndiv2,  g.toString(accs[ndiv2].p).c_str());
    // G_COPY(accs[ndiv2].p, g.zero());
//    delete[] sall;
}

#endif

template <typename Curve>
void ParallelMultiexp<Curve>::multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads) {
    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    bases = _bases;
    scalars = _scalars;
    scalarSize = _scalarSize;
    n = _n;

    ThreadLimit threadLimit (nThreads);
    printf("THREADS LIMIT: %d\n", nThreads);

    if (n==0) {
        G_COPY(r, g.zero());
        return;
    }
    if (n==1) {
        G_MUL_BY_SCALAR(r, bases[0], scalars, scalarSize);
        return;
    }
    bitsPerChunk = log2(n / PME2_PACK_FACTOR);
    if (bitsPerChunk > PME2_MAX_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MAX_CHUNK_SIZE_BITS;
    if (bitsPerChunk < PME2_MIN_CHUNK_SIZE_BITS) bitsPerChunk = PME2_MIN_CHUNK_SIZE_BITS;
    nChunks = ((scalarSize*8 - 1 ) / bitsPerChunk)+1;
    accsPerChunk = 1 << bitsPerChunk;  // In the chunks last bit is always zero.

    typename Curve::Point *chunkResults = new typename Curve::Point[nChunks];
    chunkResults = new typename Curve::Point[nChunks];
    accs = new PaddedPoint[nThreads*accsPerChunk];
    // std::cout << "InitTrees " << "\n"; 
    initAccs();

    for (uint32_t i=0; i<nChunks; i++) {
        // std::cout << "process chunks " << i << "\n"; 
        processChunk(i);
        // std::cout << "pack " << i << "\n"; 
        packThreads();
        // std::cout << "reduce " << i << "\n"; 
        reduce(chunkResults[i], bitsPerChunk);
        printf("@V ==========\n");
    }

    delete[] accs;
    cache = false;

    G_COPY(r, chunkResults[nChunks-1]);
    for  (int j=nChunks-2; j>=0; j--) {
        for (uint32_t k=0; k<bitsPerChunk; k++) g.dbl(r,r);
        G_ADD(r, r, chunkResults[j]);
    }

    delete[] chunkResults; 
}
