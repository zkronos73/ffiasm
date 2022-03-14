#ifndef PAR_MULTIEXP2_BA
#define PAR_MULTIEXP2_BA

#define PME2_PACK_FACTOR 2
#define PME2_MAX_CHUNK_SIZE_BITS 16
#define PME2_MIN_CHUNK_SIZE_BITS 2

#include "batch_accumulators.hpp"

template <typename Curve>
class ParallelMultiexpBa
{
    const typename Curve::PointAffine *bases;
    typedef BatchAccumulators<Curve> BatchAcc;

    const uint8_t* scalars;
    uint32_t scalarSize;
    uint32_t n;
    uint32_t bitsPerChunk;
    uint64_t accsPerChunk;
    uint32_t nChunks;
    Curve &g;

    int64_t resultRef;
    int64_t chunkResultRef;

    uint32_t getChunk ( uint32_t scalarIdx, uint32_t chunkIdx );
    void processChunks ( BatchAcc &ba, uint32_t idChunk );
    void reduce ( BatchAcc &ba, uint32_t idChunk);

public:
    ParallelMultiexpBa(Curve &_g): g(_g) {}
    void multiexp(typename Curve::Point &r, const typename Curve::PointAffine *_bases, const uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads=0);
};

#include "multiexp_ba.cpp"

#endif // PAR_MULTIEXP2
