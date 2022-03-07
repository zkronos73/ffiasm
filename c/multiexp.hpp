#ifndef PAR_MULTIEXP2
#define PAR_MULTIEXP2

#define PME2_PACK_FACTOR 2
#define PME2_MAX_CHUNK_SIZE_BITS 16
#define PME2_MIN_CHUNK_SIZE_BITS 2

template <typename Curve>
class ParallelMultiexp {

    struct PaddedPoint {
        typename Curve::Point p;
//        uint8_t padding[32];
    };

    typedef typename Curve::Point Point;
    typedef typename Curve::PointAffine PointAffine;
    typename Curve::PointAffine *bases;
    uint8_t* scalars;
    uint32_t scalarSize;
    uint32_t n;
    uint32_t nThreads;
    uint32_t bitsPerChunk;
    uint64_t accsPerChunk;
    uint32_t nChunks;
    Curve &g;
    PaddedPoint *accs;
    Point *chunkResults;
    bool cache;

    void initAccs();

    uint32_t getChunk(uint32_t scalarIdx, uint32_t chunkIdx);
    void processChunk(uint32_t idxChunk);
    void packThreads();
    void reduce(typename Curve::Point &res, uint32_t nBits);

public:
    ParallelMultiexp(Curve &_g): g(_g) {}
    void multiexp(typename Curve::Point &r, typename Curve::PointAffine *_bases, uint8_t* _scalars, uint32_t _scalarSize, uint32_t _n, uint32_t _nThreads=0);

    inline bool checkPos(const char *name, int line, void *data) { 
        bool result = !cache || 
               (data >= accs && data <= (accs + (nThreads*accsPerChunk))) ||
               (data >= bases && data <= (bases + n)) || 
               (data >= chunkResults && data <= (chunkResults + nChunks));
        if (!result) printf("@:%d %s %p\n", line, name, data);
        return result;
    }
    inline void gCopy(int line, Point &r, Point &a) { checkPos(__func__, line, &r); g.copy(r, a); };
    inline void gMulByScalar(int line, Point &r, PointAffine &base, uint8_t* scalar, unsigned int scalarSize) { checkPos(__func__, line, &r);g.mulByScalar(r, base, scalar, scalarSize); };
    inline void gAdd(int line, Point &p3, Point &p1, Point &p2) { checkPos(__func__, line,  &p3);g.add(p3, p1, p2); };
    inline void gAdd(int line, Point &p3, Point &p1, PointAffine &p2) { checkPos(__func__, line, &p3);g.add(p3, p1, p2); };
    inline bool gIsZero(int line, Point &p1) { return g.isZero(p1); };
    inline bool gIsZero(int line, PointAffine &p1) { return g.isZero(p1); };
    inline void gDbl(int line, Point &r, Point &a) { checkPos(__func__, line, &r); g.dbl(r, a); };
    inline void gDbl(int line, Point &r, PointAffine &a) { checkPos(__func__, line, &r);g.dbl(r, a); };


};

#include "multiexp.cpp"

#endif // PAR_MULTIEXP2