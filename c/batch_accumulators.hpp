#ifndef __FFIASM__BATCH_ACCUMULATORS__H__
#define __FFIASM__BATCH_ACCUMULATORS__H__


#define BATCH_ACCUMULATORS_STATS 

typedef struct {
    int64_t maxValues;
    int64_t singleMultiAdds;
    int64_t shortMultiAdds;
    int64_t minBlock;
    int64_t maxBlock;
    int64_t multiAddOperations;
    int64_t maxMultiAddValues;
    int64_t minMultiAddValues;
    int64_t totalMultiAddValues;
    int32_t resizes;
} BatchAccumulatorsStats;

template <typename Curve>
class BatchAccumulators
{
    public:

        BatchAccumulatorsStats stats;

        Curve &g;
        BatchAccumulators(Curve &_g);
        void setup(int64_t _initialValues, int64_t _deltaValues);

        // TODO: value as const, but need modify Curve::copy
        void addPoint(int64_t accumulatorId, const typename Curve::PointAffine &value);
        void add(int64_t accumulatorId, int64_t valueAccumulatorId);
        void dbl(int64_t accumulatorId );
        bool calculateOnlyOneLoop ( void );
        void calculate( void ) { while (!calculateOnlyOneLoop()); }
        void clear ( void );
        bool ready ( void );
        typename Curve::PointAffine isValueReady ( int64_t accumulatorId );
        const typename Curve::PointAffine &getValue ( int64_t accumulatorId );
        bool isZero ( int64_t accumulatorId );
        int64_t defineAccumulators ( int64_t count );
        void dumpStats ( void );
        void clearStats ( void );

    protected:

        typedef struct {
            typename Curve::PointAffine value;
            typename Curve::PointAffine singleValue;
            int32_t lastLoop;
            bool ready;
        } Accumulator;

        Accumulator *accumulators;
        int64_t accumulatorsCount;
        int64_t valuesCount;
        int64_t valuesSize;
        int64_t initialValues;
        int64_t deltaValues;
        int32_t currentLoop;

        typename Curve::PointAffine *leftValues;
        typename Curve::PointAffine *rightValues;
        typename Curve::PointAffine *resultValues;
        int64_t *accumulatorIds;

        typename Curve::PointAffine zero;

        void freeValues ( void );
        void resize ( void );
        void internalAdd ( int64_t accumulatorId, const typename Curve::PointAffine &value, bool internal );
        void multiAdd ( void );
        int64_t incValuesCount ( void );
};

#include "batch_accumulators.cpp"

#endif

