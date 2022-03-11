#ifndef __FFIASM__BATCH_ACCUMULATORS__H__
#define __FFIASM__BATCH_ACCUMULATORS__H__


#define BATCH_ACCUMULATORS_STATS 

template <typename Curve>
class BatchAccumulators
{
    public:
        int64_t maxValues;
        int64_t multiAddOperations;
        int64_t maxMultiAddValues;
        int64_t minMultiAddValues;
        int64_t totalMultiAddValues;

        Curve &g;
        BatchAccumulators(Curve &_g);
        void setup(int64_t _initialValues, int64_t _deltaValues);

        // TODO: value as const, but need modify Curve::copy
        void add(int64_t accumulatorId, typename Curve::PointAffine &value);
        void add(int64_t accumulatorId, int64_t valueAccumulatorId);
        void dbl(int64_t accumulatorId );
        bool calculateOnlyOneLoop ( void );
        void calculate( void ) { while (!calculateOnlyOneLoop()); }
        void clear ( void );
        bool ready ( void );
        typename Curve::PointAffine isValueReady ( int64_t accumulatorId );
        typename Curve::PointAffine getValue ( int64_t accumulatorId );
        bool isZero ( int64_t accumulatorId );
        int64_t defineAccumulators ( int64_t count );
        void dumpStats ( void );
        void clearStats ( void );

    protected:

        typedef struct {
            int64_t index;
            int64_t count;
            bool ready;
            typename Curve::PointAffine value;
        } Accumulator;

        Accumulator *accumulators;
        int64_t accumulatorsCount;
        int64_t valuesCount;
        int64_t valuesSize;
        int64_t initialValues;
        int64_t deltaValues;

        typename Curve::PointAffine *leftValues;
        typename Curve::PointAffine *rightValues;
        typename Curve::PointAffine *resultValues;
        int64_t *accumulatorIds;

        void freeValues ( void );
        void prepareSingleValues ( void );
        void prepareAccumulatorsToNextRound ( void );
        void resize ( void );
        void internalAdd ( int64_t accumulatorId, typename Curve::PointAffine &value, bool internal );
        void multiAdd ( void );
};

#include "batch_accumulators.cpp"

#endif

