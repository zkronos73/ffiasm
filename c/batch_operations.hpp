#ifndef __FFIASM__OPERATION_STACK__H__
#define __FFIASM__OPERATION_STACK__H__

#include <stdint.h>
#include <vector>

namespace BatchOperations {

typedef int64_t Reference;
typedef int64_t Index;

typedef struct {
    int64_t value;
    Index operation;
    int evaluationId;
} Operation;

typedef std::vector<Operation> OperationsList;
// typedef std::vector<OperationsList> OperationsLists;


class OperationsLists
{
    public:
        std::vector<OperationsList> lists;
        std::vector<OperationsList> evaluations;
        void add ( void ) { lists.resize(lists.size() + 1);};
};

enum OperationType
{
    op_set,
    op_add,
    op_dbl
};


template <typename Curve>
class Stack
{
    public:
        Stack ( Curve &_g );
/*          :g(_g)
        {
            Element op = { op_set, 0, 0, 0, false };

            // TODO: multithread protection
            int reference = operations.size();
            operations.push_back(op);

            references.push_back(0);
            references.push_back(0);
            references.push_back(0);
        };*/
        typedef typename Curve::PointAffine Value;

        ~Stack ( void );
        Reference define ( const std::string &label, int size = 1 );
        void copy ( Reference destRef, Reference srcRef );
        void add ( Reference destRef, Reference opRef1, Reference opRef2, const std::string &source, int line );
        void isZero ( Reference ref );
        void setZero ( Reference ref );
        void gMulByScalar ( Reference destRef, Reference baseRef, uint8_t* scalar, unsigned int scalarSize );
        void dbl ( Reference destRef, Reference srcRef, const std::string &source, int line );
        typename Curve::PointAffine resolve ( Reference reference, int maxDeepLevel = 1000000 );
        std::string getPathToResolve ( Reference reference , int maxDeepLevel = 1000000 );
        std::string getReferenceLabel ( Reference reference );
        void getAdditionValues ( std::vector<Value> &values, Index index );
        void recursiveGetAdditionValues ( std::vector<Value> &values, Index index );


        void mark ( Reference refence );
        void clearMarks ( void );

        Value zeroValue;
        typedef struct {
            Index index;
            Value value;
        } ReferenceValue;

        typedef struct {
            std::string label;
            Reference reference;
            int64_t size;
        } Define;

        struct Element
        {
            OperationType operation;
            Index oper1;
            Index oper2;
            int referenceCount;
            int operationReferenceCount;
            bool evaluate;
            bool calculated;
            Value value;
            int evaluated;
            int evaluatedListIndex;
            int64_t evaluatedFromListIndex;
            int64_t evaluatedToListIndex;
            int evaluationId;
            std::string extraInfo;
        };
    // protected:
        Curve &g;
        std::vector<Element> operations;
        std::vector<ReferenceValue> references;
        std::vector<Define> defines;
        int64_t evaluationId;
        int64_t addCounter;
        int64_t dblCounter;
        int64_t addCallCounter;
        int64_t dblCallCounter;
        
        Reference push ( Reference destRef, OperationType operation, Index opRef1 = 0, Index opRef2 = 0, const std::string &extraInfo = "" );
        Index getReferenceIndex ( Reference reference );
        const typename Curve::PointAffine &getReferenceValue ( Reference reference );
        void updateReferenceIndex ( Reference reference, Index operationReference );
        void updateReferencesCount ( Index index, int delta = 1 );
        int getReferencesCount ( Index index );
        void updateOperationReferencesCount ( Index index, int delta = 1 );
        int getOperationReferencesCount ( Index index );
        typename Curve::PointAffine iterativeCalculate ( Index index );
        typename Curve::PointAffine recursiveCalculate ( Index index, int level = 0 );
        void set ( Reference ref, typename Curve::PointAffine value );
        void dump ( void );
        std::string getOperationInfo ( Index index );
        std::string getOperationLabel ( OperationType tp );
        std::string getPathToResolveIndex ( Index index, int maxDeepLevel );
        void getOperationsLists ( OperationsLists &list, Index index, int maxDeepLevel );        
        void dumpOperationsLists ( OperationsLists &list );
        void clearEvaluated ( void );
        void normalizeOperationsLists ( OperationsLists &lists );
        std::string toString (typename Curve::PointAffine &value) { typename Curve::Point p; g.copy(p, value); return g.toString(p); };
};

} 

#include "batch_operations.cpp"

#endif