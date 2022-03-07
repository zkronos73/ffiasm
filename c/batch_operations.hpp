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

struct Element
{
    OperationType operation;
    Index oper1;
    Index oper2;
    int referenceCount;
    bool evaluate;
    int evaluated;
    int evaluatedListIndex;
    int64_t evaluatedFromListIndex;
    int64_t evaluatedToListIndex;
    int evaluationId;
};

class Stack
{
    public:
        Stack ( void );
        ~Stack ( void );
        Reference alloc ( int count = 1 );
        void copy ( Reference destRef, Reference srcRef );
        void add ( Reference destRef, Reference opRef1, Reference opRef2 );
        void isZero ( Reference ref );
        void setZero ( Reference ref );
        void gMulByScalar ( Reference destRef, Reference baseRef, uint8_t* scalar, unsigned int scalarSize );
        void dbl ( Reference destRef, Reference srcRef );
        int64_t resolve ( Reference reference, int maxDeepLevel = 1000000 );
        std::string getPathToResolve ( Reference reference , int maxDeepLevel = 1000000 );
        void mark ( Reference refence );
        void clearMarks ( void );

    protected:

        std::vector<Element> operations;
        std::vector<Index> references;
        int64_t evaluationId;
        
        Reference push ( Reference destRef, OperationType operation, Index opRef1 = 0, Index opRef2 = 0 );
        Index getReferenceIndex ( Reference reference );
        void updateReferenceIndex ( Reference reference, Index operationReference );
        void updateReferencesCount ( Index index, int delta = 1 );
        int getReferencesCount ( Index index );
        int64_t calculate ( Index index, int maxDeepLevel );
        void set ( Reference ref, int64_t value );
        void dump ( void );
        std::string getOperationLabel ( OperationType tp );
        std::string getPathToResolveIndex ( Index index, int maxDeepLevel );
        void getOperationsLists ( OperationsLists &list, Index index, int maxDeepLevel );        
        void dumpOperationsLists ( OperationsLists &list );
        void clearEvaluated ( void );
        void normalizeOperationsLists ( OperationsLists &lists );
};
} 
#endif