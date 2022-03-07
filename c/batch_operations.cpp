#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <string>

#include "batch_operations.hpp"

namespace BatchOperations {

Stack::Stack ( void )
{
    Element op = { op_set, 0, 0, 0, false };

    // TODO: multithread protection
    int reference = operations.size();
    operations.push_back(op);

    references.push_back(0);
    references.push_back(0);
    references.push_back(0);
}

Stack::~Stack ( void )
{

}

Reference Stack::alloc ( int count )
{
    auto reference = references.size();
    references.resize(references.size() + count, 0);
    updateReferencesCount(0, count);
    return reference;
}

void Stack::copy ( Reference destRef, Reference srcRef )
{
    auto srcIndex = getReferenceIndex(srcRef);
    updateReferenceIndex(destRef, srcIndex);
}

void Stack::updateReferencesCount ( Index index, int delta )
{
    printf("@DEBUG operations[%ld/%ld]\n", index, operations.size());

    operations[index].referenceCount += delta;
}

int Stack::getReferencesCount ( Index index )
{
    return operations[index].referenceCount;
}

void Stack::add ( Reference destRef, Reference ref1, Reference ref2 )
{
    auto index1 = getReferenceIndex(ref1); 
    auto index2 = getReferenceIndex(ref2); 
    push(destRef, op_add, index1, index2);
}

void Stack::isZero ( Reference ref )
{

}

void Stack::setZero ( Reference ref )
{

}

void Stack::gMulByScalar ( Reference destRef, Reference baseRef, uint8_t* scalar, unsigned int scalarSize )
{
    // TODO: optimization 1, 0
}

void Stack::dbl ( Reference destRef, Reference srcRef )
{
    // TODO: optimization 1, 0
    auto opRef1 = getReferenceIndex(srcRef); 
    push(destRef, op_dbl, opRef1);
}

Reference Stack::push ( Reference destRef, OperationType operation, Index opRef1, Index opRef2 )
{
    Element op = { operation, opRef1, opRef2, 0, false};

    // TODO: multithread protection
    int reference = operations.size();
    operations.push_back(op);
    updateReferenceIndex(destRef, reference);
    return reference;
}

Index Stack::getReferenceIndex ( Reference reference )
{
    return references[reference];
}

void Stack::updateReferenceIndex ( Reference reference, Index index )
{
    auto dstIndex = getReferenceIndex(reference);
    updateReferencesCount(index, 1);
    printf("@DEBUG references[%ld/%ld] = %ld\n", reference, references.size(), index);
    references[reference] = index;
    updateReferencesCount(dstIndex, -1);
}

int64_t Stack::resolve ( Reference reference, int maxDeepLevel )
{
    auto index = getReferenceIndex(reference);
    return calculate(index, maxDeepLevel );
}

std::string Stack::getPathToResolve ( Reference reference , int maxDeepLevel )
{
    auto index = getReferenceIndex(reference);
    return getPathToResolveIndex(index, maxDeepLevel);
}

int64_t Stack::calculate (Index index, int maxDeepLevel )
{
    Element e = operations[index];
    if (maxDeepLevel > 0) {
        switch (e.operation) {
            case op_set:
                return e.oper1;
            case op_add:
                return calculate(e.oper1, maxDeepLevel-1) + calculate(e.oper2, maxDeepLevel-1);
            case op_dbl:
                return 2 * calculate(e.oper1, maxDeepLevel-1);
        }
    }
    return -1;
}

void Stack::set ( Reference ref, int64_t value )
{
    push(ref, op_set, value);
}

void Stack::dump ( void )
{
    printf("**** REFERENCES *****\n");
    for(auto index = 0; index < references.size(); ++ index) {
        printf("[%8d] %ld\n", index, references[index]);
    }
    printf("\n**** OPERATIONS *****\n");
    for(auto index = 0; index < operations.size(); ++ index) {
        Element e = operations[index];
        printf("%8d|OP:%8s|O1:%8ld|O2:%8ld|RC:%8d|E:%c|#E:%8d|EL:%4d|EF:%8ld|ET:%8ld|EI:%8d\n", 
                index, getOperationLabel(e.operation).c_str(), e.oper1, e.oper2, e.referenceCount, e.evaluate ? 'T':'F', e.evaluated,
                e.evaluatedListIndex, e.evaluatedFromListIndex, e.evaluatedToListIndex, e.evaluationId);
    }
}

std::string Stack::getOperationLabel ( OperationType tp )
{
    switch (tp) {
        case op_dbl: return "op_dbl";
        case op_add: return "op_add";
        case op_set: return "op_set";
    }
    std::stringstream ss;
    ss << "unknown type (" << tp << ")";
    return ss.str();
}

std::string Stack::getPathToResolveIndex ( Index index, int maxDeepLevel )
{
    Element e = operations[index];
    if (maxDeepLevel > 0) {
        std::stringstream ss;
        switch (e.operation) {
            case op_set:
                ss << e.oper1;            
                return ss.str();
            case op_add:
                return getPathToResolveIndex(e.oper1, maxDeepLevel-1) + " + " + getPathToResolveIndex(e.oper2, maxDeepLevel-1);
            case op_dbl:
                return "( 2 * " +getPathToResolveIndex(e.oper1, maxDeepLevel-1)+" )";
        }
    }
    return "[UPPS!!]";    
}

void Stack::normalizeOperationsLists ( OperationsLists &lists )
{
//    auto offsetEvaluationId = lists.size();
//    lists.resize(offsetEvaluationId + evaluationId);
}

void Stack::getOperationsLists ( OperationsLists &lists, Index index, int maxDeepLevel )
{
    auto listIndex = lists.lists.size() - 1;

    if (!maxDeepLevel) {
        lists.lists[listIndex].push_back({999999, index, -1});
        return;
    }

    Element &e = operations[index];
    bool preEvaluated = e.evaluated > 0;
    
    switch (e.operation) {
        case op_set:
            if (listIndex < 0) {
                lists.add();
                ++listIndex;
            }
            lists.lists[listIndex].push_back({e.oper1, index, 0});
            return;

        case op_add:                
        {
            ++e.evaluated;
            if (preEvaluated && (e.evaluatedToListIndex - e.evaluatedFromListIndex) > 1) {
                if (!e.evaluationId) {
                    lists.evaluations.push_back(OperationsList(lists.lists[e.evaluatedListIndex].begin() + e.evaluatedFromListIndex,
                                                               lists.lists[e.evaluatedListIndex].begin() + e.evaluatedToListIndex));
                    e.evaluationId = lists.evaluations.size();                    
                }
                // TODO: provisional
                lists.lists[listIndex].push_back({0, index, e.evaluationId});
                return;
            }
            auto evaluatedFromListIndex = lists.lists[listIndex].size();
            getOperationsLists(lists, e.oper1, maxDeepLevel - 1);
            getOperationsLists(lists, e.oper2, maxDeepLevel - 1);
//            if ((lists.lists[listIndex].size() - evaluatedFromListIndex) > 3) {
                e.evaluatedListIndex = listIndex;
                e.evaluatedFromListIndex = evaluatedFromListIndex;
                e.evaluatedToListIndex = lists.lists[listIndex].size() - 1;
//            }
            return;
        }
        case op_dbl:
            getOperationsLists(lists, e.oper1, maxDeepLevel - 1);
            getOperationsLists(lists, e.oper1, maxDeepLevel - 1);
            return;
    }
}

void Stack::dumpOperationsLists ( OperationsLists &lists )
{
    std::stringstream ss;
    auto listCount = lists.lists.size();

    for (auto ilist = 0; ilist < lists.lists.size(); ++ilist) {
        ss << "list[" << ilist << "] = {";
        auto &operations = lists.lists[ilist];
        auto operCount = operations.size();
        for (auto ioper = 0; ioper < operCount; ++ioper) {
            ss << (ioper ? ",":"");
            const Operation &op = operations[ioper];
            if (op.evaluationId > 0) {
                ss << "@" << op.evaluationId;
            }
            else {
                ss << op.value;
            }
        }
        ss << "}\n";
    }
    for (auto ieval = 0; ieval < lists.evaluations.size(); ++ieval) {
        ss << "evaluations[" << ieval << "] = {";
        auto &operations = lists.evaluations[ieval];
        auto operCount = operations.size();
        for (auto ioper = 0; ioper < operCount; ++ioper) {
            ss << (ioper ? ",":"");
            const Operation &op = operations[ioper];
            if (op.evaluationId > 0) {
                ss << "@" << op.evaluationId;
            }
            else {
                ss << op.value;
            }
        }
        ss << "}\n";
    }
    printf("%s\n", ss.str().c_str());
}

void Stack::mark ( Reference reference )
{
    auto index = getReferenceIndex(reference);
    operations[index].evaluate = true;
}

void Stack::clearMarks ( void )
{
    auto size = operations.size();
    for (auto index = 0; index < size; ++index) {
        // TODO: performance operations[index].evaluate for all
        if (operations[index].evaluate) operations[index].evaluate = false;
    }
}

void Stack::clearEvaluated ( void )
{
    auto size = operations.size();
    for (auto index = 0; index < size; ++index) {
        Element &e = operations[index];
        e.evaluated = 0;
        e.evaluatedFromListIndex = 0;
        e.evaluatedListIndex = -1;
        e.evaluatedToListIndex = 0;
        e.evaluationId = 0;
    }
    evaluationId = 0;
}

}