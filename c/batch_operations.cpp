#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <stack>

#include "batch_operations.hpp"

namespace BatchOperations {

template <typename Curve>
Stack<Curve>::Stack (Curve &_g)
    :g(_g)
{
    Element op = { op_set, 0, 0, 0, false };

    // TODO: multithread protection
    int reference = operations.size();
    operations.push_back(op);

    g.copy(zeroValue, g.zero());
    references.push_back({0, zeroValue});
    defines.push_back({"ZERO", 0, 1});
    addCounter = 0;
    dblCounter = 0;
    addCallCounter = 0;
    dblCallCounter = 0;
}

template <typename Curve>
Stack<Curve>::~Stack ( void )
{
}

template <typename Curve>
Reference Stack<Curve>::define ( const std::string &label, int size )
{
    Reference reference = references.size();
    
    references.resize(references.size() + size, {0, zeroValue});
    defines.push_back({label, reference, size});
    updateReferencesCount(0, size);
    return reference;
}

template <typename Curve>
std::string Stack<Curve>::getReferenceLabel ( Reference reference )
{
    for (int index = 0; index < defines.size(); ++index) {
        const Define &_define = defines[index];
        if (_define.reference <= reference && (_define.reference + _define.size) > reference) {
            std::stringstream ss;
            ss << _define.label;
            if (index || _define.size > 1) {
                ss << "[" << (reference - _define.reference) << "]";
            }
            return ss.str();
        }
    }
    std::stringstream ss;
    ss << "UNDEFINED @" << reference ;
    return ss.str();
}

template <typename Curve>
void Stack<Curve>::copy ( Reference destRef, Reference srcRef )
{
    auto srcIndex = getReferenceIndex(srcRef);
    // printf("COPY(%ld, %ld) SRC:%ld\n", destRef, srcRef, srcIndex);
    updateReferenceIndex(destRef, srcIndex);
}

template <typename Curve>
void Stack<Curve>::updateReferencesCount ( Index index, int delta )
{
    // printf("@DEBUG operations[%ld/%ld]\n", index, operations.size());

    operations[index].referenceCount += delta;
}

template <typename Curve>
int Stack<Curve>::getReferencesCount ( Index index )
{
    return operations[index].referenceCount;
}

template <typename Curve>
void Stack<Curve>::add ( Reference destRef, Reference ref1, Reference ref2 )
{
    auto index1 = getReferenceIndex(ref1);
    auto index2 = getReferenceIndex(ref2); 
    if (!index1) {
        updateReferenceIndex(destRef, index2);
        return;
    }
    if (!index2) {
        updateReferenceIndex(destRef, index1);
        return;
    }
    ++addCounter;
    printf("@ADD(%s:%ld, %s:%ld) => %s\n", getReferenceLabel(ref1).c_str(), index1, getReferenceLabel(ref2).c_str(), index2, getReferenceLabel(destRef).c_str());
    push(destRef, op_add, index1, index2);
}

template <typename Curve>
void Stack<Curve>::isZero ( Reference ref )
{

}

template <typename Curve>
void Stack<Curve>::setZero ( Reference ref )
{

}

template <typename Curve>
void Stack<Curve>::gMulByScalar ( Reference destRef, Reference baseRef, uint8_t* scalar, unsigned int scalarSize )
{
    // TODO: optimization 1, 0
}

template <typename Curve>
void Stack<Curve>::dbl ( Reference destRef, Reference srcRef )
{
    // TODO: optimization 1, 0
    ++dblCounter;
    auto opRef1 = getReferenceIndex(srcRef); 
    printf("@DBL(%s:%ld) => %s\n", getReferenceLabel(srcRef).c_str(), opRef1, getReferenceLabel(destRef).c_str());
    if (!opRef1) {
        updateReferenceIndex(destRef, opRef1);
        return;
    }
    push(destRef, op_dbl, opRef1);
}

template <typename Curve>
Reference Stack<Curve>::push ( Reference destRef, OperationType operation, Index opRef1, Index opRef2 )
{
    Element op = { operation, opRef1, opRef2, 0, false};

    // TODO: multithread protection
    int reference = operations.size();
    operations.push_back(op);
    updateReferenceIndex(destRef, reference);
    return reference;
}

template <typename Curve>
Index Stack<Curve>::getReferenceIndex ( Reference reference )
{
    return references[reference].index;
}


template <typename Curve>
const typename Curve::PointAffine &Stack<Curve>::getReferenceValue ( Reference reference )
{
    return references[reference].value;
}

template <typename Curve>
void Stack<Curve>::updateReferenceIndex ( Reference reference, Index index )
{
    auto dstIndex = getReferenceIndex(reference);
    updateReferencesCount(index, 1);
    // printf("@DEBUG references[%ld/%ld] = %ld\n", reference, references.size(), index);
    references[reference].index = index;
    updateReferencesCount(dstIndex, -1);
}

template <typename Curve>
typename Curve::PointAffine Stack<Curve>::resolve ( Reference reference, int maxDeepLevel )
{
    auto index = getReferenceIndex(reference);
    typename Curve::PointAffine result = iterativeCalculate(index);
    return result;
}

template <typename Curve>
std::string Stack<Curve>::getPathToResolve ( Reference reference , int maxDeepLevel )
{
    auto index = getReferenceIndex(reference);
    return getPathToResolveIndex(index, maxDeepLevel);
}

template <typename Curve>
typename Curve::PointAffine Stack<Curve>::iterativeCalculate (Index index, int level )
{
    std::stack<Index> st;

    st.push(index);
    while (!st.empty()) {
        Element &e = operations[st.top()];
        if (e.evaluated) {
            st.pop();
            continue;
        }

        std::string indent((st.size() - 1) * 2, ' ');
        switch (e.operation) {
            case op_set:
                e.evaluated = true;
                e.value = getReferenceValue(e.oper1);
                printf("@C %sset(%ld)\n", indent.c_str(), e.oper1);
                st.pop();
                continue;
                
            case op_add:
                {   
                    if (!operations[e.oper1].evaluated) {
                        st.push(e.oper1);
                        continue;
                    }
                    if (!operations[e.oper2].evaluated) {
                        st.push(e.oper2);
                        continue;
                    }
                    ++addCallCounter;
                    g.add(e.value, operations[e.oper1].value, operations[e.oper2].value);
                    printf("@C %sadd(%ld, %ld) => %s\n", indent.c_str(), e.oper1, e.oper2, g.toString(e.value).c_str());
                    e.evaluated = true;
                    st.pop();
                    continue;
                }
            case op_dbl:
                {
                    if (!operations[e.oper1].evaluated) {
                        st.push(e.oper1);
                        continue;
                    }
                    ++dblCallCounter;
                    g.dbl(e.value, operations[e.oper1].value);
                    printf("@C %sdbl(%ld) => %s\n", indent.c_str(), e.oper1, g.toString(e.value).c_str());
                    e.evaluated = true;
                    st.pop();
                }
        }
    }
    return operations[index].evaluated ? operations[index].value : zeroValue;
}

template <typename Curve>
typename Curve::PointAffine Stack<Curve>::recursiveCalculate (Index index, int level )
{
    Element &e = operations[index];
    typename Curve::PointAffine ec;
    std::string indent(level*2, ' ');

    if (e.evaluated) {
        return e.value;
    }
    switch (e.operation) {
        case op_set:                            
            ec = getReferenceValue(e.oper1);
            // printf("@C %sget(%s:%ld) = %s\n", indent.c_str(), getReferenceLabel(e.oper1).c_str(), e.oper1, g.toString(ec).c_str());
            return ec;
        case op_add:
            {   
                ++addCallCounter;
                printf("@C %sadd(%ld, %ld)\n", indent.c_str(), e.oper1, e.oper2);
                typename Curve::PointAffine oper1 = calculate(e.oper1, level+1);
                typename Curve::PointAffine oper2 = calculate(e.oper2, level+1);
                g.add(e.value, oper1 ,oper2);
                e.evaluated = true;
                // printf("\n## %s\n +%s\n = %s\n", g.toString(oper1).c_str(), g.toString(oper2).c_str(), g.toString(ec).c_str());
                return e.value;
            }
        case op_dbl:
            {
                ++dblCallCounter;
                printf("@C %sdbl(%ld)\n", indent.c_str(), e.oper1);
                typename Curve::PointAffine oper1 = calculate(e.oper1, level+1);
                g.dbl(e.value, oper1);
                e.evaluated = true;
                // printf("\n## 2 * %s\n# = %s\n", g.toString(oper1).c_str(), g.toString(ec).c_str());
                return e.value;
            }
    }
    printf("ABORTING !!!!!!!\n");
    return zeroValue;
}

template <typename Curve>
void Stack<Curve>::set ( Reference ref, typename Curve::PointAffine value )
{
    references[ref].value = value;
    push(ref, op_set, ref);
}

template <typename Curve>
void Stack<Curve>::dump ( void )
{
    printf("**** REFERENCES *****\n");
    for(auto index = 0; index < references.size(); ++ index) {
        printf("[%8d] %ld\n", index, references[index].index);
        // , g.toString(references[index].value).c_str();
    }
    printf("\n**** OPERATIONS *****\n");
    for(auto index = 0; index < operations.size(); ++ index) {
        Element e = operations[index];
        printf("%8d|OP:%8s|O1:%8ld|O2:%8ld|RC:%8d|E:%c|#E:%8d|EL:%4d|EF:%8ld|ET:%8ld|EI:%8d\n", 
                index, getOperationLabel(e.operation).c_str(), e.oper1, e.oper2, e.referenceCount, e.evaluate ? 'T':'F', e.evaluated,
                e.evaluatedListIndex, e.evaluatedFromListIndex, e.evaluatedToListIndex, e.evaluationId);
    }
}

template <typename Curve>
std::string Stack<Curve>::getOperationLabel ( OperationType tp )
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

template <typename Curve>
std::string Stack<Curve>::getPathToResolveIndex ( Index index, int maxDeepLevel )
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

template <typename Curve>
void Stack<Curve>::normalizeOperationsLists ( OperationsLists &lists )
{
//    auto offsetEvaluationId = lists.size();
//    lists.resize(offsetEvaluationId + evaluationId);
}

template <typename Curve>
void Stack<Curve>::getOperationsLists ( OperationsLists &lists, Index index, int maxDeepLevel )
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

template <typename Curve>
void Stack<Curve>::dumpOperationsLists ( OperationsLists &lists )
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

template <typename Curve>
void Stack<Curve>::mark ( Reference reference )
{
    auto index = getReferenceIndex(reference);
    operations[index].evaluate = true;
}

template <typename Curve>
void Stack<Curve>::clearMarks ( void )
{
    auto size = operations.size();
    for (auto index = 0; index < size; ++index) {
        // TODO: performance operations[index].evaluate for all
        if (operations[index].evaluate) operations[index].evaluate = false;
    }
}

template <typename Curve>
void Stack<Curve>::clearEvaluated ( void )
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