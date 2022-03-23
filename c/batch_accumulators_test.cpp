#include <gmp.h>
#include <iostream>

#include "gtest/gtest.h"
#include "batch_accumulators.hpp"

namespace {

class IntAsCurvePointAffine {
    public:
        int value;
        IntAsCurvePointAffine(void) { value = 0; };
        IntAsCurvePointAffine(int _value) { value = _value; };
        operator int() { return value; };
        operator std::string() { std::stringstream ss; ss << value; return ss.str(); };
        friend std::ostream& operator<<(std::ostream& os, const IntAsCurvePointAffine& value);
};

std::ostream& operator<<(std::ostream& os, const IntAsCurvePointAffine &value)
{
    os << value.value;
    return os;
}

class IntAsCurve {
    public:
        typedef IntAsCurvePointAffine PointAffine;
        void copy (PointAffine &dst, const PointAffine &src) { dst.value = src.value; };
        void add (PointAffine &dst, const PointAffine &left, const PointAffine &right) { dst.value = left.value + right.value; };
        const IntAsCurvePointAffine &zero (void) { return _zero; };
        const IntAsCurvePointAffine &zeroAffine (void) { return _zero; };
        bool eq (PointAffine op1, const PointAffine op2) { return op1.value == op2.value; };
        void multiAdd(PointAffine *res, const PointAffine *left, const PointAffine *right, int64_t count) { 
            printf("multiAdd(,%ld)\n", count);
            for(int64_t index = 0; index < count; ++index) {
                add(res[index], left[index], right[index]);
                std::cout << "R[" << index << "] = " << left[index] << " + " << right[index] << " = " << res[index] << "\n";
            }
        };
        std::string toString (PointAffine &value) { std::stringstream ss; ss << value; return ss.str(); };
        bool isZero(PointAffine &dst) { return (dst.value == 0); };
    protected:
        IntAsCurvePointAffine _zero;
};
IntAsCurve fakeCurve;

template <typename Curve>
class _BatchAccumulators: public BatchAccumulators<Curve> {
    public:
/*        using Stack<Curve>::getReferencesCount;
        using Stack<Curve>::set;
        using Stack<Curve>::dump;
        using Stack<Curve>::getOperationsLists;
        using Stack<Curve>::getReferenceIndex;
        using Stack<Curve>::dumpOperationsLists;
        using Stack<Curve>::clearEvaluated;    
        using Stack<Curve>::normalizeOperationsLists;*/
        _BatchAccumulators():BatchAccumulators<Curve>(fakeCurve){};
        void _addPoint(int64_t accId, int64_t value) {
            IntAsCurvePointAffine v(value);
            BatchAccumulators<Curve>::addPoint(accId, v);
        };
};


typedef _BatchAccumulators<IntAsCurve> BA;

TEST(batchOperation, basic) {
    BA ba;

    int64_t r1 = ba.defineAccumulators(100);
    int64_t r2 = ba.defineAccumulators(150);
    int64_t r3 = ba.defineAccumulators(150);

    ASSERT_EQ(100, r2 - r1);
    ASSERT_EQ(250, r3 - r1);
}
/*
TEST(batchOperation, singleAccumulator) {

    int values[] = {2, 3, 5, 7, 11, 13, 17, 19};
    int results[] = {0, 2, 5, 10, 17, 28, 41, 58, 77};
    int count = 8;

    // ba.calculate();

    // ASSERT_EQ((int)ba.getValue(0), 0);

    for (int n = 0; n <= count; ++n) {
        BA ba;
        ba.defineAccumulators(1);
        ba.setup(100, 100);
        std::cout << "\n\n\x1b[32mSTART (n:" << n << ")\x1b[0m\n";
        for (int i = 0; i < n; ++i) {
            std::cout << "\x1b[32mADD(0," << values[i] << ")\x1b[0m\n";
            ba._addPoint(0, values[i]);        
        }
        ba.calculate();
        auto result = ba.getValue(0);
        std::cout << "\x1b[32mRESULT:" << result << "\x1b[0m\n";
        ASSERT_EQ(results[n], result);
    }

    BA ba;
    ba.defineAccumulators(1);
    ba.setup(100, 100);
    for (int i = 0; i < count; ++i) {
        std::cout << "\x1b[32mADD(0," << values[i] << ")\x1b[0m\n";
        ba._addPoint(0, values[i]);        
        ba.calculate();
        auto result = ba.getValue(0);
        std::cout << "\x1b[32mRESULT:" << result << "\x1b[0m\n";
        ASSERT_EQ(results[i+1], result);
    }
*/
/*    ba._addPoint(0, 2);
    ba.calculate();
    std::cout << ba.getValue(0) << "\n";
    ASSERT_TRUE(fakeCurve.eq(ba.getValue(0), 2));


    ba._addPoint(0, 3);
    ba._addPoint(0, 5);
    ba._addPoint(0, 7);
    ba._addPoint(0, 11);
    ba._addPoint(0, 13);

    ba.calculate();
    std::cout << ba.getValue(0) << "\n";

    ASSERT_TRUE(fakeCurve.eq(ba.getValue(0), 10));*/
}
/*
TEST(batchOperation, sets) {
    StackExposed st;

    Reference r1 = st.define("r1", 2);
    st.set(r1+0, 2);
    st.set(r1+1, 4);
    ASSERT_TRUE(fakeCurve.eq(2, st.resolve(r1+0)));
    ASSERT_TRUE(fakeCurve.eq(4, st.resolve(r1+1)));
}
/*
TEST(batchOperation, add) {
    StackExposed st;

    Reference r1 = st.alloc(2);
    st.set(r1+0, 2);
    st.set(r1+1, 4);
    st.add(r1+0, r1+0, r1+1);
    st.copy(r1+1, r1+0);
    st.add(r1+0, r1+0, r1+0);
    st.dump();
    printf("r1+0 = %s\n", st.getPathToResolve(r1+0).c_str());
    printf("r1+1 = %s\n", st.getPathToResolve(r1+1).c_str());
    ASSERT_EQ(12, st.resolve(r1+0, 5));
    ASSERT_EQ(6, st.resolve(r1+1, 5));
}

TEST(batchOperation, add3) {
    StackExposed st;

    Reference r1 = st.alloc(5);
    auto a = r1+0;
    auto b = r1+1;
    auto c = r1+2;
    auto d = r1+3;
    auto aux = r1+4;
    st.set(a, 2);
    st.set(b, 4);
    st.set(c, 8);
    st.set(d, 16);
    st.set(aux, 10000);
    st.add(aux, a, b);
    st.add(d, aux, c);
    st.dump();
    printf("a = %s\n", st.getPathToResolve(a).c_str());
    printf("b = %s\n", st.getPathToResolve(b).c_str());
    printf("c = %s\n", st.getPathToResolve(c).c_str());
    printf("d = %s\n", st.getPathToResolve(d).c_str());
    printf("aux = %s\n", st.getPathToResolve(aux).c_str());
    ASSERT_EQ(2, st.resolve(a, 5));
    ASSERT_EQ(4, st.resolve(b, 5));
    ASSERT_EQ(8, st.resolve(c, 5));
    ASSERT_EQ(14, st.resolve(d, 5));
    ASSERT_EQ(6, st.resolve(aux, 5));
    OperationsLists lists;
    lists.add();
    st.clearEvaluated();
    st.getOperationsLists(lists, st.getReferenceIndex(d), 1000000),
    st.dumpOperationsLists(lists);
    st.dump();
}


TEST(batchOperation, add4) {
    StackExposed st;

    Reference r1 = st.alloc(5);
    auto a = r1+0;
    auto b = r1+1;
    auto c = r1+2;
    auto d = r1+3;
    auto aux = r1+4;
    st.set(a, 2);
    st.set(b, 4);
    st.set(c, 8);
    st.set(d, 16);
    st.set(aux, 10000);
    st.add(aux, a, b);
    st.add(d, aux, c);
    st.add(c, aux, d);
    st.dump();
    printf("a = %s\n", st.getPathToResolve(a).c_str());
    printf("b = %s\n", st.getPathToResolve(b).c_str());
    printf("c = %s\n", st.getPathToResolve(c).c_str());
    printf("d = %s\n", st.getPathToResolve(d).c_str());
    printf("aux = %s\n", st.getPathToResolve(aux).c_str());
    OperationsLists lists;
    lists.add();
    st.clearEvaluated();
    st.getOperationsLists(lists, st.getReferenceIndex(d), 1000000);
    st.dumpOperationsLists(lists);
    st.dump();
    lists.add();
    st.getOperationsLists(lists, st.getReferenceIndex(c), 1000000);
    st.normalizeOperationsLists(lists);
    st.dumpOperationsLists(lists);
    st.dump();
}


}  // namespace
*/
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
