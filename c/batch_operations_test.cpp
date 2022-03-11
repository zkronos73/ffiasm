#include <gmp.h>
#include <iostream>

#include "gtest/gtest.h"
#include "batch_operations.hpp"

using namespace BatchOperations;

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
        const IntAsCurvePointAffine &zero (void) { return _zero; };
        bool eq (PointAffine op1, const PointAffine op2) { return op1.value == op2.value; };
    protected:
        IntAsCurvePointAffine _zero;
};
IntAsCurve fakeCurve;


template <typename Curve>
class _StackExposed: public Stack<Curve> {
    public:
        using Stack<Curve>::getReferencesCount;
        using Stack<Curve>::set;
        using Stack<Curve>::dump;
        using Stack<Curve>::getOperationsLists;
        using Stack<Curve>::getReferenceIndex;
        using Stack<Curve>::dumpOperationsLists;
        using Stack<Curve>::clearEvaluated;    
        using Stack<Curve>::normalizeOperationsLists;
        _StackExposed():Stack<Curve>(fakeCurve){};
};


typedef _StackExposed<IntAsCurve> StackExposed;
TEST(batchOperation, basic) {
    // StackExposed<IntAsCurve> st(fakeCurve);
    StackExposed st;

    Reference r1 = st.define("r1", 100);
    Reference r2 = st.define("r2", 150);
    Reference r3 = st.define("r3", 150);

    ASSERT_EQ(100, r2 - r1);
    ASSERT_EQ(250, r3 - r1);
}

TEST(batchOperation, zeroRefenceCount) {
    StackExposed st;

    int zeroRefCount1 = st.getReferencesCount(0);
    Reference r1 = st.define("r1", 100);
    int zeroRefCount2 = st.getReferencesCount(0);
    ASSERT_EQ(100, zeroRefCount2 - zeroRefCount1);    
}

TEST(batchOperation, sets) {
    StackExposed st;

    Reference r1 = st.define("r1", 2);
    Reference r2 = st.define("r2", 4);
    st.set(r1+0, 2);
    st.set(r1+1, 3);
    st.set(r2+0, 5);
    st.set(r2+1, 7);
    st.set(r2+2, 11);
    st.set(r2+3, 13);

    st.add(r1+0, r1+0, r2+0);
    st.add(r1+0, r1+0, r2+1);
    st.add(r1+0, r1+0, r2+2);
    st.add(r1+0, r1+0, r2+1);
    st.add(r1+0, r1+0, r2+3);

    st.add(r1+1, r1+1, r2+3);
    st.add(r1+1, r1+1, r2+2);
    st.dbl(r1+1, r1+1);
    st.dbl(r1+1, r1+1);
    st.dbl(r1+1, r1+1);
    st.dbl(r1+1, r1+1);
    st.add(r1+1, r1+0, r1+1);
    st.dbl(r1+0, r1+1);

    std::cout << st.getStringAdditionValuesOfReference(r1+0) << "\n";
    std::cout << st.getStringAdditionValuesOfReference(r1+1) << "\n";
//     ASSERT_TRUE(fakeCurve.eq(2, st.resolve(r1+0)));
//    ASSERT_TRUE(fakeCurve.eq(4, st.resolve(r1+1)));
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
*/
/*
TEST(altBn128, g1_PlusZero) {
    G1Point p1;

    G1.add(p1, G1.one(), G1.zero());

    ASSERT_TRUE(G1.eq(p1, G1.one()));
}

TEST(altBn128, g1_minus_g1) {
    G1Point p1;

    G1.sub(p1, G1.one(), G1.one());

    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g1_times_4) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);

    ASSERT_TRUE(G1.eq(p1,p2));
}


TEST(altBn128, g1_times_3) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);
    G1.sub(p2, p2, G1.one());

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_3_exp) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    mpz_t e;
    mpz_init_set_str(e, "3", 10);

    uint8_t scalar[32];
    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_5) {
    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.add(p1, p1, p1);

    G1Point p2;
    G1Point p3;
    G1Point p4;
    G1Point p5;
    G1Point p6;
    G1.dbl(p2, G1.one());
    G1.dbl(p3, p2);
    G1.dbl(p4, G1.one());
    G1.dbl(p5, p4);
    G1.add(p6, p3, p5);

    ASSERT_TRUE(G1.eq(p1,p6));
}

TEST(altBn128, g1_times_65_exp) {

    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.add(p1, p1, G1.one());

    mpz_t e;
    mpz_init_set_str(e, "65", 10);

    uint8_t scalar[32];
    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_expToOrder) {
    mpz_t e;
    mpz_init_set_str(e, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    uint8_t scalar[32];

    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p1;

    G1.mulByScalar(p1, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g2_expToOrder) {
    mpz_t e;
    mpz_init_set_str(e, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    uint8_t scalar[32];

    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    Curve<F2Field<RawFq>>::Point p1;

    G2.mulByScalar(p1, G2.one(), scalar, 32);

    ASSERT_TRUE(G2.isZero(p1));
}

TEST(altBn128, multiExp) {

    int NMExp = 40000;

    typedef uint8_t Scalar[32];

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    uint64_t acc=0;
    for (int i=0; i<NMExp; i++) {
        if (i==0) {
            G1.copy(bases[0], G1.one());
        } else {
            G1.add(bases[i], bases[i-1], G1.one());
        }
        for (int j=0; j<32; j++) scalars[i][j] = 0;
        *(int *)&scalars[i][0] = i+1;
        acc += (i+1)*(i+1);
    }

    G1Point p1;
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, 32, NMExp);

    mpz_t e;
    mpz_init_set_ui(e, acc);

    Scalar sAcc;

    for (int i=0;i<32;i++) sAcc[i] = 0;
    mpz_export((void *)sAcc, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), sAcc, 32);

    ASSERT_TRUE(G1.eq(p1, p2));

    delete[] bases;
    delete[] scalars;
}


TEST(altBn128, multiExp2) {

    int NMExp = 2;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1626275109576878988287730541908027724405348106427831594181487487855202143055");
    F1.fromString(bases[0].y, "18706364085805828895917702468512381358405767972162700276238017959231481018884");
    F1.fromString(bases[1].x, "17245156998235704504461341147511350131061011207199931581281143511105381019978");
    F1.fromString(bases[1].y, "3858908536032228066651712470282632925312300188207189106507111128103204506804");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "20187316456970436521602619671088988952475789765726813868033071292105413408473");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "9163953212624378696742080269971059027061360176019470242548968584908855004282");
    F1.fromString(ref.y, "20922060990592511838374895951081914567856345629513259026540392951012456141360");

    G1.multiMulByScalar(r, bases, (uint8_t *)scalars, 32, 2);
    G1.copy(ra, r);

    // std::cout << G1.toString(r, 10);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, fft) {
    int NMExp = 1<<10;

    AltBn128::FrElement *a = new AltBn128::FrElement[NMExp];

    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(a[i], i+1);
    }

    FFT<typename Engine::Fr> fft(NMExp);

    fft.fft(a, NMExp);
    fft.ifft(a, NMExp);

    AltBn128::FrElement aux;
    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(aux, i+1);
        ASSERT_TRUE(Fr.eq(a[i], aux));
    }

    delete[] a;
}
*/

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
