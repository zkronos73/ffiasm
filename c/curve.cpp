#include <sstream>

template <typename BaseField>
Curve<BaseField>::Curve(BaseField &aF, typename BaseField::Element &aa, typename BaseField::Element &ab, typename BaseField::Element &agx, typename BaseField::Element &agy) : F(aF) {
    F = aF;
    initCurve(aa, ab, agx, agy);
}

template <typename BaseField>
Curve<BaseField>::Curve(BaseField &aF, std::string as, std::string bs, std::string gxs, std::string gys) : F(aF) {
    F = aF;

    typename BaseField::Element aa;
    typename BaseField::Element ab;
    typename BaseField::Element agx;
    typename BaseField::Element agy;

    F.fromString(aa, as);
    F.fromString(ab, bs);
    F.fromString(agx, gxs);
    F.fromString(agy, gys);

    initCurve(aa, ab, agx, agy);
}


template <typename BaseField>
void Curve<BaseField>::initCurve(typename BaseField::Element &aa, typename BaseField::Element &ab, typename BaseField::Element &agx, typename BaseField::Element &agy) {
    F.copy(fa, aa);
    F.copy(fb, ab);
    F.copy(fone.x, agx);
    F.copy(foneAffine.x, agx);
    F.copy(fone.y, agy);
    F.copy(foneAffine.y, agy);
    F.copy(fone.zz, F.one());
    F.copy(fone.zzz, F.one());
    F.copy(fzero.x, F.one());
    F.copy(fzeroAffine.x, F.zero());
    F.copy(fzero.y, F.one());
    F.copy(fzeroAffine.y, F.zero());
    F.copy(fzero.zz, F.zero());
    F.copy(fzero.zzz, F.zero());

    if (F.isZero(fa)) {
        typeOfA = a_is_zero;
    } else if (F.eq(fa, F.one())) {
        typeOfA = a_is_one;
    } else if (F.eq(fa, F.negOne())) {
        typeOfA = a_is_negone;
    } else {
        typeOfA = a_is_long;
    }

#ifdef COUNT_OPS
    resetCounters();
#endif // COUNT_OPS

}

template <typename BaseField>
void inline Curve<BaseField>::mulByA(typename BaseField::Element &r, const typename BaseField::Element &ab) {
    switch (typeOfA) {
        case a_is_zero: F.copy(r, F.zero()); break;
        case a_is_one: F.copy(r, ab); break;
        case a_is_negone: F.neg(r, ab); break;
        case a_is_long: F.mul(r, fa, ab);
    }
}



/*
    https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-add-2008-s
    U1 = X1*ZZ2
    U2 = X2*ZZ1
    S1 = Y1*ZZZ2
    S2 = Y2*ZZZ1
    P = U2-U1
    R = S2-S1
    PP = P^2
    PPP = P*PP
    Q = U1*PP
    X3 = R^2-PPP-2*Q
    Y3 = R*(Q-X3)-S1*PPP
    ZZ3 = ZZ1*ZZ2*PP
    ZZZ3 = ZZZ1*ZZZ2*PPP
*/
template <typename BaseField>
void Curve<BaseField>::add(Point &p3, const Point &p1, const Point &p2) {
#ifdef COUNT_OPS
    cntAdd++;
#endif // COUNT_OPS

    if (isZero(p1)) {
        copy(p3, p2);
        return;
    }

    if (isZero(p2)) {
        copy(p3, p1);
        return;
    }

    typename BaseField::Element tmp;

    // U1 = X1*ZZ2
    typename BaseField::Element U1;
    F.mul(U1, p1.x, p2.zz);

    // U2 = X2*ZZ1
    typename BaseField::Element U2;
    F.mul(U2, p2.x, p1.zz);

    // S1 = Y1*ZZZ2
    typename BaseField::Element S1;
    F.mul(S1, p1.y, p2.zzz);

    // S2 = Y2*ZZZ1
    typename BaseField::Element S2;
    F.mul(S2, p2.y, p1.zzz);

    // P = U2-U1
    typename BaseField::Element P;
    F.sub(P, U2, U1);

    // R = S2-S1
    typename BaseField::Element R;
    F.sub(R, S2, S1);

    if (F.isZero(P) && F.isZero(R)) return dbl(p3, p1);

    // PP = P^2
    typename BaseField::Element PP;
    F.square(PP, P);

    // PPP = P*PP
    typename BaseField::Element PPP;
    F.mul(PPP, P, PP);

    // Q = U1*PP
    typename BaseField::Element Q;
    F.mul(Q, U1, PP);

    // X3 = R^2-PPP-2*Q
    F.square(p3.x, R);
    F.sub(p3.x, p3.x, PPP);
    F.sub(p3.x, p3.x, Q);
    F.sub(p3.x, p3.x, Q);

    // Y3 = R*(Q-X3)-S1*PPP
    F.mul(tmp, S1, PPP);
    F.sub(p3.y, Q, p3.x);
    F.mul(p3.y, p3.y, R );
    F.sub(p3.y, p3.y, tmp);

    // ZZ3 = ZZ1*ZZ2*PP
    F.mul(p3.zz, p1.zz, p2.zz);
    F.mul(p3.zz, p3.zz, PP);
    
    // ZZZ3 = ZZZ1*ZZZ2*PPP
    F.mul(p3.zzz, p1.zzz, p2.zzz);
    F.mul(p3.zzz, p3.zzz, PPP);

}


/*
    https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
    U2 = X2*ZZ1
    S2 = Y2*ZZZ1
    P = U2-X1
    R = S2-Y1
    PP = P^2
    PPP = P*PP
    Q = X1*PP
    X3 = R^2-PPP-2*Q
    Y3 = R*(Q-X3)-Y1*PPP
    ZZ3 = ZZ1*PP
    ZZZ3 = ZZZ1*PPP
*/

template <typename BaseField>
void Curve<BaseField>::add(Point &p3, const Point &p1, const PointAffine &p2) {
#ifdef COUNT_OPS
    cntAddMixed++;
#endif // COUNT_OPS


    if (isZero(p1)) {
        copy(p3, p2);
        return;
    }

    if (isZero(p2)) {
        copy(p3, p1);
        return;
    }

    typename BaseField::Element tmp;

    // U2 = X2*ZZ1
    typename BaseField::Element U2;
    F.mul(U2, p2.x, p1.zz);

    // S2 = Y2*ZZZ1
    typename BaseField::Element S2;
    F.mul(S2, p2.y, p1.zzz);

    // P = U2-X1
    typename BaseField::Element P;
    F.sub(P, U2, p1.x);

    // R = S2-Y1
    typename BaseField::Element R;
    F.sub(R, S2, p1.y);

    if (F.isZero(P) && F.isZero(R)) return dbl(p3, p2);

    // PP = P^2
    typename BaseField::Element PP;
    F.square(PP, P);

    // PPP = P*PP
    typename BaseField::Element PPP;
    F.mul(PPP, P, PP);

    // Q = X1*PP
    typename BaseField::Element Q;
    F.mul(Q, p1.x, PP);

    // X3 = R^2-PPP-2*Q
    F.square(p3.x, R);
    F.sub(p3.x, p3.x, PPP);
    F.sub(p3.x, p3.x, Q);
    F.sub(p3.x, p3.x, Q);

    // Y3 = R*(Q-X3)-Y1*PPP
    F.mul(tmp, p1.y, PPP);
    F.sub(p3.y, Q, p3.x);
    F.mul(p3.y, p3.y, R );
    F.sub(p3.y, p3.y, tmp);

    // ZZ3 = ZZ1*PP
    F.mul(p3.zz, p1.zz, PP);
    
    // ZZZ3 = ZZZ1*PPP
    F.mul(p3.zzz, p1.zzz, PPP);
}


/*
    https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
    P = X2-X1
    R = Y2-Y1
    PP = P^2
    PPP = P*PP
    Q = X1*PP
    X3 = R^2-PPP-2*Q
    Y3 = R*(Q-X3)-Y1*PPP
    ZZ3 = PP
    ZZZ3 = PPP
*/
template <typename BaseField>
void Curve<BaseField>::add(Point &p3, const PointAffine &p1, const PointAffine &p2) {
#ifdef COUNT_OPS
    cntAddAffine++;
#endif // COUNT_OPS

    if (isZero(p1)) {
        copy(p3, p2);
        return;
    }

    if (isZero(p2)) {
        copy(p3, p1);
        return;
    }

    typename BaseField::Element tmp;

    // P = X2-X1
    typename BaseField::Element P;
    F.sub(P, p2.x, p1.x);

    // R = Y2-Y1
    typename BaseField::Element R;
    F.sub(R, p2.y, p1.y);

    if (F.isZero(P) && F.isZero(R)) return dbl(p3, p2);

    // PP = P^2
    typename BaseField::Element PP;
    F.square(PP, P);

    // PPP = P*PP
    typename BaseField::Element PPP;
    F.mul(PPP, P, PP);

    // Q = X1*PP
    typename BaseField::Element Q;
    F.mul(Q, p1.x, PP);

    // X3 = R^2-PPP-2*Q
    F.square(p3.x, R);
    F.sub(p3.x, p3.x, PPP);
    F.sub(p3.x, p3.x, Q);
    F.sub(p3.x, p3.x, Q);

    // Y3 = R*(Q-X3)-Y1*PPP
    F.mul(tmp, p1.y, PPP);
    F.sub(p3.y, Q, p3.x);
    F.mul(p3.y, p3.y, R );
    F.sub(p3.y, p3.y, tmp);

    // ZZ3 = PP
    F.copy(p3.zz, PP);
    
    // ZZZ3 = PPP
    F.copy(p3.zzz, PPP);
}



/*
    https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
    U = 2*Y1
    V = U^2
    W = U*V
    S = X1*V
    M = 3*X1^2+a*ZZ1^2
    X3 = M^2-2*S
    Y3 = M*(S-X3)-W*Y1
    ZZ3 = V*ZZ1
    ZZZ3 = W*ZZZ1

*/
template <typename BaseField>
void Curve<BaseField>::dbl(Point &p3, const Point &p1) {
#ifdef COUNT_OPS
    cntDbl++;
#endif // COUNT_OPS


    if (isZero(p1)) {
        copy(p3, p1);
        return;
    }

    typename BaseField::Element tmp;

    // U = 2*Y1
    typename BaseField::Element U;
    F.add(U, p1.y, p1.y);

    // V = U^2
    typename BaseField::Element V;
    F.square(V, U);

    // W = U*V
    typename BaseField::Element W;
    F.mul(W, U, V);

    // S = X1*V
    typename BaseField::Element S;
    F.mul(S, p1.x, V);

    //M = 3*X1^2+a*ZZ1^2
    typename BaseField::Element M;
    F.square(M, p1.x);
    F.add(tmp, M, M);
    F.add(M, M, tmp);
    if (typeOfA != a_is_zero) {
        F.square(tmp, p1.zz);
        mulByA(tmp, tmp);
        F.add(M, M, tmp);
    }

    // X3 = M^2-2*S
    F.square(p3.x, M);
    F.sub(p3.x, p3.x, S);
    F.sub(p3.x, p3.x, S);

    // Y3 = M*(S-X3)-W*Y1
    F.mul(tmp, W, p1.y);
    F.sub(p3.y, S, p3.x);
    F.mul(p3.y, M, p3.y);
    F.sub(p3.y, p3.y, tmp);

    // ZZ3 = V*ZZ1
    F.mul(p3.zz, V, p1.zz);

    // ZZZ3 = W*ZZZ1
    F.mul(p3.zzz, W, p1.zzz);
}

/*
    https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
    U = 2*Y1
    V = U^2
    W = U*V
    S = X1*V
    M = 3*X1^2+a
    X3 = M^2-2*S
    Y3 = M*(S-X3)-W*Y1
    ZZ3 = V
    ZZZ3 = W
*/
template <typename BaseField>
void Curve<BaseField>::dbl(Point &p3, const PointAffine &p1) {
#ifdef COUNT_OPS
    cntDblMixed++;
#endif // COUNT_OPS

    if (isZero(p1)) {
        copy(p3, p1);
        return;
    }

    typename BaseField::Element tmp;

    // U = 2*Y1
    typename BaseField::Element U;
    F.add(U, p1.y, p1.y);

    // V = U^2   ; Already store in ZZ3
    F.square(p3.zz, U);

    // W = U*V   ; Alreadu store in ZZZ3
    F.mul(p3.zzz, U, p3.zz);

    // S = X1*V
    typename BaseField::Element S;
    F.mul(S, p1.x, p3.zz);

    // M = 3*X1^2+a
    typename BaseField::Element M;
    F.square(M, p1.x);
    F.add(tmp, M, M);
    F.add(M, tmp, M);
    F.add(M, M, fa);

    // X3 = M^2-2*S
    F.square(p3.x, M);
    F.sub(p3.x, p3.x, S);
    F.sub(p3.x, p3.x, S);

    // Y3 = M*(S-X3)-W*Y1
    F.mul(tmp, p3.zzz, p1.y);
    F.sub(p3.y, S, p3.x);
    F.mul(p3.y, M, p3.y);
    F.sub(p3.y, p3.y, tmp);

    // ZZ3 = V ; Already stored

    // ZZZ3 = W ; Already stored
}     

template <typename BaseField>
bool Curve<BaseField>::eq(const Point &p1, const Point &p2) {
#ifdef COUNT_OPS
    cntEq++;
#endif // COUNT_OPS

    if (isZero(p1)) return  isZero(p2);

    // U1 = X1*ZZ2
    typename BaseField::Element U1;
    F.mul(U1, p1.x, p2.zz);

    // U2 = X2*ZZ1
    typename BaseField::Element U2;
    F.mul(U2, p2.x, p1.zz);

    // S1 = Y1*ZZZ2
    typename BaseField::Element S1;
    F.mul(S1, p1.y, p2.zzz);

    // S2 = Y2*ZZZ1
    typename BaseField::Element S2;
    F.mul(S2, p2.y, p1.zzz);

    // P = U2-U1
    typename BaseField::Element P;
    F.sub(P, U2, U1);

    // R = S2-S1
    typename BaseField::Element R;
    F.sub(R, S2, S1);

    return (F.isZero(P) && F.isZero(R));
}


template <typename BaseField>
bool Curve<BaseField>::eq(const Point &p1, const PointAffine &p2) {
#ifdef COUNT_OPS
    cntEqMixed++;
#endif // COUNT_OPS

    if (isZero(p1)) return  isZero(p2);

    typename BaseField::Element tmp;

    // U2 = X2*ZZ1
    typename BaseField::Element U2;
    F.mul(U2, p2.x, p1.zz);

    // S2 = Y2*ZZZ1
    typename BaseField::Element S2;
    F.mul(S2, p2.y, p1.zzz);

    // P = U2-X1
    typename BaseField::Element P;
    F.sub(P, U2, p1.x);

    // R = S2-Y1
    typename BaseField::Element R;
    F.sub(R, S2, p1.y);

    return (F.isZero(P) && F.isZero(R));
}

template <typename BaseField>
bool Curve<BaseField>::eq(const PointAffine &p1, const PointAffine &p2) {
    return F.eq(p1.x, p2.x) && F.eq(p1.y, p2.y);
}


template <typename BaseField>
bool Curve<BaseField>::isZero(const Point &p1) {
    return F.isZero(p1.zz);
}

template <typename BaseField>
bool Curve<BaseField>::isZero(const PointAffine &p1) {
    return F.isZero(p1.x) && F.isZero(p1.y);
}

template <typename BaseField>
void Curve<BaseField>::copy(Point &r, const Point &a) {
    F.copy(r.x, a.x);
    F.copy(r.y, a.y);
    F.copy(r.zz, a.zz);
    F.copy(r.zzz, a.zzz);
}

template <typename BaseField>
void Curve<BaseField>::copy(Point &r, const PointAffine &a) {
    if (isZero(a)) {
        F.copy(r.x, F.one());
        F.copy(r.y, F.one());
        F.copy(r.zz, F.zero());
        F.copy(r.zzz, F.zero());
        return;
    }
    F.copy(r.x, a.x);
    F.copy(r.y, a.y);
    F.copy(r.zz, F.one());
    F.copy(r.zzz, F.one());
}

template <typename BaseField>
void Curve<BaseField>::copy(PointAffine &r, const Point &a) {
#ifdef COUNT_OPS
    cntToAffine++;
#endif // COUNT_OPS
    if (isZero(a)) {
        F.copy(r.x, F.zero());
        F.copy(r.y, F.zero());
        return;
    }
    F.div(r.x, a.x, a.zz);
    F.div(r.y, a.y, a.zzz);
}

template <typename BaseField>
void Curve<BaseField>::copy(PointAffine &r, const PointAffine &a) {
    F.copy(r.x, a.x);
    F.copy(r.y, a.y);
}

template <typename BaseField>
void Curve<BaseField>::neg(Point &r, const Point &a) {
    F.copy(r.x, a.x);
    F.neg(r.y, a.y);
    F.copy(r.zz, a.zz);
    F.copy(r.zzz, a.zzz);
}

template <typename BaseField>
void Curve<BaseField>::neg(Point &r, const PointAffine &a) {
    F.copy(r.x, a.x);
    F.neg(r.y, a.y);
    F.copy(r.zz, F.one());
    F.copy(r.zzz, F.one());
}

template <typename BaseField>
void Curve<BaseField>::neg(PointAffine &r, const Point &a) {
#ifdef COUNT_OPS
    cntToAffine++;
#endif // COUNT_OPS
    if (isZero(a)) {
        F.copy(r.x, F.zero());
        F.copy(r.y, F.zero());
        return;
    }
    F.div(r.x, a.x, a.zz);
    F.div(r.y, a.y, a.zzz);
    F.neg(r.y, r.y);
}

template <typename BaseField>
void Curve<BaseField>::neg(PointAffine &r, const PointAffine &a) {
    F.copy(r.x, a.x);
    F.neg(r.y, a.y);
}


template <typename BaseField>
std::string Curve<BaseField>::toString(const Point &p, uint32_t radix) {
    PointAffine tmp;
    copy(tmp, p);
    std::ostringstream stringStream;
    stringStream << "(" << F.toString(tmp.x, radix) << "," << F.toString(tmp.y, radix) << ")";
    return stringStream.str();
}

template <typename BaseField>
std::string Curve<BaseField>::toString(const PointAffine &p, uint32_t radix) {
    std::ostringstream stringStream;
    stringStream << "(" << F.toString(p.x, radix) << "," << F.toString(p.y, radix) << ")";
    return stringStream.str();
}

#ifdef COUNT_OPS
template <typename BaseField>
void Curve<BaseField>::resetCounters() {
    cntAddMixed = 0;
    cntAdd=0;
    cntAddAffine=0;
    cntDbl=0;
    cntDblMixed=0;
    cntEq=0;
    cntEqMixed=0;
    cntToAffine=0;
}

template <typename BaseField>
void Curve<BaseField>::printCounters() {
    printf("cntAddMixed: %'d\n", cntAddMixed);
    printf("cntAdd: %'d\n", cntAdd);
    printf("cntAddAffine: %'d\n", cntAddAffine);
    printf("cntAdd TOTAL: %'d\n", cntAddMixed + cntAdd+ cntAddAffine);
    printf("cntDbl: %'d\n", cntDbl);
    printf("cntDblMixed: %'d\n", cntDblMixed);
    printf("cntDbl TOTAL: %'d\n", cntDbl + cntDblMixed);
    printf("cntEq: %'d\n", cntEq);
    printf("cntEqMixed: %'d\n", cntEqMixed);
    printf("cntToAffine: %'d\n", cntToAffine);    
}
#endif // COUNT_OPS

template <typename BaseField>
void Curve<BaseField>::multiAdd(PointAffine *p3, const PointAffine *p1, const PointAffine *p2, u_int64_t count) {
#ifdef COUNT_OPS
    cntAddAffine++;
#endif // COUNT_OPS
    const auto p3AlreadyCalculated = -1;
//    int64_t *eqs = new int64_t [count];
    typename BaseField::Element *lambdas = (typename BaseField::Element *)malloc(sizeof(typename BaseField::Element) * count);
    u_int64_t lambdaIndex = 0;

    for (auto index = 0; index < count; ++index) {
//        __builtin_prefetch(p1 + index + 8);
//        __builtin_prefetch(p2 + index + 8);
        // auto &_p1 = p1[index];
        // auto &_p2 = p2[index];
        // auto &_eqs = eqs[index];
/*        if (isZero(p1[index])) {
            copy(p3[index], p2[index]);
            eqs[index] = p3AlreadyCalculated;
            continue;
        }
        if (isZero(p2[index])) {
            copy(p3[index], p1[index]);
            eqs[index] = p3AlreadyCalculated;
            continue;
        }*/
/*        eqs[index] = F.eq(p1[index].x, p2[index].x) && F.eq(p1[index].y, p2[index].y);
        if (eqs[index]) {*/
        if (F.eq(p1[index].x, p2[index].x) && F.eq(p1[index].y, p2[index].y)) {
            F.add(lambdas[lambdaIndex++], p1[index].y, p1[index].y);
        }
        else {
            F.sub(lambdas[lambdaIndex++], p2[index].x, p1[index].x);
        }
        // F.copy(lambdas[lambdaIndex++], _eqs ? F.add(_p1.y, _p1.y) : F.sub(_p2.x,_p1.x));
    }

    auto lambdaCount = lambdaIndex;
    F.batchInverse(lambdas, lambdas, lambdaCount);
    lambdaIndex = 0;
    for (auto index = 0; index < count; ++index) {
//        __builtin_prefetch(p1 + index + 8);
//        __builtin_prefetch(p2 + index + 8);
//        __builtin_prefetch(eqs + index + 8);
//        __builtin_prefetch(lambdas + lambdaIndex + 8);
        /* auto &_p1 = p1[index];
        auto &_p2 = p2[index];
        auto &_p3 = p3[index];
        auto &_eqs = eqs[index];
        auto &_lambda = lambdas[lambdaIndex];*/

        // if (eqs[index] == p3AlreadyCalculated) continue;

//        if (eqs[index]) {            
        if (F.eq(p1[index].x, p2[index].x) && F.eq(p1[index].y, p2[index].y)) {
            // l = l * (3 * p1.x**2) + a
            F.mul(lambdas[lambdaIndex], lambdas[lambdaIndex], F.add(F.mul(F.square(p1[index].x), 3), fa));
        }
        else {
            // l = l * (p2.y - p1.y)
            F.mul(lambdas[lambdaIndex], lambdas[lambdaIndex], F.sub(p2[index].y, p1[index].y));
        }

        // p3.x = l**2 - (p1.x + p2.x) 
        F.sub(p3[index].x, F.square(lambdas[lambdaIndex]), F.add(p1[index].x, p2[index].x));

        // p3.y = l * (p1.x - p3.x) - p1.y
        F.sub(p3[index].y, F.mul(lambdas[lambdaIndex], F.sub(p1[index].x, p3[index].x)), p1[index].y);
        
        ++lambdaIndex;
    }
//    free(eqs);
    free(lambdas);
}

#define P1(X) p[X].left
#define P2(X) p[X].right
#define P3(X) p[X].result
// #define EQS(X) p[X].eqs
#define EQS(X) eqs[X]
// #define LAMBDA(X) p[X].inverse.lambda
#define LAMBDA(X) lambdas[X]

template <typename BaseField>
void Curve<BaseField>::multiAdd2(AddPointAffine *p, u_int64_t count) {
#ifdef COUNT_OPS
    cntAddAffine++;
#endif // COUNT_OPS
    const auto p3AlreadyCalculated = -1;
    u_int64_t lambdaIndex = 0;
    int64_t *eqs = (int64_t *)malloc(sizeof(int64_t) * count);
    typename BaseField::Element *lambdas = (typename BaseField::Element *)malloc(sizeof(typename BaseField::Element) * count);

    for (auto index = 0; index < count; ++index) {
//        __builtin_prefetch(p1 + index + 8);
//        __builtin_prefetch(p2 + index + 8);
        // auto &_p1 = P1(index);
        // auto &_p2 = P2(index);
        // auto &_eqs = EQS(index);
        if (isZero(P1(index))) {
            copy(P3(index), P2(index));
            EQS(index) = p3AlreadyCalculated;
            continue;
        }
        if (isZero(P2(index))) {
            copy(P3(index), P1(index));
            EQS(index) = p3AlreadyCalculated;
            continue;
        }
        EQS(index) = F.eq(P1(index).x, P2(index).x) && F.eq(P1(index).y, P2(index).y);
        if (EQS(index)) {
            F.add(LAMBDA(lambdaIndex++), P1(index).y, P1(index).y);
        }
        else {
            F.sub(LAMBDA(lambdaIndex++), P2(index).x, P1(index).x);
        }
        // F.copy(LAMBDA(lambdaIndex++), _eqs ? F.add(_p1.y, _p1.y) : F.sub(_p2.x,_p1.x));
    }

    auto lambdaCount = lambdaIndex;
//     F.batchInverse(&(p->inverse), sizeof(p[0]), lambdaCount);
    F.batchInverse(lambdas, lambdas, lambdaCount);
    lambdaIndex = 0;
    for (auto index = 0; index < count; ++index) {
//        __builtin_prefetch(p1 + index + 8);
//        __builtin_prefetch(p2 + index + 8);
//        __builtin_prefetch(eqs + index + 8);
//        __builtin_prefetch(lambdas + lambdaIndex + 8);
        /* auto &_p1 = P1(index);
        auto &_p2 = P2(index);
        auto &_p3 = P3(index);
        auto &_eqs = EQS(index);
        auto &_lambda = LAMBDA(lambdaIndex);*/

        if (EQS(index) == p3AlreadyCalculated) continue;

        if (EQS(index)) {            
            // l = l * (3 * p1.x**2) + a
            F.mul(LAMBDA(lambdaIndex), LAMBDA(lambdaIndex), F.add(F.mul(F.square(P1(index).x), 3), fa));
        }
        else {
            // l = l * (p2.y - p1.y)
            F.mul(LAMBDA(lambdaIndex), LAMBDA(lambdaIndex), F.sub(P2(index).y, P1(index).y));
        }

        // p3.x = l**2 - (p1.x + p2.x) 
        F.sub(P3(index).x, F.square(LAMBDA(lambdaIndex)), F.add(P1(index).x, P2(index).x));

        // p3.y = l * (p1.x - p3.x) - p1.y
        F.sub(P3(index).y, F.mul(LAMBDA(lambdaIndex), F.sub(P1(index).x, P3(index).x)), P1(index).y);
        
        ++lambdaIndex;
    }
    free(lambdas);
    free(eqs);
}




