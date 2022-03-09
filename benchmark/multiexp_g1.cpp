#ifndef COUNT_OPS
#define COUNT_OPS
#endif

#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include "alt_bn128.hpp"
#include <time.h>

using namespace AltBn128;

__uint128_t g_lehmer64_state = 0xAAAAAAAAAAAAAAAALL;

// Fast random generator
// https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/

uint64_t lehmer64() {
  g_lehmer64_state *= 0xda942042e4dd58b5LL;
  return g_lehmer64_state >> 64;
}

#define EXPECTED_RESULT_1000000_32 "(7686866163780120756504704687108787598650652185649163569056142218702518519446,3906118672583014628968951493493328376867199397126296469652067564264383995501)"
#define EXPECTED_RESULT_100000_32 "(4025341918271974808785600589734026802648663429991646260424633790951720792181,12092188261120553063009105063256001824310699929784774108735665676444962251470)"
#define EXPECTED_RESULT_100_32 "(2918313335010035498813915079123510985165681786090089918491607260013248937162,17618193400648240490544433854439160718346256852721163860896583680851810317370)"
#define EXPECTED_RESULT_10_32 "(9210363875285848131508656139203480860760390189260733424246762951552226162107,5835599852467557952169717222312357036255412228954185370294061936827699274723)"
#define EXPECTED_RESULT_2_32 "(11214335563238929817636333449692455681168983913051128543176696602685434833110,19682196534069234051611777524958224662412865849089027703776500470935035076258)"
#define EXPECTED_RESULT_2_1 "(18094853483684401030214035454673779939429883778299930824975271485319162037351,13688288498368359111250679033419188395420009978839941644881073507868645204619)"
#define EXPECTED_RESULT_2_2 "(17717503517063967127010656591622579159521290730611161775508432740108962728275,8347292658269310068911278187690526007547188570782947621250465360037527320812)"
#define EXPECTED_RESULT EXPECTED_RESULT_100000_32
/*
void generateBigDataFile ( uint64_t n )
{
    int64_t maxBlock = 100000;
    uint8_t *scalars = new uint8_t[maxBlock * 32];

    int64_t count = n;
    while (count > 0) {
        int
        for (int i=0; i<N*4; i++) {
            *((uint64_t *)(scalars + i*8)) = lehmer64();
        }


    printf("%ld %ld\n", read(fd, scalars, sizeof(scalars[0]) * N * 32), sizeof(scalars[0]) * N * 32);
    printf("%ld %ld\n",read(fd, bases, sizeof(bases[0]) * N), sizeof(bases[0]) * N);
    close(fd);
    int nscalars = 32;
    G1PointAffine *bases = new G1PointAffine[];

    N = 500000; // 0000;
    printf("N=%d\n", N);

    // random scalars
    for (int i=0; i<N*4; i++) {
        *((uint64_t *)(scalars + i*8)) = lehmer64();
    }


}
*/
int main(int argc, char **argv) 
{
//     int N = atoi(argv[1]);
    setlocale(LC_ALL, "en_US.utf-8");

    int N;
    int fd = open("testdata.dat", 0666);
    printf("%ld\n",read(fd, &N, sizeof(N)));
    printf("file contents %d bases\n", N);

    uint8_t *scalars = new uint8_t[N*32];
    G1PointAffine *bases = new G1PointAffine[N];

    printf("%ld %ld\n", read(fd, scalars, sizeof(scalars[0]) * N * 32), sizeof(scalars[0]) * N * 32);
    printf("%ld %ld\n",read(fd, bases, sizeof(bases[0]) * N), sizeof(bases[0]) * N);
    close(fd);
    int nscalars = 32;

    N = 100; // 0000;
    printf("N=%d\n", N);
/*
    // random scalars
    for (int i=0; i<N*4; i++) {
        *((uint64_t *)(scalars + i*8)) = lehmer64();
    }

    G1.copy(bases[0], G1.one());
    G1.copy(bases[1], G1.one());

    for (int i=2; i<N; i++) {
        G1.add(bases[i], bases[i-1], bases[i-2]);
    }
*/    
/*
    int fd = creat("testdata.dat", 0666);
    write(fd, &N, sizeof(N));
    write(fd, scalars, sizeof(scalars[0]) * N * 32);
    write(fd, bases, sizeof(bases[0]) * N);
    close(fd);
*/
    clock_t start, end;
    double cpu_time_used;
    std::string strResult;
/*
    G1Point p1;
    
    printf("Starting multiexp. (original) \n");
    G1.resetCounters();
    start = clock();
    G1.multiMulByScalarOriginal(p1, bases, (uint8_t *)scalars, nscalars, N);
    end = clock();
    strResult = G1.toString(p1); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time used: %.2lf\n", cpu_time_used);
    printf("Avg time per exp: %.2lf us\n", (cpu_time_used*1000000)/N);
    printf("Exps per second: %.2lf\n", (N / cpu_time_used));
*/
    G1Point p2;

    printf("Starting multiexp. (mix) \n");
    G1.resetCounters();
    start = clock();
    G1.multiMulByScalarMix(p2, bases, (uint8_t *)scalars, nscalars, N);
    end = clock();
    strResult = G1.toString(p2); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());
    int64_t total = G1.cntAddMixed + G1.cntAdd + G1.cntAddAffine;
    printf("P1 cntAddTotal:%'ld\n", total);

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time used: %.2lf\n", cpu_time_used);
    printf("Avg time per exp: %.2lf us\n", (cpu_time_used*1000000)/N);
    printf("Exps per second: %.2lf\n", (N / cpu_time_used));

    G1Point p3;

    printf("Starting multiexp. \n");
    G1.resetCounters();
    start = clock();
    G1.multiMulByScalar(p3, bases, (uint8_t *)scalars, nscalars, N);
    end = clock();
    strResult = G1.toString(p3); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());
    int64_t total2 = G1.cntAddMixed + G1.cntAdd + G1.cntAddAffine;
    printf("P1 cntAddTotal:%'ld\n", total2);

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time used: %.2lf\n", cpu_time_used);
    printf("Avg time per exp: %.2lf us\n", (cpu_time_used*1000000)/N);
    printf("Exps per second: %.2lf\n", (N / cpu_time_used));
    
    printf("P1 Better %.02f%%\n", ((double)(total-total2)*100.0)/total);
}