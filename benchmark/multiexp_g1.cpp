#ifndef COUNT_OPS
#define COUNT_OPS
#endif

#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
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

#define EXPECTED_RESULT "(7686866163780120756504704687108787598650652185649163569056142218702518519446,3906118672583014628968951493493328376867199397126296469652067564264383995501)"

int main(int argc, char **argv) 
{
//     int N = atoi(argv[1]);
    int N;
    int fd = open("testdata.dat", 0666);
    printf("%ld\n",read(fd, &N, sizeof(N)));
    printf("N=%d\n", N);

    uint8_t *scalars = new uint8_t[N*32];
    G1PointAffine *bases = new G1PointAffine[N];

    printf("%ld %ld\n", read(fd, scalars, sizeof(scalars[0]) * N * 32), sizeof(scalars[0]) * N * 32);
    printf("%ld %ld\n",read(fd, bases, sizeof(bases[0]) * N), sizeof(bases[0]) * N);
    close(fd);

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

    G1Point p1;
    
    printf("Starting multiexp. \n");
    G1.resetCounters();
    start = clock();
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, 32, N, 1);
    end = clock();
    std::string strResult = G1.toString(p1); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time used: %.2lf\n", cpu_time_used);
    printf("Avg time per exp: %.2lf us\n", (cpu_time_used*1000000)/N);
    printf("Exps per second: %.2lf\n", (N / cpu_time_used));

}