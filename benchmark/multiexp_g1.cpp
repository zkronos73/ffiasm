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
#define EXPECTED_RESULT EXPECTED_RESULT_1000000_32

inline int64_t getRealTimeClockUs ( void )
{
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    return tm.tv_sec * 1000000 + tm.tv_nsec / 1000;
}

void createDataFile ( const std::string &filename, int64_t _n )
{
    int64_t bytes, bytesWritten;
    int64_t n = _n;

    int fd = creat(filename.c_str(), 0666);
    if (write(fd, &n, sizeof(n)) != sizeof(n)) {
        printf("ERROR writting header E: %d %s\n", errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    int64_t block = 10240;
    uint8_t *scalars = new uint8_t[block * 32];
    int64_t count = n*4;
    while (count > 0) {
        int64_t cblock = count > block ? block : count;
        for (int i=0; i<cblock; i++) {
            *((uint64_t *)(scalars + i*8)) = lehmer64();
        }
        bytes = cblock * 8;
        bytesWritten = write(fd, scalars, bytes);
        if (bytes != bytesWritten) {
            printf("ERROR writting %ld bytes (scalars), write only %ld bytes. E:%d %s\n", bytes, bytesWritten, errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
        count -= cblock;
    }

    G1PointAffine *bases = new G1PointAffine[block+2];

    G1.copy(bases[0], G1.one());
    G1.copy(bases[1], G1.one());
    

    count = n;
    while (count > 0) {
        int64_t cblock = count > block ? block : count;
        for (int i=2; i<(cblock + 2); i++) {
            G1.add(bases[i], bases[i-1], bases[i-2]);
        }
        bytes = cblock * sizeof(bases[0]);
        bytesWritten = write(fd, bases, bytes);
        if (bytes != bytesWritten) {
            printf("ERROR writting %ld bytes (bases), write only %ld bytes. E:%d %s\n", bytes, bytesWritten, errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
        G1.copy(bases[0], bases[cblock]);
        G1.copy(bases[1], bases[cblock+1]);
    
        count -= cblock;
    }

    close(fd);
}

void loadDataFile ( const std::string &filename, uint8_t *scalars, G1PointAffine *bases, int64_t _n )
{
    int64_t n;
    int fd = open(filename.c_str(), 0666);
    if (read(fd, &n, sizeof(n)) != sizeof(n)) {
        printf("ERROR reading header E: %d %s\n", errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (n > (1000 * 1000 * 1000)) {
        printf("ERROR number of records %'ld is too big\n", n);
        exit(EXIT_FAILURE);
    }

    if (n < _n) {
        printf("ERROR file only contents %'ld values, less than required (%'ld)\n", n, _n);
        exit(EXIT_FAILURE);
    }

    int64_t block = 10240 * 8;
    int64_t bytes = _n*32;
    int64_t bytesRead;

    while (bytes > 0) {
        bytesRead = read(fd, scalars, bytes > block ? block : bytes);
        scalars += bytesRead;
        bytes -= bytesRead;
    }

    lseek64(fd, sizeof(n) + n * 32, SEEK_SET);

    uint8_t *data = (uint8_t *)bases;
    bytes = sizeof(bases[0]) * _n;
    while (bytes > 0) {
        bytesRead = read(fd, data, bytes > block ? block : bytes);
        data += bytesRead;
        bytes -= bytesRead;
    }

    close(fd);
}


int main(int argc, char **argv) 
{
    int nscalars = 32;
    int n = atoi(argv[1]);

    setlocale(LC_ALL, "en_US.utf-8");

    uint8_t *scalars = new uint8_t[n*32];
    G1PointAffine *bases = new G1PointAffine[n];

    loadDataFile("multiexp_test_data_128000000.dat", scalars, bases, n);

    clock_t start, end;
    int64_t startT, endT;
    double cpu_time_used;
    std::string strResult;

    G1Point p1;
/*
    printf("Starting multiexp. (original) \n");
    G1.resetCounters();
    start = clock();
    startT = getRealTimeClockUs();
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, nscalars, n);
    end = clock();
    endT = getRealTimeClockUs();

    strResult = G1.toString(p1); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time real: %.4lf\n", (double)(endT - startT)/1000000);
    printf("Time used: %.4lf\n", cpu_time_used);
    printf("Avg time per exp: %.4lf us\n", (cpu_time_used*1000000)/n);
    printf("Exps per second: %.4lf\n", (n / cpu_time_used));


    G1Point p2;

    printf("cnt ==== Starting multiexp. (mix)  \n");
    G1.resetCounters();
    start = clock();
    startT = getRealTimeClockUs();
    G1.multiMulByScalarMix(p2, bases, (uint8_t *)scalars, nscalars, n);
    end = clock();
    endT = getRealTimeClockUs();
    strResult = G1.toString(p2); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());
    int64_t total = G1.cntAddMixed + G1.cntAdd + G1.cntAddAffine;
    printf("P1 cntAddTotal:%'ld\n", total);

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time real: %.4lf\n", (double)(endT - startT)/1000000);
    printf("Time used: %.4lf\n", cpu_time_used);
    printf("Avg time per exp: %.4lf us\n", (cpu_time_used*1000000)/n);
    printf("Exps per second: %.4lf\n", (n / cpu_time_used));
*/
    G1Point p3;

    printf("cnt ==== Starting multiexp. \n");
    G1.resetCounters();
    start = clock();
    startT = getRealTimeClockUs();
    G1.multiMulByScalarBa(p3, bases, (uint8_t *)scalars, nscalars, n);
    end = clock();
    endT = getRealTimeClockUs();
    strResult = G1.toString(p3); 
    printf("P1 (%s):%s\n", (strResult == EXPECTED_RESULT ? "OK":"********FAIL********"), strResult.c_str());
    int64_t total2 = G1.cntAddMixed + G1.cntAdd + G1.cntAddAffine;
    printf("P1 cntAddTotal:%'ld\n", total2);

    G1.printCounters();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time real: %.4lf\n", (double)(endT - startT)/1000000);
    printf("Time used: %.4lf\n", cpu_time_used);
    printf("Avg time per exp: %.4lf us\n", (cpu_time_used*1000000)/n);
    printf("Exps per second: %.4lf\n", (n / cpu_time_used));
    
    // printf("P1 Better %.02f%%\n", ((double)(total-total2)*100.0)/total);
}