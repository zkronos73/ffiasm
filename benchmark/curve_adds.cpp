#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <locale.h>
#include "alt_bn128.hpp"
#include <time.h>

#include <string>

using namespace AltBn128;

__uint128_t g_lehmer64_state = 0xAAAAAAAAAAAAAAAALL;
int64_t globalIndex = 0;

// Fast random generator
// https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/

uint64_t lehmer64() {
  g_lehmer64_state *= 0xda942042e4dd58b5LL;
  return g_lehmer64_state >> 64;
}

class CurveAdds {
    public:
        CurveAdds ( void );
        ~CurveAdds ( void );
        void help ( void );
        void run ( int argc, char **argv );

    protected:
        int64_t n;
        int64_t nBases;
        int64_t offsetBases;
        std::string outputFile;
        bool flgMultiSums;
        bool flgParallelMultiSums;
        bool flgSums;
        bool flgBases;
        bool flgVerifySums;
        bool flgBasesLoad;
        bool flgBasesSave;
        bool flgInverses;
        bool flgSetN;
        int times;
        std::string basesFile;

        double multiSumsSeconds;
        double sumsSeconds;
        G1PointAffine *affineBases;
        G1Point *bases;
        G1Point *sums;
        G1PointAffine *multiSums;

    public:
        void generateBases ( void );
        void freeBases ( void );
        void freeSums ( void );
        int64_t getRealTimeClockUs ( void );
        void freeMultiSums ( void );
        void calculateMultiSums ( void );
        void calculateParallelMultiSums ( void );
        void calculateSums ( void );
        void calculateBaseOperations ( void );
        void calculateInverses ( void );
        void verifySums ( void );
        void parseArguments ( int argc, char **argv );
        void loadBases ( void );
        void saveBases ( void );
        void largeIO ( const std::string &label, const std::string &filename, bool writeMode, int fd,  void *data, int64_t bytes );
        void largeRead ( const std::string &label, const std::string &filename, int fd, void *data, int64_t bytes ) { largeIO(label, filename, false, fd, data, bytes); };
        void largeWrite ( const std::string &label, const std::string &filename, int fd, void *data, int64_t bytes ) { largeIO(label, filename, true, fd, data, bytes); };
};

CurveAdds::CurveAdds ( void )
{
    n = 100000;
    times = 1;
    flgSums = false;
    flgMultiSums = false;
    flgParallelMultiSums = false;
    flgVerifySums = false;
    flgBases = false;
    flgInverses = false;
    affineBases = NULL;
    bases = NULL;
    sums = NULL;
    multiSums = NULL;
    multiSumsSeconds = 0;
    flgBasesLoad = false;
    flgBasesSave = false;
    flgSetN = false;
    sumsSeconds = 0;
    offsetBases = 0;
}

CurveAdds::~CurveAdds ( void )
{
    freeBases();
    freeSums();
    freeMultiSums();
}

void CurveAdds::help ( void )
{

}

void CurveAdds::parseArguments ( int argc, char **argv )
{
    int opt;

    while ((opt = getopt(argc, argv, "n:t:o:mbavihsp:l:")) != -1) {
        switch (opt) {
            case 'n':
                n = atoi(optarg);
                flgSetN = true;
                break;

            case 't':
                times = atoi(optarg);
                break;

            case 'o':
                offsetBases = atoi(optarg);
                break;

            case 'm':
                flgMultiSums = true;
                break;

            case 'p':
                flgParallelMultiSums = true;
                break;

            case 'i':
                flgInverses = true;
                break;

            case 'b':
                flgBases = true;
                break;

            case 'a':
                flgSums = true;
                break;

            case 'v':
                flgVerifySums = true;
                break;

            case 'h':
                help();
                exit(EXIT_SUCCESS);
            
            case 's':
                flgBasesSave = true;
                basesFile = optarg;
                break;

            case 'l':
                flgBasesLoad = true;
                basesFile = optarg;
                break;
                
            default: /* '?' */
                help();
                exit(EXIT_FAILURE);
        }
    }
/*
    if (optind >= argc) {
        fprintf(stderr, "Expected argument after options\n");
        exit(EXIT_FAILURE);
    }*/

    nBases = 2 * n;
}

void CurveAdds::freeBases ( void )
{
    if (bases) {
        delete bases;
        bases = NULL;
    }
    if (affineBases) {
        delete affineBases;
        affineBases = NULL;
    }
}

void CurveAdds::freeSums ( void )
{
    if (sums) {
        delete sums;
        sums = NULL;
    }
}

void CurveAdds::freeMultiSums ( void )
{
    if (multiSums) {
        delete multiSums;
        multiSums = NULL;
    }
}

void CurveAdds::generateBases ( void )
{
    freeBases();

    if (flgBasesLoad) {
        loadBases();
        return;
    }

    affineBases = new G1PointAffine[nBases];
    bases = new G1Point[nBases];

    printf("Generation of %'ld affine bases (%p, %p)\n", nBases, affineBases, bases);

    G1.copy(affineBases[0], G1.one());
    G1.copy(affineBases[1], G1.one());

    for (int64_t i=2; i < nBases; i++) {
        if ( i%100000 == 0) printf(" · %6.02f%% affine bases generated\n", (double)i * 100 / nBases);
        G1.add(affineBases[i], affineBases[i-1], affineBases[i-2]);
    }

    printf(" · 100.00%% affine bases generated\n\n");

    printf("Generation of %'ld bases\n", nBases);

    G1.copy(bases[0], G1.one());
    G1.copy(bases[1], G1.one());

    for (int64_t i=2; i < nBases ; i++) {
        if ( i%100000 == 0) printf(" · %6.02f%% bases generated\n", (double)i * 100 / nBases);
        G1.add(bases[i], bases[i-1], bases[i-2]);
    }
    printf(" · 100.00%% bases generated\n\n");
    if (flgBasesSave) {
        saveBases();
    }

}

void CurveAdds::largeIO (const std::string &label, const std::string &filename, bool writeMode, int fd,  void *data, int64_t bytes) 
{
    int64_t bytesIO;
    uint8_t *_data = (uint8_t *)data;
    uint64_t steps = 0;

    printf("prepare to %s %'ld bytes of %s %s file %s\n", writeMode ? "write":"read", bytes, label.c_str(), writeMode ? "to":"from", filename.c_str());
    while (bytes > 0) {        
        bytesIO = writeMode ? write(fd, _data, bytes) : read(fd, _data, bytes);
        if (bytesIO < 0) {
            printf("ERROR %d (%s) %s %s %'ld bytes %s\n", errno, strerror(errno), writeMode ? "saving":"loading", label.c_str(), bytes, filename.c_str());
            close(fd);
            exit(EXIT_FAILURE);
        }
        if (bytesIO == 0) {
            printf("ERROR %d (%s) %s (no %s bytes, but %'ld bytes remainding) %s\n", errno, strerror(errno), writeMode ? "saving":"loading", label.c_str(), bytes, filename.c_str());
            close(fd);
            exit(EXIT_FAILURE);            
        }
        bytes -= bytesIO;
        _data += bytesIO;
        ++steps;
    }
/*
    if (steps > 1) {
        printf("WARNING: %s %s in %'ld operations\n", writeMode ? "write":"read", label.c_str(), steps);
    }*/
}

void CurveAdds::saveBases()
{
    int fd = creat(basesFile.c_str(), 0666);
    if (fd < 0) {
        printf("ERROR %d (%s) saving bases file %s\n", errno, strerror(errno), basesFile.c_str());
        exit(EXIT_FAILURE);
    }
    if (write(fd, &nBases, sizeof(nBases)) < 0) {
        printf("ERROR %d (%s) saving bases file %s\n", errno, strerror(errno), basesFile.c_str());
        close(fd);
        exit(EXIT_FAILURE);
    }

    largeWrite("affineBases", basesFile, fd, affineBases, sizeof(affineBases[0]) * nBases);
    largeWrite("bases", basesFile, fd, bases, sizeof(bases[0]) * nBases);

    close(fd);
    printf("Saved %'ld bases (affine + non-affine) to %s\n", nBases, basesFile.c_str());
}


void CurveAdds::loadBases()
{
    int fd = open(basesFile.c_str(), O_RDONLY);
    if (fd < 0) {
        printf("ERROR %d (%s) loading bases file %s\n", errno, strerror(errno), basesFile.c_str());
        exit(EXIT_FAILURE);
    }
    
    int64_t count;
    if (read(fd, &count, sizeof(count)) != sizeof(count)) {
        printf("ERROR %d (%s) loading bases file %s\n", errno, strerror(errno), basesFile.c_str());
        exit(EXIT_FAILURE);
    }

    if (!flgSetN || count < nBases) {
        nBases = count;
        n = nBases >> 1;
    }
    affineBases = new G1PointAffine[nBases];    
    if (!affineBases) {
        printf("ERROR allocating memory for %'ld affineBases\n", nBases);
        exit(EXIT_FAILURE);
    }

    bases = new G1Point[nBases];
    if (!bases) {
        printf("ERROR allocating memory for %'ld bases\n", nBases);
        exit(EXIT_FAILURE);
    }

    printf("Loading %'ld bases from %s\n", nBases, basesFile.c_str());

    if (offsetBases > 0) {
        if (lseek64(fd, sizeof(nBases) + sizeof(affineBases[0]) * offsetBases, SEEK_SET) < 0) {
            printf("ERROR %d (%s) loading (lseek %'ld) bases file %s\n", errno, strerror(errno), sizeof(nBases) + sizeof(affineBases[0]) * offsetBases, basesFile.c_str());
            exit(EXIT_FAILURE);
        }
    }
    largeRead("affineBases", basesFile, fd, affineBases, sizeof(affineBases[0]) * nBases);

    if (lseek64(fd, sizeof(nBases) + sizeof(affineBases[0]) * count + sizeof(bases[0]) * offsetBases, SEEK_SET) < 0) {
        printf("ERROR %d (%s) loading (lseek %'ld) bases file %s\n", errno, strerror(errno), sizeof(nBases) + sizeof(affineBases[0]) * offsetBases, basesFile.c_str());
        close(fd);
        exit(EXIT_FAILURE);
    }

    largeRead("bases", basesFile, fd, bases, sizeof(bases[0]) * nBases);
    close(fd);

    printf("Loaded %'ld bases (affine + non-affine) from %s\n", nBases, basesFile.c_str());
}

/*
void CurveAdds::calculateInverses ( void )
{
    clock_t start, end;

    F1Element *values = new F1Element[n];
    F1Element *results = new F1Element[n];
    F1Element *results2 = new F1Element[n];
    if (!values || !results) {
        printf("ERROR allocating memory for %'ld F1Elements (values %p, results %p)\n", n, values, results);
        exit(EXIT_FAILURE);
    }
    printf("\nPreparing %'ld values ....\n", n);
    for (int i = 0; i < n; ++i) {
        values[i] = affineBases[i].x;
    }
    G1.F.batchInverse(results, values, n);
    printf("\nStarting %'ld batchInverse I ....\n", n);
    start = clock();
    G1.F.batchInverse(results, values, n);
    end = clock();
    multiSumsSeconds += ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Seconds used on %'ld batchInverse I: %.2lf\n", n, multiSumsSeconds);

    G1.F.batchInverse_2(results2, values, n);
    printf("\nStarting %'ld batchInverse II ....\n", n);
    start = clock();
    G1.F.batchInverse_2(results2, values, n);
    end = clock();
    multiSumsSeconds += ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Seconds used on %'ld batchInverse II: %.2lf\n", n, multiSumsSeconds);
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        if (!G1.F.eq(results[i], results2[i])) {
            ++cnt;
        }
    }
    printf("%d differents values\n", cnt);
}
*/

void CurveAdds::calculateInverses ( void )
{
    clock_t start, end;

    typedef struct {
        F1Element value;
        F1Element result;
    } TwoF1Elements;

    int loops = 1; // 500000;
    int loopSize = n / loops;
    TwoF1Elements *elements = new TwoF1Elements[n];
    F1Element *results = new F1Element[n];
    F1Element *values = new F1Element[n];
    // F1Element *results2 = new F1Element[n];
    if (!values || !results) {
        printf("ERROR allocating memory for %'ld F1Elements (values %p, results %p)\n", n, values, results);
        exit(EXIT_FAILURE);
    }
    printf("\nPreparing %'ld values ....\n", n);
    for (int i = 0; i < n; ++i) {
        if (G1.F.isZero(affineBases[i].x)) {
            printf("ERROR on %d\n",i );
            exit(EXIT_FAILURE);
        }
        G1.F.copy(values[i], affineBases[i].x);
        G1.F.copy(elements[i].value, affineBases[i].x);
    }
    G1.F.batchInverse(results, values, n);
    printf("\nStarting %'ld batchInverse (loops:%d block:%d) I ....\n", n, loops, loopSize);
    start = clock();
    for (int loop = 0; loop < loops; ++loop) {
        G1.F.batchInverse(results + loop * loopSize, values + loop * loopSize, loopSize);
    }
    end = clock();
    printf("Seconds used on %'ld batchInverse I: %.2lf\n", n, ((double) (end - start)) / CLOCKS_PER_SEC);

    // G1.F.batchInverse_3(results2, sizeof(results2[0]), values, sizeof(values[0]), k);
    G1.F.batchInverse_3(&elements[0].result, sizeof(elements[0]), &elements[0].value, sizeof(elements[0]), n);
    printf("\nStarting %'ld batchInverse II ....\n", n);
    start = clock();
    for (int loop = 0; loop < loops; ++loop) {
        // G1.F.batchInverse_3(&elements[loopSize * loop].result, 0 /*, sizeof(elements[0])*/, &elements[loopSize * loop].value, 0 /*sizeof(elements[0])*/, loopSize);
        G1.F.batchInverse_3(&elements[0].result, 0 /*, sizeof(elements[0])*/, &elements[0].value, 0 /*sizeof(elements[0])*/, loopSize);
    }
    end = clock();
    printf("Seconds used on %'ld batchInverse II: %.2lf\n", n, ((double) (end - start)) / CLOCKS_PER_SEC);
    int cntOk = 0;
    int cntFail = 0;
    for (int i = 0; i < n; ++i) {
//        printf("#V %s\n#1 %s\n#2 %s\n", G1.F.toString(values[i]).c_str(), G1.F.toString(results[i]).c_str(), G1.F.toString(results2[i]).c_str());
        if (!G1.F.eq(results[i], elements[i].result)) {
            ++cntFail;
        }
        else {
            ++cntOk;
        }
    }
    printf("%d differents values of %d\n", cntFail, cntOk);
}

void CurveAdds::calculateMultiSums ( void )
{
    clock_t start, end;

    freeMultiSums();
    multiSums = new G1PointAffine[n];
    if (!multiSums) {
        printf("ERROR allocating memory for %'ld G1PointAffine\n", n);
        exit(EXIT_FAILURE);
    }
    printf("\nStarting %'ld multi sums ....\n", n);
    #ifdef FFIASM_FR_COUNTERS
    F1Element k;
    RawFq::Stats stats = G1.F.stats;
    #endif
    start = clock();
    G1.multiAdd(multiSums, affineBases, affineBases + n, n);
    end = clock();
    #ifdef FFIASM_FR_COUNTERS
    RawFq::Stats stats2 = G1.F.stats;
    printf("cntAdd:%'ld\n",  stats2.cntAdd-stats.cntAdd);
    printf("cntSub:%'ld\n",  stats2.cntSub-stats.cntSub);
    printf("cntMMul:%'ld\n", stats2.cntMMul-stats.cntMMul);
    printf("cntSquare:%'ld\n", stats2.cntSquare-stats.cntSquare);
    printf("cntMul1:%'ld\n", stats2.cntMul1-stats.cntMul1);
    #endif
    multiSumsSeconds += ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Seconds used on %'ld multi sums: %.2lf\n", n, multiSumsSeconds);
}

void CurveAdds::calculateParallelMultiSums ( void )
{
    clock_t start, end;
    int64_t startUs, endUs;


    freeMultiSums();
    multiSums = new G1PointAffine[n];
    if (!multiSums) {
        printf("ERROR allocating memory for %'ld G1PointAffine\n", n);
        exit(EXIT_FAILURE);
    }
    printf("\nStarting %'ld multi sums ....\n", n);
    #ifdef FFIASM_FR_COUNTERS
    printf("\x1b[31mWARNING active FFIASM_FR_COUNTERS !!!\x1b[0m\n");
    #endif
    #ifdef COUNT_OPS
    printf("\x1b[31mWARNING active COUNT_OPS !!!\x1b[0m\n");
    #endif
    start = clock();
    startUs = getRealTimeClockUs();
    G1.multiAdd(multiSums, affineBases, affineBases + n, n);
    end = clock();
    endUs = getRealTimeClockUs();
    printf("Seconds used on %'ld multi sums [real:%.4lf cpu:%.4lf]\n", n, ((double) (endUs - startUs)) / 1000000.0, ((double) (end - start)) / CLOCKS_PER_SEC);
}

void CurveAdds::calculateSums ( void )
{
    clock_t start, end;

    freeSums();
    sums = new G1Point[n];
    if (!sums) {
        printf("ERROR allocating memory for %'ld G1Point\n", n);
        exit(EXIT_FAILURE);
    }

    printf("\nStarting %'ld sums ....\n", n);
    #ifdef FFIASM_FR_COUNTERS
    F1Element k;
    RawFq::Stats stats = G1.F.stats;
    #endif

    auto startUs = getRealTimeClockUs();
    start = clock();
    for (int64_t i=0; i < n; ++i) {
        // globalIndex = i;
        G1.add(sums[i], affineBases[i], bases[n + i]);
    }
    end = clock();
    auto endUs = getRealTimeClockUs();    
    #ifdef FFIASM_FR_COUNTERS
    RawFq::Stats stats2 = G1.F.stats;
    printf("cntAdd:%'ld\n",  stats2.cntAdd-stats.cntAdd);
    printf("cntSub:%'ld\n",  stats2.cntSub-stats.cntSub);
    printf("cntMMul:%'ld\n", stats2.cntMMul-stats.cntMMul);
    printf("cntSquare:%'ld\n", stats2.cntSquare-stats.cntSquare);
    printf("cntMul1:%'ld\n", stats2.cntMul1-stats.cntMul1);
    #endif
    sumsSeconds += ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Seconds used on %'ld sums: %.4lf (%.4lf)\n", n, sumsSeconds, (double)(endUs - startUs)/1000000.0);
}

inline int64_t CurveAdds::getRealTimeClockUs ( void )
{
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    return tm.tv_sec * 1000000 + tm.tv_nsec / 1000;
}

void CurveAdds::calculateBaseOperations ( void )
{
    clock_t start[6], end[6];

    freeSums();
    auto sums = new G1PointAffine[n];
    if (!sums) {
        printf("ERROR allocating memory for %'ld G1Point\n", n);
        exit(EXIT_FAILURE);
    }

    #ifdef FFIASM_FR_COUNTERS
    F1Element k;
    RawFq::Stats stats = G1.F.stats;
    #endif
    #define R_IND (8192 + (i & 0x001))
    #define R_OP1 (8192 + (i & 0x001))
    #define R_OP2 (8912 + ((i+1) & 0x001))
    #define PREFETCH_N 1
    auto ndiv2 = n / 2;
    for (int loop = 0; loop < 10; ++loop) {
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.add(sums[R_IND].x, affineBases[R_OP1].x, affineBases[R_OP2].x);
            G1.F.add(sums[R_IND].y, affineBases[R_OP1].y, affineBases[R_OP2].y);
        }
        int timerIndex = 0;
        printf("\nStarting %'ld sums ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.add(sums[R_IND].x, affineBases[R_OP1].x, affineBases[R_OP2].x);
            G1.F.add(sums[R_IND].y, affineBases[R_OP1].y, affineBases[R_OP2].y);
        }
        end[timerIndex++] = clock();
        printf("\nStarting %'ld sub ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.sub(sums[R_IND].x, affineBases[R_OP1].x, affineBases[R_OP2].x);
            G1.F.sub(sums[R_IND].y, affineBases[R_OP1].y, affineBases[R_OP2].y);
        }
        end[timerIndex++] = clock();
        printf("\nStarting %'ld mul ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.mul(sums[R_IND].x, affineBases[R_OP1].x, affineBases[R_OP2].x);
            G1.F.mul(sums[R_IND].y, affineBases[R_OP1].y, affineBases[R_OP2].y);
        }
        end[timerIndex++] = clock();
        printf("\nStarting %'ld squares ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.square(sums[R_IND].x, affineBases[R_OP1].x);
            G1.F.square(sums[R_IND].y, affineBases[R_OP1].y);
        }
        end[timerIndex++] = clock();
        printf("\nStarting %'ld inv ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            G1.F.inv(sums[R_IND].x, affineBases[R_OP1].x);
            G1.F.inv(sums[R_IND].y, affineBases[R_OP1].y);
        }
        end[timerIndex++] = clock();
        printf("\nStarting %'ld moves ....\n", n);
        start[timerIndex] = clock();
        for (int64_t i=0; i < ndiv2; ++i) {
            __builtin_prefetch(affineBases + R_OP1 + PREFETCH_N, 0);
            __builtin_prefetch(affineBases + R_OP2 + PREFETCH_N, 0);
            __builtin_prefetch(sums + R_IND + PREFETCH_N, 1);            
            sums[R_IND].y = affineBases[R_OP1].x;
            sums[R_IND].x = affineBases[R_OP2].x;
            sums[R_IND].x = affineBases[R_OP1].y;
            sums[R_IND].y = affineBases[R_OP2].y;
        }
        end[timerIndex++] = clock();
        #ifdef FFIASM_FR_COUNTERS
        RawFq::Stats stats2 = G1.F.stats;
        printf("cntAdd:%'ld\n",  stats2.cntAdd-stats.cntAdd);
        printf("cntSub:%'ld\n",  stats2.cntSub-stats.cntSub);
        printf("cntMMul:%'ld\n", stats2.cntMMul-stats.cntMMul);
        printf("cntSquare:%'ld\n", stats2.cntSquare-stats.cntSquare);
        printf("cntMul1:%'ld\n", stats2.cntMul1-stats.cntMul1);
        #endif
        printf("Seconds used on %'ld x add: %.4lf\n", n, ((double) (end[0] - start[0])) / CLOCKS_PER_SEC);
        printf("Seconds used on %'ld x sub: %.4lf\n", n, ((double) (end[1] - start[1])) / CLOCKS_PER_SEC);
        printf("Seconds used on %'ld x mmul: %.4lf\n", n, ((double) (end[2] - start[2])) / CLOCKS_PER_SEC);
        printf("Seconds used on %'ld x square: %.4lf\n", n, ((double) (end[3] - start[3])) / CLOCKS_PER_SEC);
        printf("Seconds used on %'ld x inv: %.4lf\n", n, ((double) (end[4] - start[4])) / CLOCKS_PER_SEC);
        printf("Seconds used on %'ld x move: %.4lf\n", n, ((double) (end[5] - start[5])) / CLOCKS_PER_SEC);
    }
}

void CurveAdds::verifySums ( void )
{
    if (!sums) return;
    if (!multiSums) return;

    for (auto i=0; i<n; i++) {
        if (!G1.eq(multiSums[i], sums[i])) {
            printf("Index %d different:\n", i);
            exit(EXIT_FAILURE);
        }
    }   
    printf("verify Sums OK!\n");
}

void CurveAdds::run( int argc, char **argv )
{
    parseArguments(argc, argv);
    generateBases();

    multiSumsSeconds = 0;
    sumsSeconds = 0;

    for (int i = 0; i < times; ++i) {
        if (flgSums) {
            calculateSums();
        }
        if (flgMultiSums) {
            calculateMultiSums();
        }
        if (flgParallelMultiSums) {
            calculateParallelMultiSums();
        }
        if (flgBases) {
            calculateBaseOperations();
        }
        if (flgInverses) {
            calculateInverses();
        }
    }

    if (flgSums && flgMultiSums) {
        printf("\nratio: %.04f (multi/sum)\n", multiSumsSeconds/sumsSeconds);
        if (multiSumsSeconds < sumsSeconds) {
            printf("multisum is %.02f%% more faster than sums\n", (sumsSeconds - multiSumsSeconds) / sumsSeconds);
        }
        if (multiSumsSeconds > sumsSeconds) {
            printf("WARNING: sum is %.02f%% more faster than multisum\n", (multiSumsSeconds - sumsSeconds) / sumsSeconds);
        }
    }
    if (flgVerifySums) {
        verifySums();
    }
}

void segfault_sigaction(int signal, siginfo_t *si, void *arg)
{
    printf("Caught segfault at address %p\n", si->si_addr);
    printf("value of globalIndex: %'ld\n", globalIndex);
    exit(0);
}

int main(int argc, char **argv) 
{
    CurveAdds ca;
    struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = segfault_sigaction;
    sa.sa_flags   = SA_SIGINFO;

    sigaction(SIGSEGV, &sa, NULL);

    setlocale(LC_ALL, "en_US.utf-8");
    ca.run(argc, argv);
}