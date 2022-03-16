#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
// #include <malloc.h>
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

class MultiExpG1
{
    public:

        MultiExpG1 ( void );
        ~MultiExpG1 ( void );
        void run (int argc, char **argv);

    protected:

        int nscalars;
        int loop;
        int64_t n;
        int times;
        const int modes = 2;
        bool flgSaveDataFile;
        bool flgLoadDataFile;
        bool flgOriginalMode;
        bool flgMultiMode;
        std::string dataFilename;
        double cpuTime [2];
        double realTime [2];
        uint8_t *scalars;
        G1PointAffine *bases;

        inline int64_t getRealTimeClockUs ( void );
        void createDataFile ( const std::string &filename, int64_t _n );
        void loadDataFile ( const std::string &filename, uint8_t *scalars, G1PointAffine *bases, int64_t _n );
        void allocateData ( void );
        void freeData ( void );
        void executeBenchmark ( int mode = 0 );
        void executeBenchmarks ( void );    
        const std::string &getExpectedResult ( void );
        void parseArguments ( int argc, char *argv[] );
        void help ( const std::string &prgname );
        void summary ( void );
        int64_t showMemoryUsage ( int64_t previous = 0 );
};

inline int64_t MultiExpG1::getRealTimeClockUs ( void )
{
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    return tm.tv_sec * 1000000 + tm.tv_nsec / 1000;
}

void MultiExpG1::createDataFile ( const std::string &filename, int64_t _n )
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

void MultiExpG1::loadDataFile ( const std::string &filename, uint8_t *scalars, G1PointAffine *bases, int64_t _n )
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

int64_t MultiExpG1::showMemoryUsage ( int64_t previous )
{
/*    struct mallinfo2 minfo;
    minfo = mallinfo2();
    int64_t used = minfo.uordblks;
    int64_t factor = 1;
    std::string units = "bytes";

    if (used > 1000000000) {
        factor = 1000000000;
        units = "GB";
    }
    else if (used > 1000000) {
        factor = 1000000;
        units = "MB";
    }
    else if (used > 1000) {
        factor = 1000;
        units = "KB";
    }

    printf("[ mem: %0.03f %s", (double)used/(double)factor, units.c_str());
    if (previous) {
        uint64_t diff = used - previous;
        if (!diff) {
            printf(" = \x1b[32mOK - no memory lost\x1b[0m");
        }
        else {
            printf(" \x1b[31m%s%0.03f %s WARNING found memory difference\x1b[0m", (diff > 0.0 ? "+":""), (double)diff/(double)factor, units.c_str());
        }
    }
    printf(" ]\n");
    return used;*/
    return 0;
}

void MultiExpG1::executeBenchmark ( int mode )
{
    G1Point p1;
 
    clock_t start, end;
    int64_t startT, endT;
    double cpu_time_used;

    double total1, total3;
    double cpu1, cpu3;

    uint64_t startMem = showMemoryUsage();

    #ifdef COUNT_OPS
    G1.resetCounters();
    #endif
    start = clock();
    startT = getRealTimeClockUs();
    switch (mode) {
        case 0:
            G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, nscalars, n);
            break;

        case 1:
            G1.multiMulByScalarBa(p1, bases, (uint8_t *)scalars, nscalars, n);
            break;
    }
    end = clock();
    endT = getRealTimeClockUs();
    #ifdef COUNT_OPS
    printf("AddM: %'d | Add: %'d | AddA: %'d | AddT: %'d | Dbl: %'d | DblM: %'d | DblT: %'d | Eq: %'d | EqM: %'d | ToA: %'d\n",
        G1.cntAddMixed, G1.cntAdd, G1.cntAddAffine, G1.cntAddMixed + G1.cntAdd + G1.cntAddAffine,
        G1.cntDbl, G1.cntDblMixed, G1.cntDbl + G1.cntDblMixed, G1.cntEq, G1.cntEqMixed, G1.cntToAffine);
    #endif

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("real: %.4lf | cpu: %.4lf | avg/exp: %.4lf | exp/seg: %.4lf\n", (double)(endT - startT)/1000000 , cpu_time_used,
        (cpu_time_used*1000000)/n, (n / cpu_time_used));

    std::string expectedResult = getExpectedResult();
    std::string result = G1.toString(p1); 
    printf("P1 ");
    bool match = true;
    if (!expectedResult.empty()) {
        match = (expectedResult == result);
        printf("%s ",  match ? "[\x1b[32mOK\x1b[0m]  ":"[\x1b[31mFAIL\x1b[0m]");
    }
    printf("%s", result.c_str());
    if (!match) {
        printf("\nEXPECTED  %s", expectedResult.c_str());
    }
    printf("\n");

    showMemoryUsage(loop ? startMem : 0);

    realTime[mode] +=  (double)(endT - startT)/1000000;
    cpuTime[mode] += cpu_time_used;
}

MultiExpG1::MultiExpG1 ( void )
{
    nscalars = 32;
    n = 1000000;
    times = 1;
    scalars = NULL;
    bases = NULL;
    flgSaveDataFile = false;
    flgLoadDataFile = true;
    flgOriginalMode = false;
    flgMultiMode = false;
    dataFilename = "multiexp_test_data_128000000.dat";

    setlocale(LC_ALL, "en_US.utf-8");

    for (int index = 0; index < modes; ++index) {
        realTime[index] = 0.0;
        cpuTime[index] = 0.0;
    }
}

MultiExpG1::~MultiExpG1 ( void )
{
    freeData();
}

void MultiExpG1::help ( const std::string &prgname )
{
    std::string cmd = prgname.substr(prgname.find_last_of("/\\") + 1);

    printf("usage: %s [-h] [-1] [-2] [-a] [-s <filename>] [-l <filename>] [-n <#points>] [-t <#loops>]\n", cmd.c_str());
    printf(" where:\n");
    printf("  -h show this help.\n");
    printf("  -1 benchmark using one-by-one add.\n");
    printf("  -2 benchmark using batch adds.\n");
    printf("  -a make all benchmarks.\n");
    printf("  -s <filename> generate a data file <filename> with #points.\n");
    printf("  -l <filename> load data from file <filename>.\n");
    printf("  -n <#points> benchmark with #points, could use M or K suffix.\n");
    printf("  -t <#loops> number of repetitions.\n\n");
}

void MultiExpG1::parseArguments( int argc, char **argv ) 
{
    int opt;

    while ((opt = getopt(argc, argv, "n:t:a12hs:l:")) != -1) {
        switch (opt) {
            case 'n':
            {
                int len = strlen(optarg);
                n = 1;
                if (len > 0) {
                    if (optarg[len-1] == 'M') n = 1000000;
                    else if (optarg[len-1] == 'K') n = 1000;
                }
                n = n * atoi(optarg);
                break;
            }
            case 't':
                times = atoi(optarg);
                break;

            case 'a':
                flgMultiMode = true;
                flgOriginalMode = true;
                break;

            case '2':
                flgMultiMode = true;
                break;

            case '1':
                flgOriginalMode = true;
                break;

            case 'h':
                help(argv[0]);
                exit(EXIT_SUCCESS);

            case 's':
                flgSaveDataFile = true;
                dataFilename = optarg;
                break;

            case 'l':
                flgLoadDataFile = true;
                dataFilename = optarg;
                break;
                
            default: /* '?' */
                help(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
}    

void MultiExpG1::allocateData ( void )
{
    printf("allocating %'ld scalars (%ld MB) ...\n", n, (n * 32)/1000000);
    scalars = new uint8_t[n*32];

    printf("allocating %'ld bases (%ld MB) ...\n", n, (n * sizeof(G1PointAffine))/1000000);
    bases = new G1PointAffine[n];
}

void MultiExpG1::freeData ( void )
{
    if (scalars) {
        free(scalars);
        scalars = NULL;
    }
    if (bases) {
        free(bases);
        bases = NULL;
    }
}

void MultiExpG1::run (int argc, char **argv) 
{
    parseArguments(argc, argv);
    
    if (flgSaveDataFile) {
        createDataFile(dataFilename, n);
        exit(EXIT_SUCCESS);
    }

    allocateData();

    if (!flgOriginalMode && !flgMultiMode) {
        help(argv[0]);
        exit(EXIT_SUCCESS);
    }

    if (flgLoadDataFile) {
        loadDataFile(dataFilename, scalars, bases, n);
    }

    executeBenchmarks();

    summary();
}

void MultiExpG1::summary ( void )
{
    printf("\n\n=== SUMMARY ===\n\n");
    printf("n: %'ld x %'d times\n", n, times);
    if (flgOriginalMode) {
        printf("original avg(real): %.4lf  avg(cpu): %.4lf seconds\n", realTime[0]/(double) times, cpuTime[0]/(double)times );
    }
    if (flgMultiMode) {
        printf("multiadd avg(real): %.4lf  avg(cpu): %.4lf seconds\n", realTime[1]/(double) times, cpuTime[1]/(double)times );
    }
    if (flgOriginalMode && flgMultiMode) {
        printf("multi vs original real: %0.2lf%%  cpu: %0.2lf%%\n", realTime[1]*100.0/realTime[0], cpuTime[1]*100.0/cpuTime[0]);
    }
}

void MultiExpG1::executeBenchmarks ( void )
{
    for (loop = 0; loop < times; ++loop) {
        if (flgOriginalMode) {
            printf("\n====> BENCHMARK (normal) %d/%d <====\n\n", loop + 1, times);
            executeBenchmark(0);
        }

        if (flgMultiMode) {
            printf("\n====> BENCHMARK (multi-add) %d/%d <====\n\n", loop + 1, times);
            executeBenchmark(1);
        }
    }
}

const std::string &MultiExpG1::getExpectedResult ( void )
{
    static const struct {
        int64_t n;
        std::string result;
    } results[] = { 
        { 128000000, "(7619599290849605139674112334648777838192883052106003259536319974986231832328,9581711967996342227363571221298777964748202435299455322059359416458683753394)"},
        { 100000000,  "(4657314544325888627957202934312756489824753756071102142771983332701648414869,3316764596681176473417459795061529664032079299847944460908660216822682297676)"},
        { 64000000, "(14027679857912970260632559145576090903498424321393515460717381462887634325449,9867510856157004017513284373910286041492679570957103391987912936413016686938)"},
        { 32000000, "(3725627854558833746836478303362535384290726728799886813929673107562122721569,17171034889966285705425666698235454513929560386782128996039157896572606207470)"},
        { 16000000, "(6543188983016917526542824388468771425679360998945436057822982100797713736596,17753380454259490450850315541625409008054761226270337663717967071441276427171)"},
        { 10000000, "(6086097182414925163347122340328541596364360534174103734671785045522740508610,14647465025473290936996772034695920750119607499710914676572169653112563974073)"},
        { 8000000, "(1541909324057016050073316337252903946297640000394358240240890919105074030131,5704284559889083593836388165452987109010626349939979695851698579856715926171)"},
        { 4000000, "(16311154981726638560478224772782498352996355389530878784848285224044621066290,12856304264105132790446629774993624637149500326272925810477620682150834845440)"},
        { 2000000, "(18673661583012128740346417472298503774180827300104501275450680830677636626882,13148789420577468670826155995867989841989392416473967554675927069569216523691)"},
        { 1000000, "(7686866163780120756504704687108787598650652185649163569056142218702518519446,3906118672583014628968951493493328376867199397126296469652067564264383995501)"},
        { 100000, "(4025341918271974808785600589734026802648663429991646260424633790951720792181,12092188261120553063009105063256001824310699929784774108735665676444962251470)"},
        { 100, "(2918313335010035498813915079123510985165681786090089918491607260013248937162,17618193400648240490544433854439160718346256852721163860896583680851810317370)"},
        { 10, "(9210363875285848131508656139203480860760390189260733424246762951552226162107,5835599852467557952169717222312357036255412228954185370294061936827699274723)"},
        { 2, "(11214335563238929817636333449692455681168983913051128543176696602685434833110,19682196534069234051611777524958224662412865849089027703776500470935035076258)"},
        { 0, "" }
    };

    int index = 0;
    while (results[index].n && results[index].n != n) {
        ++index;
    }
    return results[index].result;
}

int main(int argc, char **argv) 
{
    MultiExpG1 meg1;

    meg1.run(argc, argv);

    return EXIT_SUCCESS;
}
    
