# TODO and NOTES

- Size of memory, chunks and his organization, it's very important. For example, with threads, if we active light operation counters, real time is 3.5 times higher, because all threads access in read/write to the same region, and it's bottleneck.
  
- Make optimizations and measurements on big machine as ryzen, better on amd to use all available cores.
  
- Start from low level parts, as batchInverse to detect the best performance way. If source array an destination array are different, not need allocate memory for temporal array (this temporal array could send by the caller, to reuse it)
  
- Comparative measurements should never be done in cold, need to do first access to load data in cache. For example, at beginning we did 4 tests (add, sub, mmul, square), but we did't do a previous useless reading, the result was that the add operation (the first test) was the less performance.

- To known real cost of operation, do it using same memory region (same index), thus the performance losses due to memory accesses will be insignificant.

- Increase memory space to detect what is the optimal space with current machine cache level sizes. L1 and L2 cache memory was by core, but L3 was common.

- Be careful with copy operations, or call that use internally a copy, because copy to pass from a Point to a AffinePoint implies an expensive division.
  
- TODO: self-test to detect best performance sizes, these tests could be used to compare servers to avoid surprises.
  
- TODO: measure and compare different versions of batchInverse and perhaps call one or another, in function of parameters (size, destination =? source)

- TODO: measure and compare different versions, passing elements in differents arrays (\*a, \*b, \*c) where a[i] is far b[i], or using structure arrays (\*{a,b,c}) where a[i] is near b[i].

## Base operations count on 128M of adds (multiexp)
  
  |operation|add-by-add |batch|%|
  |---:|---:|---:|---:|
  |add|0|128,000,000||
  |sub|896,000,000|640,000,000||
  |**add+sub**|**896,000,000**|**768,000,000**|85.71%|
  ||||
  |mmul|1,024,000,000|639,999,997||
  |square|256,000,000|128,000,000||
  |**mmul+square**|**1,280,000,000**|**767,999,997**|60.00%|


  |operation|add by add|batch|
  |---:|---:|---:|
  |**theoretical muls**|10|6|
  ||**1,280,000,000**|**768,000,000**|

  The ratio of finite field multiplication on add-by-add method and batch method, is the same from the empirical and theoretical point of view.

## Cost of basic operations
Results was in seconds to do 128M of basic operations

MTA = memory acces time

  |operation|wo/MTA|w/MTA|MTA|%MTA|
  |---:|---:|---:|---:|---:|
  |add|0.2696|1.8518|1.5822|85.44%
  |sub|0.2453|1.8526|1.6073|86.75%
  |mmul|1.6613|2.2868|0.6255|27.35%
  |square|1.6234|2.1499|0.5375|25.00%

Perhaps, how add and sub are fast, when them have finished, memory isn't available yet.

## Conclusion

- With **with MTA** times and count of operation of the first table, we calculate that minimum time of batch method was 69.45% of add-by-add, **means a reduction near 30%.**

- With **without MTA** times and count of operation of the first table, we calculate that minimum time of batch method was 62.61% of add-by-add, **means a reduction near 37%.**

- These reductions was considered in non parallel situation, in parallel situation times depends that how we could divide data throught threads.

- For batch method, it's necessary **build a structure to collect add operations**. I first measurements, this part take near of 30% of total time, it's a bottleneck.

## Sources 
- **c/multiexp_ba.cpp/.hpp** has been implementation of batch method, in this way we could use add-by-add method or batch method. In curve.hpp was defined multiMulByScalarBa to call batch method. 
- **benchmark/curve_adds.cpp** has been implemented to make performance tests, in this file has been implemented base operations test.
- **benchmark/multiexp_g1.cpp** has been updated, to allow make multiexp tests.
- **batch_accumulators.cpp/.hpp** has been implemented a class to collect adds, instances of this class are used by multiexp_ba.
- **batch_accumulators_test.cpp** has been implemented some tests using google test library.

**NOTICE:** be carefull with counters (**COUNT_OPS**) specially in multithreading, because it seriously affects performance.