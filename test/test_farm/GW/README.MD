# GW Testing

## Test: GW/ZnO-PBESOL

ZnO appears to be the only system we have identified for which GW QP energies do not vary as function
of threads or MPI processes used for performing the calculation.

As such, this is the only GW test that is enabled in the test suite. Whilst this is not ideal, ZnO 
provides test coverage for much of the main GW code. One should take care to ensure that the code
they are touching is covered by this test case prior to making changes.

Results are consistent for **all** build types locally. In the CI pipeline, all GCC 8 build types 
(serial, SMP, MPIANDSMP) and Intel 2019 serial and SMP build types give consistent, reproducible numbers.
The one exception is Intel 2019 MPIANDSMP, however this gives reproducible results if the number of
OMP threads is set to 1[^1]. As such, a special case is made:

**Application test run 1:** 
* GW test is disabled in the CI by adding a GROUP tag.
* All other application tests are run with the default MPIANDSMP settings.

**Application test run 2:**
* GW test is re-enabled.
* The test suite is run only for GW, explicitly setting MKL threads _> 1_ (which appears to have no
  effect on the result), _N_ MPI processes to 2 (the default) and _N_ OMP threads to 1 (not the default).

Note, one could also:  
* a) Not run the GW test in the Intel 2019 MPIANDSMP pipeline.  This means slightly less automated coverage, but
     as MPI is not used in most of the quartic implementation, it would not be too detrimental.

* b) Or, an easier way to implement this (instead of using _sed_ and _cmp_) would be to have the test framework
    accept an argument `-except GW` which defines tests *not* to run. This would be easy to implement
    using the existing regex functions in the test framework.

Additionally, the tolerance for `n_q_cycles` was loosened, as one gets different numbers for this when running
with different numbers of MPI processes. Alex did not think too hard about this, so one should verify what is
being reported, and if this is reasonable. Perhaps the code reporting needs improving, or this key does not
need to be tested.

[^1] Alex notes that he does not understand why this test behaves fine locally, but in what 
should be an identical environment in the CI, it's not reproducible for all OMP thread settings.
This is even more unusual as the SMP build gives consistent numbers w.r.t. OMP threads.

## Other Test Cases

* Si-PBESOL
* ZrO2-PBESOL

These test cases can give consistent results if run locally, but have not proved to be reproducible in 
the CI pipeline. One should spend more time looking at their numbers when run with serial and SMP builds.

They rely on restarting the calculation from the binary STATE.OUT file, and fixing the number of MKL 
threads for the GW's SCF call to the ground state. As the latter is not possible with openBLAS, these
tests are only run with Intel MKL.

In order to use these GW tests, one needs to initially generate some benchmark data with
the input settings:

```
<groundstate
    do="fromscratch"

<gw
    skipgnd="false"
```

The GW calculation can then be repeated with the settings:

```
<groundstate
    do="skip"

<gw
    skipgnd="false"
```

In order to ensure the result is independent of the choice of threads or MPI processes,
all test executions (following generation of reference data, including `STATE.OUT`) should be run with:

```
mpirun -np <NP> ~/exciting/bin/exciting_mpismp -gw-scf-threads <THREADS>
```

where `<THREADS>` is equal to the number of threads used by MKL for the
ground state calculation. That is, in order to get QP eigenvalues that are invariant
w.r.t. the choice of threads or MPI processes, the SCF call in the GW must be performed
with the same number of threads as was used by MKL in the ground state calculation.

This is achieved in the GW by temporarily setting MKL threads to `<THREADS>` while performing
the SCF call. For the rest of the GW execution, the MKL threads will equal whatever is set
via `export MKL_NUM_THREADS`, therefore having a negligible impact on performance.

As such, we:
* a) Store STATE.OUT files
* b) Only test with MKL, as one cannot do the same thread setting in the code with
       openBLAS (see `src/src_gw/gw_scf.F0` for further discussion), and `STATE.OUT`
       is compiler-specific. One could move towards addressing this situation by
       converting `STATE.OUT` to HDF5.

To ensure consistent settings are used, the test suite parses a `spec.yml` and checks if the runtime 
settings against what the tests were generated with. This file was generated manually by Alex.

## Final Thoughts

These statements are not complete. One should investigate for what MPI processes, MKL and OMP threads the ZnO test case fails
(if any), and whether one can get consistent results with the other tests when linking to vanilla BLAS/LAPACK
(i.e. a non-threaded version). 

One could also look at the tests when using GCC + MKL.

However, one reproducible test case with tight tolerances already enables the GW developers to make significant
progress in regard to refactoring, which previously was not possible (in a systematic way).
