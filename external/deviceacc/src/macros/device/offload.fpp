! Copyright (C) 2024 DEVICEACC developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied. See the License for the specific language governing
! permissions and limitations under the License.

! This include file set up macros for device offloading
! Be aware that as opposed to the elements in the register
! this does not directly work with derived types elements
! except if you hide those by an association 

use omp_lib

! This macro instructs the compiler to create in the device
! a counterpart of X, where X is a procedure and/or a 
! module variable within the module scope

#define DECLARE_IN_DEVICE(X) !$omp declare target (X)

! The same but declares the object itself
#define DECLARE_THIS_IN_DEVICE !$omp declare target

! This macro maps X into the device by allocating the memory but without
! memory transfer
#define DEVICE_MAP_ALLOC(X) !$omp target enter data map(always, alloc: X)

! This macro maps X into the device by allocating the memory and copying
! data from the host to the device
#define DEVICE_MAP_TO(X) !$omp target enter data map(always, to: X)

! This macro releases an association if it was the last association
! it also deletes the object from the device memory
#define DEVICE_MAP_RELEASE(X) !$omp target exit data map(always, release: X)

! This macro deletes the associated object from the device
#define DEVICE_MAP_DELETE(X) !$omp target exit data map(always, delete: X)

! This macro transfers data from host to device
#define DEVICE_UPDATE_TO(X) !$omp target update to(X)

! This macro transfers data from device to host
#define DEVICE_UPDATE_FROM(X) !$omp target update from(X)

! This returns the team number
#define DEVICE_GET_TEAM_ID omp_get_team_num() + 1

! This macro inits the target block, that is offloadable 
#define DEVICE_BEGIN_BLOCK !$omp target

! This macro ends the target block, that is offloadable
#define DEVICE_END_BLOCK !$omp end target

! This macro inits the target block, that is offloadable with a device ptr
#define DEVICE_BEGIN_BLOCK_HAS_DEVICE_ADDR(X) !$omp target has_device_addr(X)

! This macro creates a parallel workshare
#define PARALLEL_WORKSHARE !$omp parallel workshare

! This macro ends the parallel workshare
#define END_PARALLEL_WORKSHARE !$omp end parallel workshare

! This macro gives the current thread id 
#define DEVICE_GET_THREAD_ID omp_get_thread_num() + 1

! This macro returns the number of teams
#define DEVICE_GET_NUM_TEAMS omp_get_num_teams()

! This macro returns the number of threads
#define DEVICE_GET_NUM_THREADS omp_get_num_threads()

! This macro indicates that X can be executed asynchronously
#define ASYNCHRONOUS(X) X nowait

! This macro establishes that X has a write-only dependency
#define WRITE_DEPENDENCY(X) depend(out: X)

! This macro establishes a read-only dependency, i.e. it will wait write dependencies on X to finish
#define READ_DEPENDENCY(X)  depend(in: X)

! This macro establishes a read/write dependency on X (no mutex exist between tasks of the same construct)
#define READ_WRITE_DEPENDENCY(X)  depend(inout: X)

! This macro establishes a read/write dependency on X (there is a mutex between tasks of the same construct)
#define MUTEX_DEPENDENCY(X)  depend(mutexinoutset: X)

! This macro syncronizes the device tasks
#define DEVICE_OMP_KERNELS_SYNCHRONIZE !$omp taskwait

! This macro is a safe do simd 
#define DEVICE_BEGIN_THREAD_WORK !$omp do simd

! This macro is a safe end do simd
#define DEVICE_END_THREAD_WORK !$omp end do simd

! This macro is to add safely add conditionals to the macros
#define WHEN(X) if(X)
