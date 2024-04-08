# Deviceacc
**Device acc**elerated is a modern Fortran library presenting developers with tools to develop their device accelerated
codes independently from the major vendors (NVIDIA, AMD, and Intel). The library is MPI aware, associating each MPI-process to
an independent device. Additionally, it provides a CPU-only backend, so no modification is needed when compiling for systems without
devices.

# Supported operations 

* __Host-device management__: The code offers a __device_world_t__ which performs proper association between host and device (i.e. associates each
  MPI process to a device). Moreover, this class provides __register__ object, which allows to control the memory allocation in the device, its association
  with host elements, as well as the transfer. To do so, each element is associated to a user-defined __tag__ (string) which uniquely represents it in the given MPI process. 
* __Linear algebra__: basic linear algebra are supported (not all operations are supported, further work is done here to extend the number of supported operations).
* __FFT__: Fast Fourier Transforms complex-complex transforms are supported in both single and double precision, through __fft_device_t__ class.
# Dependencies:

The dependencies vary depending on the vendor.

# Limits

The number of devices must match that of the devices, otherwise the performance can be severely affected. Consequently,
the code fails under any other condition.
