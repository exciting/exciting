include ../../../make.inc

.SUFFIXES: .cpp .cu

NVCC = nvcc

OBJ = lapw.o lapw_init.o lapw_load_global.o lapw_load_species.o lapw_set_sv.o \
      lapw_load_kpoint.o lapw_seceqn_init.o lapw_fft.o lapw_eigen_states.o \
      timer.o lapw_timers.o lapw_execute.o lapw_density.o linalg_cpu.o
      
ifneq (,$(findstring D_GPU_,$(CPP_OPTS)))
	OBJ := $(OBJ) gpu_interface.o
endif
	
.cpp.o:
	$(CXX) $(CXX_OPTS) -c $<

.cu.o:
	$(NVCC) -DADD_ -I/project/s299/magma/magma-1.1_my-xk6-gcc-4.6.2-libsci/include -c $<

all: $(OBJ)
	ar -r libapw.a $(OBJ)

clean:
	rm *.o *.a
