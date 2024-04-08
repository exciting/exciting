# Local components
if(NOT CPUBACKEND)
    if (AMD)
        set(SRC_FFT ${PROJECT_SOURCE_DIR}/external/hipfort/hipfort_rocfft_enums.f90
                    ${PROJECT_SOURCE_DIR}/external/hipfort/hipfort_rocfft.f90)    
    endif() 
    set(SRC_FFT ${SRC_FFT} ${CMAKE_CURRENT_SOURCE_DIR}/device/fft_device_t.f90)
else()
    set(SRC_FFT ${CMAKE_CURRENT_SOURCE_DIR}/host/fft_device_t.f90)
endif()
# Push back to SRC
set(SRC ${SRC} ${SRC_FFT} PARENT_SCOPE)
