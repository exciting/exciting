/// Copyright (C) 2024 DEVICEACC developers
///
/// Licensed under the Apache License, Version 2.0 (the "License");
/// you may not use this file except in compliance with the License.
/// You may obtain a copy of the License at
///
///   http://www.apache.org/licenses/LICENSE-2.0
///
/// Unless required by applicable law or agreed to in writing, software
/// distributed under the License is distributed on an "AS IS" BASIS,
/// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
/// implied. See the License for the specific language governing
/// permissions and limitations under the License.
///
/// @file
/// @brief Contains C++ functions for generating batched arrays required by batching routines.
///
/// This file provides C++ utility functions to set up pointers for batched array data.
/// The functions assume that the data is organized contiguously, with each batch segment
/// representing a distinct, sequential portion of the array.
///
/// @note This file uses `extern "C"` to prevent name mangling of C++.

#include <cstddef>
#include <memory>

extern "C" {

    /// @brief Generates pointers for batched routines.
    ///
    /// This function initializes an array of pointers, `batched_ptr`, for use in batched operations.
    /// Each pointer in `batched_ptr` corresponds to a batch in `data`, assuming a contiguous layout
    /// where the batch identifier is the slowest (outermost) index.
    ///
    /// @param[in] data           Pointer to the contiguous block of data.
    /// @param[in] nbatch         Number of batches.
    /// @param[in] batch_byte_size Size (in bytes) of each batch.
    /// @param[out] batched_ptr   Array of pointers, where each entry points to the start of a batch in `data`.
    ///
    /// @details
    /// Given that `data` is arranged with the batch as the outermost index, this function
    /// computes the starting address of each batch by offsetting from the base pointer `data`.
    /// For each batch `ibatch`, `batched_ptr[ibatch]` is set to point to `data + ibatch * batch_byte_size`.
    ///
    void generate_batched_array(void* data, std::size_t nbatch, std::size_t batch_byte_size, void** batched_ptr) {
        for (std::size_t ibatch = 0; ibatch < nbatch; ++ibatch) {
            batched_ptr[ibatch] = reinterpret_cast<unsigned char*>(data) + ibatch * batch_byte_size;
        }
    }
}
