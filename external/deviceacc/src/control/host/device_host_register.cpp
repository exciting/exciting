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
///  @file
///  This file contains a host-device register, which can in principle can be used to
///  create a fine map host-device map address (even for derived types). 
///  I have exposed C functions for Fortran interfacing.
///  Limits: each host must have its own register, though it can connect to several devices.
///  Limits: does not work with Unified Shared Memory (USM)
///  Dependencies: OMP4.5 or higher
///  This is the CPU backend
///
#include <unordered_map>
#include <tuple>
#include <string>
#include <stdexcept>
#include <iostream>

/// @brief for: {Device id, size, host_ptr, device_ptr}
using ref_info      = std::tuple<std::size_t, std::size_t, void*, void*>;
/// @brief key for the dictionary
using identifier    = std::string;
/// @brief the register type is nothing else than a hash map
using register_type = std::unordered_map<identifier, ref_info>;
/// @brief the index for device id for give host
const std::size_t idx_device_id  = 0;
/// @brief the index for object size in bytes
const std::size_t idx_obj_size   = 1;
/// @brief  the index for the device pointer
const std::size_t idx_device_ptr = 2;
/// @brief  the index for the host pointer
const std::size_t idx_host_ptr   = 3;

class device_host_register {
private:
    /// Internal register
    register_type rg;
    /// Host id, this is obtained through OMP call
    const int host_id;
public:
    /// Delete the copy constructor.
    device_host_register(const device_host_register& source) = delete;
    /// Delete assignment operator.
    device_host_register& operator=(const device_host_register& source) = delete;

    /// Inits the register (the host identity is obtained through an omp call)
    device_host_register() :  host_id(0) {}

    /// When out of scope or explicitly deallocated the pointer in C compliant calls
    /// It checks the registered addresses and cleans them
    ~device_host_register() {
        for (auto &entry : this->rg) {
            if (std::get<idx_host_ptr>(entry.second) != nullptr) {
                std::free(std::get<idx_host_ptr>(entry.second));
            }
        }
    }

    /// Allocs memory on the device.
    /// In the CPU backend this only sets the size
    /// @param[in] id        - name identifying the variable
    /// @param[in] size      - the bytes that should be allocated in the device
    /// @param[in] device_id - device id in which the memory should be allocated
    void alloc_device(identifier id, std::size_t size, std::size_t device_id) {
        if (rg.count(id) != 0) {
            throw std::runtime_error("Error in alloc_device: " + id + " is already in use");
        }

        rg[id] = ref_info(device_id, size, nullptr, nullptr);
    }

    /// Associates a host pointer to a device pointer
    /// In the CPU backend this simply do nothing but fill the device and host vector with the same
    /// @param[in] id        - name identifying the variable
    /// @param[in] host_ptr  - the pointer in the host which should be associated to the id
    void associate(identifier id, void* host_ptr) {
        // Check existence
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in associate: " + id + " is not in the list");
        }
        // Check that host is not nullptr
        if (host_ptr == nullptr) {
            throw std::runtime_error("Error in associate: " + id + " the host pointer is c_null_ptr");
        }
        auto& mrg = rg[id];
        // Check if the device pointer is already associated
        if (std::get<idx_device_ptr>(mrg) != nullptr) {
            throw std::runtime_error("Error in associate: " + id + " device pointer is already associated to a host pointer");
        }
        std::get<idx_device_ptr>(mrg) = host_ptr;
    }

    /// Diassociates a host pointer to a device pointer
    /// In the CPU backend does nothing
    /// @param[in] id        - name identifying the variable which should be disassociated
    void disassociate(identifier id) {}

    /// Send data between the host and the device, namely between memory addresses associated to the id
    /// In the CPU backend this does nothing
    /// @param[in] id        - name identifying the variable that should be transferred from the host to the device
    void host_to_device(identifier id) {}

    /// Send data between the device and the host, namely Between memory addresses associated to the id
    /// In the CPU backend this does nothing
    /// @param[in] id  - name identifying the variable that should be transferred from the device to the host
    void device_to_host(identifier id) {}

    /// Removes id from the register, it cleans and disassociate if required
    /// @param[in] id        - name identifying the variable that should be transferred from the host to the device
    void remove(identifier id) {
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in remove: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        if (std::get<idx_host_ptr>(mrg) != nullptr) {
            std::free(std::get<idx_host_ptr>(mrg));
        }
        rg.erase(id);
    }

    /// Get device pointer associated to an entry
    /// In CPU backend, simply provides address to host memory
    /// Moreover, in case no allocation is found it is done here (there is a critical section to ensure thread safety)
    /// @param[in] id        - name identifying the variable for which the device pointer will be retrieved
    void* get_device_ptr(identifier id) {
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in get_device_ptr: " + id + " is not in the list");
        }
        auto &mrg = rg[id];

        #pragma omp critical
        {
            if (std::get<idx_device_ptr>(mrg) == nullptr) {
                std::get<idx_device_ptr>(mrg) = std::malloc(std::get<idx_obj_size>(mrg));
                std::get<idx_host_ptr>(mrg) = std::get<idx_device_ptr>(mrg);
            }
        }

        return std::get<idx_device_ptr>(mrg);
    }

};

/// Expose C++ class and members for Fortran interfacing
extern "C" {

    void _constructor_device_host_register(device_host_register** rg) {
        try {
            *rg = new device_host_register;
        } catch (const std::bad_alloc& e) {
            throw std::runtime_error("Error in create_register : allocation of the register object failed");
        }
    }

    void _destructor_device_host_register(device_host_register* rg) {
        delete rg;
    }

    void _alloc_device_device_host_register(device_host_register* rg, const int id_len, const char* id, std::size_t size, int device_id) {
        std::string cxx_id(id, id_len);
        rg->alloc_device(cxx_id, size, device_id);
    }

    void _associate_device_device_host_register(device_host_register* rg, const int id_len, const char* id, void* host_ptr) {
        std::string cxx_id(id, id_len);
        rg->associate(cxx_id, host_ptr);
    }

    void _disassociate_device_device_host_register(device_host_register* rg, const int id_len, const char* id) {
        std::string cxx_id(id, id_len);
        rg->disassociate(cxx_id);
    }

    void _host_to_device_device_device_host_register(device_host_register* rg, const int id_len, const char* id){
        std::string cxx_id(id, id_len);
        rg->host_to_device(cxx_id);
    }

    void _device_to_host_device_device_host_register(device_host_register* rg, const int id_len, const char* id){
        std::string cxx_id(id, id_len);
        rg->device_to_host(cxx_id);
    }

    void _remove_device_device_host_register(device_host_register* rg, const int id_len, const char* id){
        std::string cxx_id(id, id_len);
        rg->remove(cxx_id);
    }
    void*  _get_device_ptr_device_device_host_register(device_host_register* rg, const int id_len, const char* id){
        std::string cxx_id(id, id_len);
        auto d_ptr = rg->get_device_ptr(cxx_id);
        return d_ptr;
    }

}

