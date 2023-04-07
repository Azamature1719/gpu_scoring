#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <chrono>
#include <algorithm>

const size_t NUM_OF_PEPTIDES = 2000000;  //number of peptides in document 
const size_t MAX_LENGTH_PEPTIDE = 12;
const size_t SIZE_PEPTIDES = NUM_OF_PEPTIDES * MAX_LENGTH_PEPTIDE;
const size_t SIZE_SPECTRA = 20;
const size_t SIZE_INNER_PRODUCT = NUM_OF_PEPTIDES; //ip is counted for each peptide

const size_t BLOCK_SIZE = 1024;
const size_t GRID_SIZE = NUM_OF_PEPTIDES / BLOCK_SIZE;

__constant__ int spectra_const[SIZE_SPECTRA];
__global__ void kernel(const int* peptides, int* inner_product, int length)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    for (int j = 0; j < length; j++)
    {
        inner_product[i] += spectra_const[peptides[i * length + j]];
    }
}

int main(int argc, char* argv[])
{
    auto all_begin = std::chrono::steady_clock::now();
    std::ofstream time_res("time_results.txt", std::ios::app);

    // Write the time of execution
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
    time_res
        << "GLOB 2: "
        << (now->tm_year + 1900) << '-'
        << now->tm_mday << ' '
        << now->tm_hour << ':'
        << now->tm_min << ':'
        << now->tm_sec << std::endl;

    // Read files and fill arrays
    auto read_beg = std::chrono::steady_clock::now();

    std::fstream spectra_file("../spectra.txt");
    std::fstream peptides_file("../2m_pep.txt");
    std::ofstream inner_product_file("results.txt");

    // Allocation memory on host (spectra, peptides)
    auto malloc_host_beg = std::chrono::steady_clock::now();

    int* h_spectra = new int[SIZE_SPECTRA];
    int* h_peptides = new int[SIZE_PEPTIDES];
    
    auto malloc_host_end = std::chrono::steady_clock::now();
    time_res << "MallocH: " << std::chrono::duration_cast<std::chrono::milliseconds>(malloc_host_end - malloc_host_beg).count() << std::endl;

    // Read data from files
    for (int i = 0; i < SIZE_SPECTRA; ++i) {
        int spec_val;
        spectra_file >> spec_val;
        h_spectra[i] = spec_val;
    }

    for (int i = 0; i < NUM_OF_PEPTIDES; i++)
    {
        for (int j = 0; j < MAX_LENGTH_PEPTIDE; ++j) {
            int pep_val;
            peptides_file >> pep_val;
            h_peptides[i * MAX_LENGTH_PEPTIDE + j] = pep_val;
        }
    }

    auto read_end = std::chrono::steady_clock::now();
    time_res << "Read: " << std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_beg).count() << std::endl;

     // Allocatte memory on the device
     int * d_peptides, * d_inner_product;
     auto malloc_beg = std::chrono::steady_clock::now();
     cudaMalloc((void**)&d_peptides, SIZE_PEPTIDES * sizeof(int));
     auto malloc_end_peptides = std::chrono::steady_clock::now();
     time_res << "Malloc end peptides: " << std::chrono::duration_cast<std::chrono::milliseconds>(malloc_end_peptides - malloc_beg).count() << std::endl;
     cudaMalloc((void**)&d_inner_product, SIZE_INNER_PRODUCT * sizeof(int));
     auto malloc_end = std::chrono::steady_clock::now();
     time_res << "Malloc end ip: " << std::chrono::duration_cast<std::chrono::milliseconds>(malloc_end - malloc_end_peptides).count() << std::endl;

     // Copy spectra to const memory
     auto htc_beg = std::chrono::steady_clock::now();
     cudaMemcpyToSymbol(spectra_const, h_spectra, SIZE_SPECTRA * sizeof(int));
     auto htc_end = std::chrono::steady_clock::now();
     time_res << "SpecHTC: " << std::chrono::duration_cast<std::chrono::milliseconds>(htc_end - htc_beg).count() << std::endl;

     // Copy peptides to global memory
     auto pep_copy_beg = std::chrono::steady_clock::now();
     cudaMemcpy(d_peptides, h_peptides, SIZE_PEPTIDES * sizeof(int), cudaMemcpyHostToDevice);
     auto pep_copy_end = std::chrono::steady_clock::now();
     time_res << "PepHTD: " << std::chrono::duration_cast<std::chrono::milliseconds>(pep_copy_end - pep_copy_beg).count() << std::endl;

     // Count IP
     auto kernel_begin = std::chrono::steady_clock::now();
     kernel << <GRID_SIZE, BLOCK_SIZE >> > (d_peptides, d_inner_product, MAX_LENGTH_PEPTIDE);
     auto kernel_end = std::chrono::steady_clock::now();
     time_res << "CountIp: " << std::chrono::duration_cast<std::chrono::microseconds>(kernel_end - kernel_begin).count() << std::endl;

     // Transmit IP to host
     std::vector<int> h_inner_product(SIZE_INNER_PRODUCT);
     auto ip_copy_beg = std::chrono::steady_clock::now();
     cudaMemcpy(h_inner_product.data(), d_inner_product, SIZE_INNER_PRODUCT * sizeof(int), cudaMemcpyDeviceToHost);
     auto ip_copy_end = std::chrono::steady_clock::now();
     time_res << "IpDTH: " << std::chrono::duration_cast<std::chrono::milliseconds>(ip_copy_end - ip_copy_beg).count() << std::endl;

     // Save IP in the file
     auto to_res_beg = std::chrono::steady_clock::now();
     for (int i = 0; i < SIZE_INNER_PRODUCT; i++)
         inner_product_file << h_inner_product[i] << std::endl;
     auto to_res_end = std::chrono::steady_clock::now();
     time_res << "ResFile: " << std::chrono::duration_cast<std::chrono::milliseconds>(to_res_end - to_res_beg).count() << std::endl;

     cudaFree(d_peptides);
     cudaFree(d_inner_product);
    
     delete[] h_peptides;
     delete[] h_spectra;

     auto all_end = std::chrono::steady_clock::now();
     time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
     << "----------------" << std::endl;

    return 0;
}
