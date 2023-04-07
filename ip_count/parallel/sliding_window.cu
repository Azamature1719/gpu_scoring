#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <chrono>
#include <limits>
#include <algorithm>

const int ALIGNMENT_ELEM = -1;

const int delta = 10000; // Difference between starting points
const int match_range = 200000; // Number of peptides that should be matched with a single spectra

// -- To develop a custom configurable program -- 
const int number_of_launches = 56; // Calculated manually
const int spectra_at_once = 180; // Calculated manually

__global__ void ip_count(int* peptides, const int* spectra, int* inner_product, int peptide_length, int launch_id) {
 
    int t_id = blockDim.x * blockIdx.x + threadIdx.x;
    int spec_id = t_id + launch_id * spectra_at_once; // Order number of spectra

    __shared__ int s[20];
    for(int i = 0; i < 20; ++i){
	s[i] = spectra[spec_id * 20 + i];	
    }   
    __syncthreads();

    if (spec_id < 10000 && t_id < 180) {
        for (int match_pep_id = 0; match_pep_id < match_range; ++match_pep_id) { // Num of the peptide in a 200K peptide sequence
            int match_pep_val_id = 0; // Order num of the peptide value in the concrete peptide 
            int cur_pep_val = t_id * delta * peptide_length + match_pep_id * peptide_length + match_pep_val_id; // Order num of the peptide value in the peptides array. t_id is used because peptides are updated each launch
            int ip_val_id = spec_id * match_range + match_pep_id;
            inner_product[ip_val_id] = 0;
            while (peptides[cur_pep_val] != ALIGNMENT_ELEM) {
                int spec_value = s[peptides[cur_pep_val]];
                inner_product[ip_val_id] += spec_value;
                ++match_pep_val_id;
                cur_pep_val = t_id * delta * peptide_length + match_pep_id * peptide_length + match_pep_val_id;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    auto all_begin = std::chrono::steady_clock::now();
    std::ofstream time_res("time_results.txt", std::ios::app);
    std::ofstream res("results.txt");

    cudaSetDevice(0); // RTX 3090ti

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

    // Get vars from command line
    //std::string spectra_filename = "../../" + std::string(argv[1]);
    //std::string peptides_filename = "../../" + std::string(argv[2]);
    //int peptide_length = std::atoi(argv[3]);
    //int block_size = std::atoi(argv[4]);

    std::string spectra_filename = "../../10k_spectra.txt";
    std::string peptides_filename = "../../2m_pep.txt";
    int peptide_length = 13;

    std::fstream spectra_file(spectra_filename);
    std::fstream peptides_file(peptides_filename);

    std::vector<int> h_peptides;
    std::vector<int> h_spectra;

    //// Fill spectras
    /*size_t max_size_spectra = 0;*/
    std::string temp;
    //while(std::getline(spectra_file, temp)){
    //    std::vector<int> cur_spectra;
    //    std::stringstream row(temp);
    //    int spec_val;
    //    while (row >> spec_val) {
    //        cur_spectra.push_back(spec_val);
    //    }
    //    if (cur_spectra.size() > max_size_spectra)
    //        max_size_spectra = cur_spectra.size();
    //    h_spectra.push_back(cur_spectra);
    //}
    //size_t num_of_spectra = h_spectra.size();

    auto read_begin = std::chrono::steady_clock::now();

    // Fill spectra
    int num_of_spectra = 0;
    while (std::getline(spectra_file, temp)) {
        int spec_val;
        std::stringstream row(temp);
        while (row >> spec_val) {
            h_spectra.push_back(spec_val);
        }
        ++num_of_spectra;
    }

    // Fill peptides
    size_t num_of_peptides = 0;
    while (std::getline(peptides_file, temp)) {
        int pep_size = 0;
        int pep_val;
        std::stringstream row(temp);
        while (row >> pep_val) {
            h_peptides.push_back(pep_val);
            ++pep_size;
        }
        for (; pep_size < peptide_length; ++pep_size) {
            h_peptides.push_back(ALIGNMENT_ELEM);
        }
        ++num_of_peptides;
    }
    auto read_end = std::chrono::steady_clock::now();
    time_res << "Read: " << std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_begin).count() << std::endl;

    size_t records_num = match_range * num_of_spectra;

    // Allocate memory on the device
    size_t size_peptides = h_peptides.size() * sizeof(int);
    size_t size_spectra = h_spectra.size() * sizeof(int);
    size_t size_inner_product = records_num * sizeof(int);

    auto malloc_beg = std::chrono::steady_clock::now();
    int* d_spectra, * d_peptides, * d_inner_product;
    cudaMalloc((void**)&d_peptides, size_peptides);
    cudaMalloc((void**)&d_spectra, size_spectra);
    cudaMalloc((void**)&d_inner_product, size_inner_product);
    cudaMemcpy(d_peptides, h_peptides.data(), size_peptides, cudaMemcpyHostToDevice);
    cudaMemcpy(d_spectra, h_spectra.data(), size_spectra, cudaMemcpyHostToDevice);

    auto malloc_end = std::chrono::steady_clock::now();
    time_res << "Malloc: " << std::chrono::duration_cast<std::chrono::milliseconds>(malloc_end - malloc_beg).count() << std::endl;

    auto ann_beg = std::chrono::steady_clock::now();

    // Spectrum annotation algorithm
    int block_size = 32;
    size_t grid_size = spectra_at_once / block_size + 1; // Rounding up to get results from all the threads 
    std::vector<int> h_inner_product(records_num);

    for (int launch_id = 0; launch_id < number_of_launches; ++launch_id) {
        ip_count <<<grid_size, block_size >>> (d_peptides, d_spectra, d_inner_product, peptide_length, launch_id);
	cudaMemcpy(d_peptides, h_peptides.data(), size_peptides, cudaMemcpyHostToDevice);
        cudaError er = cudaDeviceSynchronize();
        std::cout << er << "\n";
    }
    cudaMemcpy(h_inner_product.data(), d_inner_product, size_inner_product, cudaMemcpyDeviceToHost);

    auto ann_end = std::chrono::steady_clock::now();
    time_res << "Annotation: " << std::chrono::duration_cast<std::chrono::milliseconds>(ann_end - ann_beg).count() << std::endl;
   
    for(int i = 0; i < 200000; ++i){
        res << h_inner_product[i] << "\n";
    }

    cudaFree(d_spectra);
    cudaFree(d_peptides);
    cudaFree(d_inner_product);

    auto all_end = std::chrono::steady_clock::now();
    time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
        << "----------------" << std::endl;

    return 0;
}
