//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <string>
//#include <stdlib.h>
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
//#include <chrono>
//#include <algorithm>
//
//const size_t NUM_OF_PEPTIDES = 2000000; // number of peptides in document 
//const size_t MAX_LENGTH_PEPTIDE = 10;
//const size_t SIZE_PEPTIDES = NUM_OF_PEPTIDES * MAX_LENGTH_PEPTIDE;
//const size_t SIZE_SPECTRA = 20;
//const size_t SIZE_INNER_PRODUCT = NUM_OF_PEPTIDES; // ip is counted for each peptide
//
//const size_t BLOCK_SIZE = 1024;
//const size_t GRID_SIZE = NUM_OF_PEPTIDES / BLOCK_SIZE;
//
//__constant__ int spectra_const[SIZE_SPECTRA];
//__global__ void kernel(int* peptides, int * ip, int length)
//{
//    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    int sum = 0;
//    for (int j = 0; j < length; j++)
//    {
//        sum += spectra_const[peptides[i * length + j]];
//    }
//    ip[i] = sum;
//}
//
//int main(int argc, char* argv[])
//{
//    auto all_begin = std::chrono::steady_clock::now();
//    std::ofstream time_res("time_results.txt", std::ios::app);
//
//    // Write the time of execution
//    std::time_t t = std::time(0);
//    std::tm* now = std::localtime(&t);
//    time_res
//        << "GLOB 2: "
//        << (now->tm_year + 1900) << '-'
//        << now->tm_mday << ' '
//        << now->tm_hour << ':'
//        << now->tm_min << ':'
//        << now->tm_sec << std::endl;
//
//    // Read files and fill arrays
//    std::fstream spectra_file("../../spectra.txt");
//    std::fstream peptides_file("../../2m_size.txt");
//    std::ofstream inner_product_file("results.txt");
//
//    // Allocation spectra and peptides on Host
//    int* h_spectra, * h_peptides, * inner_product;
//    cudaError_t spectra_stat = cudaMallocHost((void**)&h_spectra, SIZE_SPECTRA * sizeof(int));
//    cudaError_t peptides_stat = cudaMallocHost((void**)&h_peptides, SIZE_PEPTIDES * sizeof(int));
//    cudaError_t ip_stat = cudaMallocHost((void**)&inner_product, SIZE_INNER_PRODUCT * sizeof(int));
//   
//    // Reading data
//    for (int i = 0; i < SIZE_SPECTRA; ++i) {
//        int spec_val;
//        spectra_file >> spec_val;
//        h_spectra[i] = spec_val;
//    }
//
//    for (int i = 0; i < NUM_OF_PEPTIDES; i++)
//    {
//        for (int j = 0; j < MAX_LENGTH_PEPTIDE; ++j) {
//            int pep_val;
//            peptides_file >> pep_val;
//            h_peptides[i * MAX_LENGTH_PEPTIDE + j] = pep_val;
//        }
//    }
//
//    // Copy spectra to const memory
//    cudaMemcpyToSymbol(spectra_const, h_spectra, SIZE_SPECTRA * sizeof(int));
//
//    // Streams configuration
//    const int nStreams = 4;
//    const int streamSize = SIZE_PEPTIDES / nStreams;
//    const int streamBytes = streamSize * sizeof(int);
//    cudaStream_t stream[nStreams];
//
//    // Device mem configuration
//    int* d_1, * d_2; // *d_3, * d_4;         // -- Peptides --
//    // int* d_s1, * d_s2, * d_s3, * d_s4;     // -- Spectras --
//    int* d_ip1, * d_ip2; // * d_ip3, * d_ip4; // -- Ips --
//
//    // --- 1 ----
//    // Creating sequence of operations
//    // cudaMalloc((void**)&d_s1, SIZE_SPECTRA * sizeof(int));
//    cudaMalloc((void**)&d_1, streamBytes);
//    cudaMalloc((void**)&d_ip1, (NUM_OF_PEPTIDES) * sizeof(int));
//    
//    int offset = 0;
//    cudaMemcpyAsync(&d_1, &h_peptides, streamBytes, cudaMemcpyHostToDevice, stream[0]);
//    /*cudaMemcpyAsync(&d_s1, &h_spectra, SIZE_SPECTRA * sizeof(int), cudaMemcpyHostToDevice, stream[0]);*/
//    kernel << <streamSize/BLOCK_SIZE, BLOCK_SIZE, 0, stream[0] >>> (d_1, d_ip1, MAX_LENGTH_PEPTIDE);
//    
//    // Transmit IP to host
//    std::vector<int> h_ip1(SIZE_INNER_PRODUCT / nStreams);
//    cudaMemcpy(h_ip1.data(), d_1, (SIZE_INNER_PRODUCT / nStreams) * sizeof(int), cudaMemcpyDeviceToHost);
//    
//  
//    // --- 2 ----
//    // Creating sequence of operations
//    // cudaMalloc((void**)&d_s1, SIZE_SPECTRA * sizeof(int));
//    cudaMalloc((void**)&d_2, streamBytes);
//    cudaMalloc((void**)&d_ip2, (SIZE_INNER_PRODUCT/ nStreams) * sizeof(int));
//
//
//    offset = SIZE_INNER_PRODUCT / nStreams;
//    cudaMemcpyAsync(&d_2, &h_peptides[offset], streamBytes, cudaMemcpyHostToDevice, stream[1]);
//    /*cudaMemcpyAsync(&d_s1, &h_spectra, SIZE_SPECTRA * sizeof(int), cudaMemcpyHostToDevice, stream[0]);*/
//    kernel << <streamSize / BLOCK_SIZE, BLOCK_SIZE, 0, stream[1] >> > (d_2, d_ip2, MAX_LENGTH_PEPTIDE);
//
//
//    // Transmit IP to host
//    std::vector<int> h_ip2(SIZE_INNER_PRODUCT / nStreams);
//    cudaMemcpy(h_ip2.data(), d_2, (SIZE_INNER_PRODUCT / nStreams) * sizeof(int), cudaMemcpyDeviceToHost);
//    
//
//
//    // Save IP in the file
//    for (int i = 0; i < SIZE_INNER_PRODUCT / nStreams; i++)
//        inner_product_file << h_ip1[i] << std::endl;
//
//    // Save IP in the file
//    for (int i = 0; i < SIZE_INNER_PRODUCT / nStreams; i++)
//        inner_product_file << h_ip2[i] << std::endl;
//
//    cudaFree(d_1);
//    cudaFree(d_ip1);
//    cudaFree(d_2);
//    cudaFree(d_ip2);
//
//    cudaFreeHost(h_peptides);
//    cudaFreeHost(h_spectra);
//
//    auto all_end = std::chrono::steady_clock::now();
//    time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
//        << "----------------" << std::endl;
//
//    return 0;
//}