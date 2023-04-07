//#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <vector>
//#include <stdlib.h>
//#include <math.h>
//#include <chrono>
//#include "common_vars.cuh"
//#include <device_launch_parameters.h>
////#include "kernel.cuh"
//#include "dev_array.cu"
//
//using namespace std;
//
//__constant__ int spectra_const[SIZE_SPECTRA];
//
//__global__ void comp_ip_kernel(int* peptides, int* inner_product) {
//    int current_peptide = blockDim.x * blockIdx.x + threadIdx.x;
//    int result_ip = 0;
//    for (size_t i = 0, j = 0; i < MAX_LENGTH_PEPTIDE; ++i, ++j) {
//        int peptide_val = peptides[current_peptide * MAX_LENGTH_PEPTIDE + i];
//        if (peptide_val != -1) {
//            result_ip += spectra_const[peptide_val];
//        }
//    }
//    inner_product[current_peptide] = result_ip;
//}
//
//void comp_ip(int* peptides, int* inner_product, int num_of_peptides) {
//    dim3 blockSize(BLOCK_SIZE, 1, 1);
//    dim3 gridSize(num_of_peptides / BLOCK_SIZE, 1, 1);
//
//    comp_ip_kernel <<<gridSize, blockSize >>> (peptides, inner_product);
//}
//
//void print_vector(std::vector<int> data) {
//    for (int var : data) {
//        std::cout << var << ' ';
//    }
//}
//
//bool is_end_of(std::string& line) {
//    return(line.find('\n') != std::string::npos);
//}
//
//int main()
//{
//    while (cin) {
//        int NUM_OF_PEPTIDES = 0;
//        std::string peptides_file = "";
//        int input = 0;
//        std::cin >> input;
//
//        switch (input) {
//        case 50: {
//            peptides_file = "../../50k_size";
//            NUM_OF_PEPTIDES = 500;
//        }
//               break;
//        case 100: {
//            peptides_file = "../../100k_size";
//            NUM_OF_PEPTIDES = 100000;
//        }
//                break;
//        case 250: {
//            peptides_file = "../../250k_size";
//            NUM_OF_PEPTIDES = 250000;
//        }
//                break;
//        case 500: {
//            peptides_file = "../../500k_size";
//            NUM_OF_PEPTIDES = 500000;
//        }
//                break;
//        case 1000: {
//            peptides_file = "../../1m_size";
//            NUM_OF_PEPTIDES = 1000000;
//        }
//                 break;
//        case 2000: {
//            peptides_file = "../../2m_size";
//            NUM_OF_PEPTIDES = 2000000;
//        }
//                 break;
//        case 4000: {
//            peptides_file = "../../4m_size";
//            NUM_OF_PEPTIDES = 4000000;
//        }
//                 break;
//        }
//
//        int SIZE_PEPTIDES = NUM_OF_PEPTIDES * MAX_LENGTH_PEPTIDE;
//        // Allocate memory on the host
//        vector<int> h_peptides(SIZE_PEPTIDES);
//        //vector<int> h_spectra(SIZE_SPECTRA);
//        int h_spectra[SIZE_SPECTRA];
//        vector<int> h_inner_product(SIZE_INNER_PRODUCT);
//
//        std::fstream spectra("../../spectra");
//        std::string temp_line = "";
//
//        // Fullfill spectra data    
//        for (size_t i = 0; i < 20; ++i) {
//            if (is_end_of(temp_line))
//                break;
//            std::getline(spectra, temp_line, ' ');
//            h_spectra[i] = std::stoi(temp_line);
//        }
//
//        cudaMemcpyToSymbol(spectra_const, h_spectra, SIZE_SPECTRA * sizeof(int));
//
//        // Fullfill peptide data
//        std::fstream peptides(peptides_file);
//        size_t pept_in_seq = 0;
//        size_t pos_n = 0;
//        temp_line = "";
//
//        for (size_t i = 0; i < SIZE_PEPTIDES; ++i) {
//            if (peptides.eof())
//            {
//                temp_line = "-1";
//            }
//            else
//            {
//                std::getline(peptides, temp_line, ' ');
//            }
//            if (is_end_of(temp_line)) {
//                pos_n = temp_line.find('\n'); // position of '\n' character
//                int num1 = std::stoi(temp_line.substr(0, pos_n)); // number before '\n'
//                int num2 = std::stoi(temp_line.substr(pos_n));    // number after '\n'
//               // std::cout << "1" << "\n";
//                h_peptides[i] = num1;   // assign i element
//                //std::cout << "2" << "\n";
//                for (size_t last_peptide = 1; last_peptide < MAX_LENGTH_PEPTIDE - pept_in_seq; ++last_peptide) {
//                    //  std::cout << "3" << "\n";
//                    h_peptides[i + last_peptide] = -1;  // fullfill last elements with 0
//                    //std::cout << "4" << "\n";
//                }
//                i += MAX_LENGTH_PEPTIDE - pept_in_seq; // add last elements to i
//                if (i < SIZE_PEPTIDES) {
//                    //std::cout << "5" << "\n";
//                    h_peptides[i] = num2; // assign element of the next peptide
//                    //std::cout << "6" << "\n";
//                }
//                pept_in_seq = 1;      // increase number of peptides in a sequence
//            }
//            else {
//                //std::cout << "7" << "\n";
//                h_peptides[i] = std::stoi(temp_line);
//                //std::cout << "8" << "\n";
//                ++pept_in_seq;
//            }
//        }
//
//        // Allocate memory on the device
//        dev_array<int> d_peptides(SIZE_PEPTIDES);
//        //dev_array<int> d_spectra(SIZE_SPECTRA);
//        dev_array<int> d_inner_product(SIZE_INNER_PRODUCT);
//
//        // Compute IP
//        auto begin = std::chrono::steady_clock::now();
//
//        d_peptides.set(&h_peptides[0], SIZE_PEPTIDES);
//        //d_spectra.set(&h_spectra[0], SIZE_SPECTRA);
//
//        comp_ip(d_peptides.getData(), d_inner_product.getData(), NUM_OF_PEPTIDES);
//
//        cudaDeviceSynchronize();
//        d_inner_product.get(&h_inner_product[0], SIZE_INNER_PRODUCT);
//        cudaDeviceSynchronize();
//
//        auto end = std::chrono::steady_clock::now();
//        std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;
//
//        ofstream results("results");
//        for (int i = 0; i < SIZE_INNER_PRODUCT; i++) {
//            results << h_inner_product[i] << "\n";
//        }
//    }
//}
