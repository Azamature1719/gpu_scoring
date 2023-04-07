#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <chrono>

const size_t NUM_OF_PEPTIDES = 2000000;  // number of peptides in document 
const size_t MAX_LENGTH_PEPTIDE = 10;
const size_t SIZE_PEPTIDES = NUM_OF_PEPTIDES * MAX_LENGTH_PEPTIDE;
const size_t SIZE_SPECTRA = 20;
const size_t SIZE_INNER_PRODUCT = NUM_OF_PEPTIDES;  // ip is counted for each peptide


int main()
{
    auto all_begin = std::chrono::steady_clock::now();
    std::ofstream time_res("time_results.txt", std::ios::app);

    // Write the time of execution
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
    time_res
        << "SEQ: "
        << (now->tm_year + 1900) << '-'
        << (now->tm_mon + 1) << '-'
        << now->tm_mday << ' '
        << now->tm_hour << ':'
        << now->tm_min << ':'
        << now->tm_sec << std::endl;

    auto read_beg = std::chrono::steady_clock::now();

    // Read values from files
    std::fstream file_spectra("../../spectra.txt");
    std::fstream file_peptides("../../2m_pep.txt");

    // Allocate memory on the host
    int* h_peptides = new int[SIZE_PEPTIDES];
    int* h_spectra = new int[SIZE_SPECTRA];
    int* h_inner_product = new int[SIZE_INNER_PRODUCT];

    for (int i = 0; i < SIZE_SPECTRA; ++i) {
        int spec_val;
        file_spectra >> spec_val;
        h_spectra[i] = spec_val;
    }

    for (int i = 0; i < NUM_OF_PEPTIDES; i++)
    {
        for (int j = 0; j < MAX_LENGTH_PEPTIDE; ++j) {
            int pep_val;
            file_peptides >> pep_val;
            h_peptides[i * MAX_LENGTH_PEPTIDE + j] = pep_val;
        }
    }
    auto read_end = std::chrono::steady_clock::now();
    time_res << "Read: " << std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_beg).count() << std::endl;

    // Count IP
    auto count_ip_beg = std::chrono::steady_clock::now();
    
    for (size_t i = 0; i < NUM_OF_PEPTIDES; ++i) {
        for (size_t j = 0; j < MAX_LENGTH_PEPTIDE; ++j) {
            int peptide_val = h_peptides[MAX_LENGTH_PEPTIDE * i + j];
            h_inner_product[i] += h_spectra[peptide_val];
        }
    }

    auto count_ip_end = std::chrono::steady_clock::now();
    time_res << "CountIp: " << std::chrono::duration_cast<std::chrono::milliseconds>(count_ip_end - count_ip_beg).count() << std::endl;

    // // Save results to the file
    // auto to_res_beg = std::chrono::steady_clock::now();
    // std::ofstream inner_product_file("results.txt");
    // for (int i = 0; i < SIZE_INNER_PRODUCT; i++)
    //     inner_product_file << h_inner_product[i] << std::endl;
    // auto to_res_end = std::chrono::steady_clock::now();
    // time_res << "ResFile: " << std::chrono::duration_cast<std::chrono::milliseconds>(to_res_end - to_res_beg).count() << std::endl;

    auto all_end = std::chrono::steady_clock::now();
    time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
        << "----------------" << std::endl;

    return 0;
}
