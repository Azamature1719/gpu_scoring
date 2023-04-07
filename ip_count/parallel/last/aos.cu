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

const int NUM_OF_PEPTIDES = 2000000;
const int LENGTH_PEPTIDE = 10;
struct peptide {
    int A;
    int B;
    int C;
    int D;
    int E;
    int F;
    int G;
    int H;
    int I;
    int J;
};

__global__ void kernel(const int* spectra, peptide* peps, int* inner_product)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    inner_product[id] = spectra[peps[id].A] + spectra[peps[id].B] + spectra[peps[id].C] + spectra[peps[id].D] + spectra[peps[id].E] + spectra[peps[id].F] + spectra[peps[id].G] + spectra[peps[id].H] + spectra[peps[id].I] + spectra[peps[id].J];
}

void compute_inner_product(const std::vector<int>& spectra, peptide* peps, std::vector<int>& inner_product)
{
    int spectra_size = spectra.size();
    inner_product.resize(NUM_OF_PEPTIDES);
    int* d_spectra, * d_inner_product;
    peptide* d_peptides = new peptide[NUM_OF_PEPTIDES];

    auto begin = std::chrono::steady_clock::now();

    cudaMalloc(&d_spectra, spectra_size * sizeof(int));
    cudaMalloc(&d_peptides, NUM_OF_PEPTIDES * sizeof(peptide));
    cudaMalloc(&d_inner_product, NUM_OF_PEPTIDES * sizeof(int));

    cudaMemcpy(d_spectra, spectra.data(), spectra_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaError_t err = cudaMemcpy(d_peptides, peps, NUM_OF_PEPTIDES * sizeof(peptide), cudaMemcpyHostToDevice);

    int block_size = 1024;
    int grid_size = NUM_OF_PEPTIDES / block_size;

    kernel <<<grid_size, block_size>>> (d_spectra, d_peptides, d_inner_product);
    cudaMemcpy(inner_product.data(), d_inner_product, NUM_OF_PEPTIDES * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(d_spectra);
    cudaFree(d_peptides);
    cudaFree(d_inner_product);

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed for GPU: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;

}

bool is_end_of(std::string& line) {
    return(line.find('\n') != std::string::npos);
}

// -- Test on 2M hardcode --
int main(int argc, char* argv[])
{
    auto begin = std::chrono::steady_clock::now();

    std::vector<int> spectra;
    std::fstream spectra_file("../../spectra");
    std::fstream peptides_file("../../2m_size");
    std::ofstream inner_product_file("results_aos");

    std::string temp_line = "";
    for (size_t i = 0; i < 20; ++i) {
        if (is_end_of(temp_line))
            break;
        std::getline(spectra_file, temp_line, ' ');
        spectra.push_back(std::stoi(temp_line));
    }

    peptide *peps = new peptide[NUM_OF_PEPTIDES];
    for (int i = 0; i < NUM_OF_PEPTIDES; ++i) {
        peptides_file >> peps[i].A;
        peptides_file >> peps[i].B;
        peptides_file >> peps[i].C;
        peptides_file >> peps[i].D;
        peptides_file >> peps[i].E;
        peptides_file >> peps[i].F;
        peptides_file >> peps[i].G;
        peptides_file >> peps[i].H;
        peptides_file >> peps[i].I;
        peptides_file >> peps[i].J;
    }

    //for (int i = 0; i < NUM_OF_PEPTIDES; ++i) {
    //        std::cout << peps[i].A << " " <<
    //        peps[i].B << " " <<
    //        peps[i].C << " " <<
    //        peps[i].D << " " <<
    //        peps[i].E << " " <<
    //        peps[i].F << " " <<
    //        peps[i].G << " " <<
    //        peps[i].H << " " <<
    //        peps[i].I << " " <<
    //        peps[i].J;
    //}

    std::vector<int> inner_product;
    compute_inner_product(spectra, peps, inner_product);

    for (int i = 0; i < NUM_OF_PEPTIDES; i++)
    {
        inner_product_file << inner_product[i] << std::endl;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed for the program: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;

    return 0;
}