#include "kernel.cuh"

using namespace std;

 Global memory
__global__ void compute_ip_kernel(int* spectra, int* peptides, int* inner_product) {

    int current_peptide = blockDim.x * blockIdx.x + threadIdx.x;
    int temp_inner_product = 0;
    for (size_t i = 0, j = 0; i < MAX_LENGTH_PEPTIDE; ++i, ++j) {
        int peptide_val = peptides[current_peptide * MAX_LENGTH_PEPTIDE + i];
        if (peptide_val != -1) {
            temp_inner_product += spectra[peptide_val];
        }
    }
    inner_product[current_peptide] = temp_inner_product;
}

__global__ void find_ip_kernel(int* spectra, int* peptides, int* inner_product) {

    int current_peptide = blockDim.x * blockIdx.x + threadIdx.x;
    int temp_inner_product = 0;
    for (size_t i = 0, j = 0; i < MAX_LENGTH_PEPTIDE; ++i, ++j) {
        int peptide_val = peptides[current_peptide * MAX_LENGTH_PEPTIDE + i];
        if (peptide_val != -1) {
            temp_inner_product += spectra[peptide_val];
        }
    }
    inner_product[current_peptide] = temp_inner_product;
}

void find_ip(int* spectra, int* peptides, int* inner_product) {
    // Threads per block == Num of peptides; Num of blocks = 1 in a grid

    dim3 blockSize(BLOCK_SIZE, 1, 1);
    dim3 gridSize(NUM_OF_PEPTIDES / BLOCK_SIZE, 1, 1);

    find_ip_kernel <<<gridSize, blockSize >>> (spectra, peptides, inner_product);
}