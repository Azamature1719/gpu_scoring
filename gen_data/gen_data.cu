#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

#define NUM_OF_PEPTIDES_IN_FILE 5

using namespace std;

void generate_peptides(size_t num_of_peptides, std::string filename){
    std::ofstream big_peptides(filename);
    std::string pep = "";
    num_of_peptides /= NUM_OF_PEPTIDES_IN_FILE;
    for(size_t i = 0; i < num_of_peptides; ++i){
        std::fstream peptides("spectra.txt");
        pep = "";
        while(std::getline(peptides, pep)){
            big_peptides << pep;
            big_peptides << '\n';
        }
    }
}

int main()
{
   //generate_peptides(10000, "peptides_10K.txt");
   //generate_peptides(100000, "100k_size.txt");
   //generate_peptides(250000, "250k_size.txt");
   //generate_peptides(500000, "500k_size.txt");
   //generate_peptides(1000000, "1m_size.txt");
   generate_peptides(10000, "10k_spectra.txt");
   //generate_peptides(4000000, "4m_size.txt");
   //generate_peptides(5000000, "5m_size.txt");
//    generate_peptides(10000000, "10m_size.txt");
//    generate_peptides(25000000, "25m_size.txt");
//    generate_peptides(50000000, "50m_size.txt");
//   generate_peptides(1000000000, "1b_size.txt");
}
