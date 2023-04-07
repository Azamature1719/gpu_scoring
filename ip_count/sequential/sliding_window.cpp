#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <chrono>
#include <limits>
#include <algorithm>

const int ALIGNMENT_ELEM = -1;

const int peptide_length = 13;
const int spectra_length = 20;
const int delta = 10000; // Difference between starting points
const int match_range = 200000; // Number of peptides that should be matched with a single spectra

// -- To develop a custom configurable program -- 
const int number_of_launches = 56; // Calculated manually
const int spectra_at_once = 180; // Calculated manually

int main(int argc, char* argv[]) {
    auto all_begin = std::chrono::steady_clock::now();
    std::ofstream time_res("time_results.txt", std::ios::app);
    std::ofstream res("spec_results.txt", std::ios::app);

    // Write the time of execution
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
    time_res
        << "SEQ 2: "
        << (now->tm_year + 1900) << '-'
	<< now->tm_mon << "-"
        << now->tm_mday << ' '
        << now->tm_hour << ':'
        << now->tm_min << ':'
        << now->tm_sec << std::endl;

    // Get vars from command line
    //std::string spectra_filename = "../../" + std::string(argv[1]);
    //std::string peptides_filename = "../../" + std::string(argv[2]);
    //int peptide_length = std::atoi(argv[3]);
    //int block_size = std::atoi(argv[4]);

    std::string spectra_filename = "../10k_spectra.txt";
    std::string peptides_filename = "../2m_pep.txt";

    std::fstream spectra_file(spectra_filename);
    std::fstream peptides_file(peptides_filename);

    std::vector<int> h_peptides;
    std::vector<int> h_spectra;
    std::string temp;

    auto read_begin = std::chrono::steady_clock::now();

    // Fill spectras
    int num_of_spectra = 0;
    while (std::getline(spectra_file, temp)) {
        int spec_val;
        std::stringstream row(temp);
        while (row >> spec_val) {
            h_spectra.push_back(spec_val);
        }
        ++num_of_spectra;
    }

    /*/ Fill peptides
    size_t num_of_peptides = 0;
    size_t j = 0;
    while (std::getline(peptides_file, temp)) {
        int pep_size = 0;
        int pep_val;
        std::stringstream row(temp);
        while (row >> pep_val) {
            h_peptides[j] = pep_val;
            ++pep_size;
	    ++j;
        }
        for (; pep_size < peptide_length; ++pep_size) {
            h_peptides[j] = ALIGNMENT_ELEM;
	    ++j;
        }
        ++num_of_peptides;
    }*/

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

    std::vector<int> inner_product(num_of_spectra * match_range);

    auto ann_begin = std::chrono::steady_clock::now();
    int match_pep_id = 0;
     for(int j = 0; j < 56; ++j){
	    for (int spec_id = 0; spec_id < 180; ++spec_id) {
	        for (match_pep_id = 0; match_pep_id < match_range; ++match_pep_id) { // Num of the peptide in a 200K peptide sequence
	            int match_pep_val_id = 0; // Order num of the peptide value in the concrete peptide 
	            int cur_pep_val = spec_id * delta * peptide_length + match_pep_id * peptide_length + match_pep_val_id; // Order num of the peptide value in the peptides array. t_id is used because peptides are updated each launch
	            int ip_val_id = spec_id * match_range + match_pep_id;
	            inner_product[ip_val_id] = 0;
	            while (h_peptides[cur_pep_val] != ALIGNMENT_ELEM) {
	                inner_product[ip_val_id] += h_spectra[spec_id * spectra_length + h_peptides[cur_pep_val]]; // spec_id * match_range - each spectra is matched with 200K peptides
	                ++match_pep_val_id;
	                cur_pep_val = spec_id * delta * peptide_length + match_pep_id * peptide_length + match_pep_val_id;
	            }
	//	   std::cout << "cur_pep_val " << cur_pep_val << " ";
	        }
	       //std::cout << "spec_id: " <<spec_id << " " << "pep_id" << match_pep_id;
	    }
	
     }
    
    auto ann_end = std::chrono::steady_clock::now();
   
    time_res << "Annotation: " << std::chrono::duration_cast<std::chrono::milliseconds>(ann_end - ann_begin).count() << std::endl;


/*    for(int i = 0; i < num_of_spectra * match_range; ++i){
	res << inner_product[i] << "\n";
    }*/

    auto all_end = std::chrono::steady_clock::now();
    time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
        << "----------------" << std::endl;

    return 0;
}

