//#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <vector>
//#include <stdlib.h>
//#include <math.h>
//#include <ctime>
//#include <chrono>
//
//using namespace std;
//
//const int ALIGNMENT_ELEM = -1;
//
//int main()
//{
//    auto all_begin = std::chrono::steady_clock::now();
//    std::ofstream time_res("time_results.txt", std::ios::app);
//
//    // Write the time of execution
//    std::time_t t = std::time(0);
//    std::tm* now = std::localtime(&t);
//    time_res
//        << "SEQ: "
//        << (now->tm_year + 1900) << '-'
//        << (now->tm_mon + 1) << '-'
//        << now->tm_mday << ' '
//        << now->tm_hour << ':'
//        << now->tm_min << ':'
//        << now->tm_sec << std::endl;
//
//    auto read_beg = std::chrono::steady_clock::now();
//
//    // Read from files
//    std::string spectra_filename = "../../spectra_10K.txt";
//    std::string peptides_filename = "../../peptides_2M.txt";
//    int peptide_length = 13;
//
//    std::fstream spectra_file(spectra_filename);
//    std::fstream peptides_file(peptides_filename);
//    std::ofstream inner_product_file("cpu_results.txt");
//
//    std::vector<int> h_peptides;
//    std::vector<std::vector<int>> h_spectra;
//    std::vector<int> h_inner_product;
//
//    // Fill spectras
//    size_t max_size_spectra = 0;
//    std::string temp;
//    while (std::getline(spectra_file, temp)) {
//        std::vector<int> cur_spectra;
//        std::stringstream row(temp);
//        int spec_val;
//        while (row >> spec_val) {
//            cur_spectra.push_back(spec_val);
//        }
//        if (cur_spectra.size() > max_size_spectra)
//            max_size_spectra = cur_spectra.size();
//        h_spectra.push_back(cur_spectra);
//    }
//    size_t num_of_spectra = h_spectra.size();
//
//    // Fill peptides
//    size_t num_of_peptides = 0;
//    while (std::getline(peptides_file, temp)) {
//        int pep_size = 0;
//        int pep_val;
//        std::stringstream row(temp);
//        while (row >> pep_val) {
//            h_peptides.push_back(pep_val);
//            ++pep_size;
//        }
//        for (; pep_size < peptide_length; ++pep_size) {
//            h_peptides.push_back(ALIGNMENT_ELEM);
//        }
//        ++num_of_peptides;
//    }
//    auto read_end = std::chrono::steady_clock::now();
//    time_res << "Read: " << std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_beg).count() << std::endl;
//
//    // Count IP
//    auto count_ip_beg = std::chrono::steady_clock::now();
//
//    for (size_t spec_id = 0; spec_id < num_of_spectra; ++spec_id) {
//        std::vector<int> cur_spectra = h_spectra[spec_id];
//        std::vector<int> h_inner_product;
//        size_t pep_id = 0;
//        while (pep_id < num_of_peptides) {
//            int pep_val = 0;
//            int cur_res_ip = 0;
//            while (h_peptides[pep_id * peptide_length + pep_val] != ALIGNMENT_ELEM) {
//                cur_res_ip += cur_spectra[h_peptides[pep_id * peptide_length + pep_val]];
//                ++pep_val;
//            }
//            h_inner_product.push_back(cur_res_ip);
//            ++pep_id;
//        }
//        // operation write results to file
//        //for (size_t i = 0; i < num_of_peptides; ++i) {
//        //    inner_product_file << h_inner_product[i] << '\n';
//        //}
//    }
//
//    auto count_ip_end = std::chrono::steady_clock::now();
//    time_res << "CountIp: " << std::chrono::duration_cast<std::chrono::microseconds>(count_ip_end - count_ip_beg).count() << std::endl;
//
//    // Save results to the file
//    // auto to_res_beg = std::chrono::steady_clock::now();
//    // ofstream inner_product_file("results.txt");
//
//    // auto to_res_end = std::chrono::steady_clock::now();
//    // time_res << "ResFile: " << std::chrono::duration_cast<std::chrono::milliseconds>(to_res_end - to_res_beg).count() << std::endl;
//
//    auto all_end = std::chrono::steady_clock::now();
//    time_res << "All: " << std::chrono::duration_cast<std::chrono::milliseconds>(all_end - all_begin).count() << std::endl
//        << "----------------" << std::endl;
//
//    return 0;
//}
