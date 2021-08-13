/*
  Sparse kmer count matrix filter
*/

#include <iostream>
#include <string>
#include <fstream>
#include <armadillo>
#include <vector>
#include <omp.h>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/mapped_file.hpp>
#include <unordered_map>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "BSlogger/src/BSlogger.hpp"

#define K 31

uint64_t get_nkmers(char const * path)
{
  struct stat buffer;
  stat(path, &buffer);
  return buffer.st_size / (K + 1);
}

int main(int argc, char const ** argv)
{

  logger log(std::cerr, "sparse_vec");
  std::ifstream ifs(argv[1]);
  std::string kmer;
  jellyfish::mer_dna::k(K);

  std::vector<std::string> jellyfish_files;
  uint64_t nk = get_nkmers(argv[1]);

  std::unordered_map<std::string, uint64_t> k_count;

  log(LOG_INFO) << "Reading kmers from union\n";
  uint64_t km_cnt = 0;
  std::vector<jellyfish::mer_dna> km_vec;
  km_vec.reserve(nk);

  progbar_fancy<uint64_t> pb(std::cerr, nk);
  while (std::getline(ifs, kmer))
  {
    km_vec.push_back(jellyfish::mer_dna(kmer));
    ++km_cnt;
    ++pb;
  }
  pb.finalize();

  for (int i = 2; i < argc; i++)
  {
    jellyfish_files.push_back(std::string(argv[i]));
  }

  uint64_t nz_k = 0;
  log(LOG_INFO) << "Adding files...\n";
  #pragma omp parallel for
  for (uint64_t fi = 0; fi < jellyfish_files.size(); fi++)
  {
    auto const & file = jellyfish_files[fi];
    #pragma omp critical
    {
      log(LOG_INFO) << "Thread " << omp_get_thread_num() << " - File: "
                    << jellyfish_files[fi]
                    << " (" << fi << "/" << jellyfish_files.size() << ")\n";
    }
    std::ifstream ifj(file.c_str(), std::ios::in | std::ios::binary);
    jellyfish::file_header header(ifj);
    jellyfish::mapped_file f(file.c_str());
    f.load();
    binary_query bq(f.base() + header.offset(),
                    header.key_len(), header.counter_len(), header.matrix(),
                    header.size() - 1, f.length() - header.offset());

    std::vector<uint64_t> ybuffer;
    ybuffer.resize(nk);

    for (uint64_t i = 0; i < nk; i++)
    {
      uint64_t val = bq.check(km_vec[i]);
      ybuffer[i] = val;
    }
    #pragma omp critical
    {
      log(LOG_INFO) << "Thread " << omp_get_thread_num() << " - File: "
                    << jellyfish_files[fi]
                    << " flushing to file\n";
    }
    std::ofstream ofs((jellyfish_files[fi] + ".count").c_str());
    for (uint64_t i = 0; i < nk; i++)
    {
      ofs << ybuffer[i] << (i == (nk - 1) ? '\n' : '\t');
    }
  } // end omp parallel for
  return 0;
}
