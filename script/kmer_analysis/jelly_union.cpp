#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <chrono>
#include <cmath>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/cpp_array.hpp>

#include "robin_hood.h"

#define XXH_INLINE_ALL
#include "xxhash.h"
//typedef uint64_t XXH64_hash_t;
#include "BSlogger/src/BSlogger.hpp"

/*
 Why are these values hard coded instead of program options?
 Knowing these values at compile time allows for some optimizations that would not
 be possible in a dynamic program, i.e. it will be faster this way.
*/
#define BUFSIZE 1000000000
#define K 31
#define GENOME 480000000

class kmer {
public:
  kmer() {  }
  kmer(std::string const && str) { data = str; }
  bool operator==(kmer const & other) const {
    return data == other.data;
  }
  std::string data;
};

class fasthash_k {
public:
  XXH64_hash_t operator()(kmer const & km) const
  {
    return XXH64(&km.data[0], K, 31415);
  }
};

void add_to_set(std::string const &  f, // Jellyfish file
                robin_hood::unordered_map<kmer, XXH64_hash_t, fasthash_k>& set, // map target
                uint64_t& index, // enumeration of files to read
                uint64_t& cov, // total coverage of the fastq file in BP
                std::vector<jellyfish::mer_dna>& m_vec, // vector with kmer sequences (ram preallocated)
                std::vector<uint64_t>& c_vec) // vector with kmer counts (ram preallocated)
{
  logger log(std::cerr, "add_to_set");
  std::ifstream in_stream(f.c_str());
  jellyfish::file_header header(in_stream);
  jellyfish::mer_dna::k(header.key_len() / 2);
  binary_reader reader(in_stream, &header);
  uint64_t lo_pass = (((cov / GENOME) + 1) / 4);
  lo_pass = lo_pass < 1 ? 1 : lo_pass;
  log(LOG_INFO) << "Low pass filter: " << lo_pass << '\n';

  // Faster to first buffer all kmers and then start hashing all at once
  log(LOG_INFO) << "Reading into memory... " << header.key_len() << '\n';
  while (reader.next())
  {
    m_vec.push_back(reader.key());
    c_vec.push_back(reader.val());
  }

  log(LOG_INFO) << "Read data, now hashing at most " << m_vec.size() << " elements\n";
  progbar_fancy<uint64_t> progbar(std::cerr, m_vec.size(), 1000, 80, "Kmer");
  for (uint64_t i = 0; i < m_vec.size(); i++)
  {
    uint64_t const & val = c_vec[i];
    progbar++;
    if (val > lo_pass)
    {
      set[m_vec[i].to_str()]++; // on ket miss, default initialized to 0
    }
  }
  progbar.finalize();
}

uint64_t read_cov(std::string const & x)
{
  std::ifstream ifs(x.c_str());
  uint64_t ret;
  ifs >> ret;
  return ret;
}

int main(int argc, char ** argv)
{
  logger log(std::cerr, "jelly_union");
  robin_hood::unordered_map<kmer, XXH64_hash_t, fasthash_k> set;
  log(LOG_INFO) << "Reserving: " << BUFSIZE << " elements in HashMap\n";
  set.reserve(BUFSIZE);

  log(LOG_INFO) << "Allocating space for I/O buffers...\n";
  std::vector<jellyfish::mer_dna> m_vec;
  std::vector<uint64_t> c_vec;
  m_vec.reserve(BUFSIZE);
  c_vec.reserve(BUFSIZE);

  uint64_t index = 0;
  double secs = 0;
  for (int i = 1; i < argc; i++)
  {
    m_vec.clear();
    c_vec.clear();
    auto start = std::chrono::high_resolution_clock::now();
    log(LOG_INFO) << "Adding: " << argv[i] << '\n';
    std::string cov_path(argv[i]);
    cov_path += ".cov";
    uint64_t cov = read_cov(cov_path);
    log(LOG_INFO) << "Number of bases in raw input: " << cov << '\n';
    add_to_set(std::string(argv[i]), set, index, cov, m_vec, c_vec);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    log(LOG_INFO) << "Took " << duration.count() << " seconds\n";
    log(LOG_INFO) << "Set size is " << set.size() << '\n';
    secs += duration.count();
    double eta = secs * static_cast<double>((argc - i));
    std::chrono::duration<uint64_t> dur =
      std::chrono::duration_cast<std::chrono::duration<uint64_t>>(stop - start);
    log(LOG_INFO) << "ETA: "
                  << format_duration<uint64_t>(dur.count()) << '\n';
    index++;
  }

  double nsamples = static_cast<double>(argc);
  uint64_t cutoff = std::ceil(nsamples / 5);

  for (auto const & entry : set)
  {
    if (entry.second >= cutoff)
    {
      std::cout << entry.first.data << '\n';
    }
  }
  return 0;
}
