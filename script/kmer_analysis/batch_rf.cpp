#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "BSlogger/src/BSlogger.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct param {
  std::string kmer_file;
  std::string count_dir;
  std::string meta_file;
  uint64_t batch_size;
  uint64_t features;
};

std::unordered_map<std::string, uint64_t> get_metadata(std::string const & path)
{
  logger log(std::cerr, "get_metadata");
  std::unordered_map<std::string, uint64_t> ret;
  std::ifstream ifs(path.c_str(), std::ios::in);
  std::string line;
  while (std::getline(ifs, line))
  {
    std::stringstream ss(line);
    std::string key, value;
    ss >> key;
    ss >> value;
    if (value == "F")
    {
      ret[key] = 0;
    }
    else if (value == "M")
    {
      ret[key] = 1;
    }
    else
    {
      log(LOG_WARN) << "Skipping key " << key << " with non M/F value: " << value << '\n';
    }
  }
  return ret;
}

std::vector<double> get_normfactors(std::unordered_map<std::string, uint64_t> const & map,
    std::string const & path)
{
  std::vector<double> ret;
  for (auto const & entry : map)
  {
    fs::path p(path);
    fs::path suf(entry.first + ".cov");
    p /= suf;
    std::ifstream ifs(p.string().c_str(), std::ios::in);
    std::string l;
    std::getline(ifs, l);
    ret.push_back(std::stod(l));
  }
  return ret;
}

std::vector<std::ifstream> get_filehandles(std::unordered_map<std::string, uint64_t> const & map,
    std::string const & path)
{
  std::vector<std::ifstream> ret;
  for (auto const & entry : map)
  {
    fs::path p(path);
    fs::path suf(entry.first + ".count.lcount");
    p /= suf;
    ret.push_back(std::ifstream(p.string().c_str(), std::ios::in));
  }
  return ret;
}

std::vector<uint64_t> get_classes(std::unordered_map<std::string, uint64_t> const & map)
{
  std::vector<uint64_t> ret;
  for (auto const & entry : map)
  {
    ret.push_back(entry.second);
  }
  return ret;
}

int main(int argc, char const ** argv)
{
  logger log(std::cerr, "Batch_RF");
  param opts;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show this help")
    ("kmers,k", po::value<std::string>(&opts.kmer_file)->required(),
     "File with kmers")
    ("meta,m", po::value<std::string>(&opts.meta_file)->required(),
     "File with metadata")
    ("dir,d", po::value<std::string>(&opts.count_dir)->required(),
     "Directory with counts")
    ("features,n", po::value<uint64_t>(&opts.features)->required(),
     "Number of features")
    ("batch_size,b", po::value<uint64_t>(&opts.batch_size)->default_value(10000),
     "Number of kmers to fit at once")
    ;

  po::variables_map vm;

  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  if (vm.count("help") > 0)
  {
    std::cerr << desc << '\n';
    return 0;
  }

  try
  {
    po::notify(vm);
  }
  catch(std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  try
  {
    opts.meta_file = fs::canonical(opts.meta_file).string();
    opts.kmer_file = fs::canonical(opts.kmer_file).string();
    opts.count_dir = fs::canonical(opts.count_dir).string();
  }
  catch(std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  auto meta = get_metadata(opts.meta_file);
  auto normfactors = get_normfactors(meta, opts.count_dir);
  auto file_handles = get_filehandles(meta, opts.count_dir);
  auto file_classes = get_classes(meta);

  std::ifstream kmer_handle(opts.kmer_file.c_str(), std::ios::in);
  progbar_fancy<uint64_t> progbar(std::cerr, opts.features, 1000, 30, "kmer");
  while (true)
  {
    uint64_t rec_read = 0;
    fs::path tempfile = fs::unique_path();
    std::ofstream ofs(tempfile.string().c_str(), std::ios::out);
    std::string line;
    while (std::getline(kmer_handle, line) && (rec_read < opts.batch_size))
    {
      ofs << line << ' ';
      rec_read++;
    }
    ofs << "Sex\n";
    if (rec_read == 0)
    {
      ofs.close();
      fs::remove(tempfile);
      break;
    }
    for (uint64_t i = 0; i < file_handles.size(); i++)
    {
      rec_read = 0;
      while (std::getline(file_handles[i], line) && (rec_read < opts.batch_size))
      {
        double d = std::stod(line);
        d /= (normfactors[i] / 1e9);
        ofs << d << ' ';
        rec_read++;
      }
      ofs << file_classes[i] << '\n';
    }
    ofs << std::flush;
    ofs.close();
    std::stringstream cmd;
    cmd << "ranger"
        << " --file " << tempfile.string()
        << " --depvarname Sex"
        << " --impmeasure 5"
        << " --ntree 10000"
        << " --outprefix " << tempfile.string();
    //log(LOG_INFO) << cmd.str() << '\n';
    if (system(cmd.str().c_str()) > 0)
    {
      break;
    }
    progbar += rec_read;
    fs::remove(tempfile);
    //break;
  }
  progbar.finalize();

  return 0;
}
