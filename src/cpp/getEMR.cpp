#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>

#define INDEL 0
#define SNP 1
#define MAX_DISTANCE 250

typedef struct {
  size_t position;
  unsigned short length;
  char type;
} variant;

bool compare_variants(variant a, variant b)
{
  return a.position < b.position;
}

int main(int argc, char ** argv)
{

  std::string in_snp = argv[1];
  std::string in_indel = argv[2];
  std::map<std::string, std::vector<variant> > variants;

  std::ifstream ifs;
  ifs.open(in_snp.c_str(), std::ios::in);
  std::string line;
  bool first = true;
  while(std::getline(ifs, line))
    {
      if(first)
	{
	  first = false; continue;
	}
      std::stringstream ss(line);
      std::string field;
      unsigned short col = 0;

      std::string chrom;
      size_t pos;
      std::string seq;

      while(std::getline(ss, field, '\t'))
	{
	  if(col == 0) chrom = field;
	  if(col == 1) pos = std::stoul(field);
	  if(col == 2) seq = field;
	  col++;
	}
      size_t p1 = seq.find("(");
      size_t p2 = seq.find("/");
      size_t p3 = seq.find(")");
      
      std::string ref = seq.substr(p1 + 1, p2 - p1 - 1);
      std::string alt = seq.substr(p2 + 1, p3 - p2 - 1);
      
      variant x;
      x.position = pos;
      x.length = ref.length() > alt.length() ? ref.length() : alt.length();
      x.type = INDEL;

      variants[chrom].push_back(x);
      
    }

  ifs.close();

  ifs.open(in_indel.c_str(), std::ios::in);
  first = true;
  while(std::getline(ifs, line))
    {
      if(first)
        {
          first = false; continue;
        }
      std::stringstream ss(line);
      std::string field;
      unsigned short col = 0;

      std::string chrom;
      size_t pos;
      std::string seq;

      while(std::getline(ss, field, '\t'))
        {
          if(col == 0) chrom = field;
          if(col == 1) pos = std::stoul(field);
          if(col == 2) seq = field;
          col++;
        }
      size_t p1 = seq.find("(");
      size_t p2 = seq.find("/");
      size_t p3 = seq.find(")");

      std::string ref = seq.substr(p1 + 1, p2 - p1 - 1);
      std::string alt = seq.substr(p2 + 1, p3 - p2 - 1);

      variant x;
      x.position = pos;
      x.length = ref.length() > alt.length() ? ref.length() : alt.length();
      x.type = INDEL;

      variants[chrom].push_back(x);

    }

  ifs.close();

  for(auto it = variants.begin(); it != variants.end(); it++)
    {
      std::vector<size_t> starts;
      std::vector<size_t> ends;
      std::vector<size_t> errors;
      std::vector<size_t> lengths;
      std::sort(it->second.begin(), it->second.end(), compare_variants);
      size_t sum_errors = 0;
      std::vector<variant> current = it->second;
      
      auto start = current.begin();
      for(auto it2 = current.begin(); it2 != current.end(); it2++)
	{
	
	  if(it2->position - start->position < MAX_DISTANCE)
	    {
	      sum_errors += it2->length;
	    }
	  else
	    {
	      --it2;
	      if(sum_errors == 0)
		sum_errors += it2->length;
	      starts.push_back(start->position);
	      ends.push_back(it2->position);
	      errors.push_back(sum_errors);
	      lengths.push_back(it2->position - start->position);
	      start = ++it2;
	      sum_errors = 0;
	    }
	  
	}

      for(size_t i = 0; i < starts.size(); i++)
	{
	  std::cout << it->first << '\t' << starts[i] << '\t' << ends[i] << 
	    '\t' << lengths[i]+1 << '\t'  << errors[i] << std::endl;
	}

    }



  return 0;
}
