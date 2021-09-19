#include <sam.h>
#include <bgzf.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <map>
#include <tclap/CmdLine.h>

#define TOLERANCE 250 //Tolerance around a region to look for het SNPs

typedef unsigned short ushort;
typedef unsigned long ulong;

static std::map<std::string, unsigned long> order;


struct samtools_sort : public std::binary_function <bool, std::string, std::string> {   
  bool operator() (const std::string& a, const std::string& b)   {
    return order[a] < order[b];
  } 
};

typedef std::map<std::string, std::vector<std::pair<ulong, ulong> >, samtools_sort> crange;
typedef std::map<std::string, std::vector<std::pair<ulong, char> >, samtools_sort> refmap;

class region {
public:
  region(void){start = 0; stop = 0;}
  region(std::string c, ulong sa, ulong so){chrom = c; start = sa; stop = so;}
  friend bool operator== (const region& a, const region& b);
  friend bool operator!= (const region& a, const region& b);
  friend bool operator< (const region& a, const region& b);
  friend bool operator> (const region& a, const region& b);
  std::string chrom;
  ulong start;
  ulong stop;
};

// Class to hold homozygous records
class hom_record {
public:
  void set_chrom(std::string x){chrom = x;}
  void set_start(ulong x){start = x;}
  void set_end(ulong x){end = x;}
  void print(std::ostream& o){o << chrom << '\t' << start << '\t' << end << std::endl;}
  hom_record(void){start = 0; end = 0;}
  friend bool operator== (const hom_record& a, const hom_record& b);
  friend bool operator!= (const hom_record& a, const hom_record& b);
  friend bool operator< (const hom_record& a, const hom_record& b);
  friend bool operator> (const hom_record& a, const hom_record& b);
  std::string chrom;
  ulong  start;
  ulong  end;
};

bool operator== (const region& a, const region& b){
  if(a.chrom == b.chrom && a.start == b.start && a.stop == b.stop)
    return true;

  return false;
}

bool operator!= (const region& a, const region& b){
  if( a == b)
    return false;

  return true;
}

bool operator< (const region& a, const region& b){
  if(order[a.chrom] < order[b.chrom])
    return true;
  else if(order[a.chrom] > order[b.chrom])
    return false;
  else if(a.start < b.start)
    return true;
  else if(a.start > b.start)
    return false;
  else if(a.stop < b.stop)
    return true;
  
  return false;
}

bool operator> (const region& a, const region& b){
  if(a < b)
    return false;
  return true;
}


// Comparison operators to sort hom_records
bool operator== (const hom_record& a, const hom_record& b){
  if(a.chrom == b.chrom && a.start == b.start && a.end == b.end)
    return true;

  return false;
}

bool operator!= (const hom_record& a, const hom_record& b){
  if( a == b)
    return false;

  return true;
}

bool operator< (const hom_record& a, const hom_record& b){
  if(order[a.chrom] < order[b.chrom])
    return true;
  else if(order[a.chrom] > order[b.chrom])
    return false;
  else if(a.start < b.start)
    return true;
  else if(a.start > b.start)
    return false;
  else if(a.end < b.end)
    return true;
  
  return false;
}

bool operator> (const hom_record& a, const hom_record& b){
  if(a < b)
    return false;
  return true;
}

// Binary search to check if a SNP is within a region
bool in(crange& map, std::string& chr, ulong& p){
  if(map.find(chr) == map.end())
    return false;

  std::vector<std::pair<ulong, ulong> > range = map[chr];
  long low = 0;
  long high = range.size() - 1;

  while(true)
    {
      if(high < low)
	return false;
      long mid = low + (high - low) / 2;
      if(p > range[mid].second + TOLERANCE)
	{
	  low = mid + 1; continue;
	}
      if(p < range[mid].first - TOLERANCE)
	{
	  high = mid - 1; continue;
	}
      if(p <= range[mid].second + TOLERANCE && p >= range[mid].first - TOLERANCE)
	return true;
    }
}

// Class to hold heterozygous SNP record
class het_record {
public:
  het_record(){pos=0;ref=0;alt=0;}
  void print(std::ostream& o) {o << chrom << '\t' << pos << 
      '\t' << ref << '/' << alt << std::endl;}
  friend bool operator== (const het_record& a, const het_record& b);
  friend bool operator!= (const het_record& a, const het_record& b);
  friend bool operator< (const het_record& a, const het_record& b);
  friend bool operator> (const het_record& a, const het_record& b);
  std::string chrom;
  ulong pos;
  char ref;
  char alt;
};

bool operator== (const het_record& a, const het_record& b){
  if(a.chrom == b.chrom && a.pos == b.pos)
    return true;

  return false;
}

bool operator!= (const het_record& a, const het_record& b){
  if( a == b)
    return false;

  return true;
}

bool operator< (const het_record& a, const het_record& b){
  if(order[a.chrom] < order[b.chrom])
    return true;
  else if(order[a.chrom] > order[b.chrom])
    return false;
  else if(a.pos < b.pos)
    return true;
  else if(a.pos > b.pos)
    return false;
  
  return false;
}

bool operator> (const het_record& a, const het_record& b){
  if(a < b)
    return false;
  return true;
}


char bam2nuc(int enc)
{
  if(enc & 1) return 'A';
  if(enc & 2) return 'C';
  if(enc & 4) return 'G';
  if(enc & 8) return 'T';
  return 'N';
}

int main(int argc, char ** argv)
{

  std::string in_bam;
  std::string in_hets;
  std::string in_homs;
  std::set<region> regions;
  std::string in_order;

  try 
    {
      TCLAP::CmdLine cmd("Extract alignments matching ref and alt genotypes", ' ', "0.1");

      TCLAP::ValueArg<std::string> arg_bam("b", "bam", "BAM file of alignments.", true, "", "string");
      TCLAP::ValueArg<std::string> arg_hets("e", "het", "File containing heterozygous SNPs.", true, "", "string");
      TCLAP::ValueArg<std::string> arg_homs("o", "hom", "File containing the clusters of homozygous SNPs.", true, "", "string");
      TCLAP::ValueArg<std::string> arg_order("s", "order", "File containing the sort order from the BAM file.", true, "", "string");

      cmd.add(arg_bam);
      cmd.add(arg_hets);
      cmd.add(arg_homs);
      cmd.add(arg_order);

      cmd.parse(argc, argv);

      in_bam = arg_bam.getValue();
      in_hets = arg_hets.getValue();
      in_homs = arg_homs.getValue();
      in_order = arg_order.getValue();
    } 
  catch(TCLAP::ArgException &e) 
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

      

  crange region_map;
  std::ifstream ifs_order(in_order.c_str());
  std::string line;
  unsigned long l = 0; 
  while(std::getline(ifs_order, line)){
    order[line] = l++;
  }

  // Create map of interesting regions from homo snps
  crange hom_map;
  std::ifstream ifs_homs(in_homs.c_str());

  while(std::getline(ifs_homs, line))
    {
      std::string col;
      std::stringstream ss(line);
      //unsigned short fn = 0;
      hom_record r;
      for(ushort i = 0; i < 5; i++)
	{
	  std::getline(ss, col, '\t');
	  if(i == 0) r.set_chrom(col);
	  if(i == 1) r.set_start(std::stoul(col));
	  if(i == 2) r.set_end(std::stoul(col));
	}
      
      if(static_cast<long>(r.start) - TOLERANCE < 1)
	{
	  regions.emplace(region(r.chrom, 1, r.end + TOLERANCE));
	}
      else
	{
	  regions.emplace(region(r.chrom, r.start - TOLERANCE, r.end + TOLERANCE));
	}

      hom_map[r.chrom].push_back(std::pair<ulong, ulong>(r.start, r.end));

    }

  for(auto it = regions.begin(); it != regions.end(); it++)
    {
      region_map[it->chrom].push_back(std::pair<ulong, ulong>(it->start, it->stop));
    }
  
  for(auto it = hom_map.begin(); it != hom_map.end(); it++)
    {
      std::sort(it->second.begin(), it->second.end());
    }

  // Get useful heterozygous SNPs
  std::ifstream ifs_hets(in_hets.c_str());
  bool first = true;
  refmap snp_map;
  while(std::getline(ifs_hets, line))
    {
      if(first)
	{
	  first = false;
	  continue;
	}
      std::string col;
      std::stringstream ss(line);
      //unsigned short fn = 0;
      het_record r;
      for(ushort i = 0; i < 4; i++)
        {
	  std::getline(ss, col, '\t');
          if(i == 0) r.chrom = col;
          if(i == 1) r.pos = std::stoul(col);
          if(i == 2) 
	    {
	      ulong p = col.find_first_of('(');
	      r.ref = col[p+1];
	      r.alt = col[p+3];
	    }
	  
        }

      if(in(hom_map, r.chrom, r.pos))
	{
	  snp_map[r.chrom].push_back(std::pair<ulong, char>(r.pos, r.ref));
	}
    }

  for(auto it = snp_map.begin(); it != snp_map.end(); it++)
    {
      std::sort(it->second.begin(), it->second.end());
    }

  /*  
  for(auto it = snp_map.begin(); it != snp_map.end(); it++)
    {
      std::string chr = it->first;
      for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
	{
	  std::cout << chr << '\t' << it2->first << '\t' << it2->second << std::endl;
	}
    }
  
  for(auto it = regions.begin(); it != regions.end(); it++){
    std::cout << it->chrom << '\t' << it->start << '\t' << it->stop << std::endl;
    }
  */


  BGZF* bamf = bgzf_open(in_bam.c_str(), "r");
  htsFile* out_ref = hts_open("out_ref.sam", "w");
  htsFile* out_alt = hts_open("out_alt.sam", "w");
  bam_hdr_t * header = bam_hdr_init();   
  header = bam_hdr_read(bamf);    
  //Read 1 read and loop   
  int bam_read_ret; 
  //Return value (length of read if good, -1 otherwise)   
  bam1_t * curr_read = bam_init1();   
  while((bam_read_ret = bam_read1(bamf, curr_read)) != -1){

    std::string chrom(header->target_name[curr_read->core.tid]);
    unsigned long pos = curr_read->core.pos + 1;

    bool in_region = false;
    
    for(auto it = region_map[chrom].begin(); it != region_map[chrom].end(); it++)
      {
	if(pos >= it->first && pos <= it->second)
	  {
	    in_region = true;
	    break;
	  }
      }

    if(! in_region) continue;

    int32_t l_qseq = curr_read->core.l_qseq;
    uint8_t * seq = bam_get_seq(curr_read);

    for(auto it = snp_map[chrom].begin(); it != snp_map[chrom].end(); it++)
      {
	if(it->first >= pos && it->first <= pos + l_qseq)
	  {
	    ulong diff = it->first - pos;
	    if(it->second == bam2nuc(bam_seqi(seq, diff)))
	      {
		if(sam_write1(out_ref, header, curr_read) != 0)
		  {
		    //std::cerr << "Error writing record" << std::endl;
		  }
	      }
	    else 
	      {
		if(sam_write1(out_alt, header, curr_read) != 0) 
		  {
		    //std::cerr << "Error writing record" << std::endl;
		  }
	      }
	    break;
	  }
      }
  
  }

  hts_close(out_ref);
  hts_close(out_alt);
  
}
