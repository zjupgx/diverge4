#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <string>
#include <vector>

class sequence_t {
public:
  sequence_t() {}
  sequence_t(std::string name, std::string seq) : label(name), sequence(seq) {}
  sequence_t(const sequence_t &s) { label = s.label; sequence = s.sequence; }

  std::string label;
  std::string sequence;

  sequence_t &operator = (const sequence_t &s) {
    label = s.label;
    sequence = s.sequence;
    return *this;
  }
};

bool reorder_sequences(const std::vector<std::string> &tree_leaf_names, const std::vector<sequence_t> &orig_sequences,
					   std::vector<sequence_t> &new_sequences);

bool load_sequences(std::string filename, std::vector<sequence_t> &sequences);
bool load_clustal_sequences(std::string filename,
			    std::vector<sequence_t> &sequences);
bool load_fasta_sequences(std::string filename,
			  std::vector<sequence_t> &sequences);

void clean_gaps(std::vector<sequence_t> &sequences, std::vector<int> &removed);
int pos_w_gaps(const std::vector<int> &removed, int pos);
void clean_gaps(std::vector<sequence_t> &sequences);

bool cut_sequences(const std::vector<sequence_t> &seqs, int start, int len,
		   std::vector<sequence_t> &cut_seqs);
//void print_seqs(const std::vector<sequence_t> &seqs);

//add by zyw
void clean_gaps(std::vector<sequence_t> &sequences, std::vector<int> &removed, std::vector<int> &kept);
	
#endif
