#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>

#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>

#ifndef WIN32
#include <ctype.h>
#include <unistd.h>
#endif
#include <cstring>
#include "sequence.h"



//----------------------------------------------------------------------

using namespace std;


extern vector<int> kept; //the non-gap position

//----------------------------------------------------------------------


//This file contains the defintions of methods which are related to 
//the operations on gene sequence, such as loading sequence, gap cleaning, etc.


//load_sequences() is used to load the specific gene file into the 
//variable - sequences. The working procedure of this method is as follows:
// 1. clean the sequences first
// 2. load the first 1024 byte of file to determine the type of sequence file
//    (FASTA or CLUSTAL). The detailed information of file format can be found
//    in the following link:
//           http://www.umdnj.edu/rcompweb/PA/Notes/GenbankFF.htm
// 3. According to the file type, call the correspoding method (load_fasta_sequences
//    or load_clustal_sequences).
// 4. Double check the size of each sequence, and make sure they are the same.
/*
The purpose of this function is to generate a new, re-ordered sequences based on 
the input: tree_leaf_names. The function is mainly used for tree reroot.
*/
bool reorder_sequences(const std::vector<std::string> &tree_leaf_names, const std::vector<sequence_t> &orig_sequences,
					   std::vector<sequence_t> &new_sequences) {
	int j;

	new_sequences.clear();

	if(tree_leaf_names.size() != orig_sequences.size())
		return false;

	for( j = 0; j < tree_leaf_names.size(); j++ ) {
		vector<sequence_t>::const_iterator i;
		for(i=orig_sequences.begin(); i!=orig_sequences.end(); i++) {
			if ( tree_leaf_names[j] == i->label ) {
				new_sequences.push_back(sequence_t(*i));
				break;
			}
		}

		if(new_sequences.size() <= j) //did not find match sequence name
			return false;
	}

	return true;
}



 
bool load_sequences(string filename, vector<sequence_t> &sequences) {
  sequences.clear();

  //open the file
  FILE *file = fopen(filename.c_str(), "r");
  if(!file) {
    return false;
  }

  //load the first 1024 bytes of the file into variable - line
  char line[1024];
  char *tmp;
  tmp = fgets(line, 1024, file);
  fclose(file);

  //based on the line information to determine file type
  if(line[0] == '>') { //if first character is '>', then it is FASTA file format
    if(!load_fasta_sequences(filename, sequences)) {
      return false;
    }
  } else if(strstr(line, "CLUSTAL")) { //otherwise, it is a CLUSTAL format.
    if(!load_clustal_sequences(filename, sequences)) {
      return false;
    }
  } else {
    return false;
  }

  //check the length of each sequence, if they are not the same,
  //return false
  for(int i=1; i<(int)sequences.size(); i++) {
    if(sequences[0].sequence.size() != sequences[i].sequence.size()) {
      return false;
    }
  }
  
  return true;
}

//----------------------------------------------------------------------
//The method is used to load the sequences of CLUSTAL format. The CLUSTAL
//format has the following characteristics:
//1. The word "CLUSTAL" should be the first word on the first line of the file.
//2. The alignment is displayed in blocks of fixed length.
//3. Each line in the block corresponds to one sequence.
//4. Each sequence line starts with a sequence name followed by at least one space and then the sequence.
//
//BUG:
//Kent did not implement this method very well. Potentially, it will cause some 
//error while loading CLUSTAL format file. We need to improve this method a little
//bit.

bool
load_clustal_sequences(string filename, vector<sequence_t> &sequences) {
  sequences.clear();

  //open the file
  FILE *file = fopen(filename.c_str(), "r");
  if(!file) {
    return false;
  }


  const int line_len = 1024;
  char line[line_len];
  char * tmp;
  tmp = fgets(line, line_len, file);
  while(fgets(line, line_len, file)) {
    if(isspace(line[0])) continue;
    char seq_name[line_len], seq[line_len];
    if(sscanf(line, "%s %s", seq_name, seq) == 2) {
      char *p;
      do {
		  p = strchr(seq_name, '_');
		  if(p) *p = ' ';
      } while(p);
      
      vector<sequence_t>::iterator i;
      for(i = sequences.begin(); i != sequences.end(); i++) {
		  if(i->label == seq_name) {
			  i->sequence += seq;
			  break;
		  }
      }
      if(i == sequences.end()) {
		  sequences.push_back(sequence_t(seq_name, seq));
      }
    } //if
  } //while
  
  fclose(file);

  return true;
}

//----------------------------------------------------------------------
//This method is used to load the FASTA format sequence file into variable - 
//sequences. This method is used by the load_sequences() above. A demo file (CASP.fasta)
//can be found in the /demo directory. It uses a Finite Stat Automata (FSA) to 
//identify the format of sequence data
bool
load_fasta_sequences(string filename, vector<sequence_t> &sequences) {
	//clean the sequences first
  sequences.clear();
  
  //open the sequence file
  FILE *file = fopen(filename.c_str(), "r");
  if(!file) {
    return false;
  }

  enum { unknown, header, data } state = unknown;

  char c;
  char prev_c = '\n';

  sequence_t seq;
  
  int n;
  do { //Repeatly do the following loop for each character in the sequence file
    n = fgetc(file); 
    if(n != EOF) {
      c = (char)n; //read one character
      switch(c) {
	  case '>': //if current character is '>', then we encounter a new sequence
        if(prev_c == '\n') {
			if(seq.sequence.size() > 0) { //if previous sequence is not empty, push it
				sequences.push_back(seq); //into sequences and begin to process current
				seq.label = "";           //sequence
				seq.sequence = "";
			}
			state = header;
        }
        break;
      case '\n':
      case '\r': //encounter a newline, the state becomes data
        if(state == header || state == data) {
          state = data;
        }
        break;
	/*
      case ' ':
        if(state == data) {
          state = data;
        }
        break;
	*/
      default:
		  if(state == header) { //if state is header, concatenate the c to sequence 
			  if(c == '_') c = ' '; //label
			  seq.label += c;
		  }
		  if(state == data) { //otherwise, put it into the sequence data
			  seq.sequence += toupper(c);
		  }
		  break;
	  }
	  prev_c = c;
    }
  } while(n != EOF);

  //process the last sequence
  if(seq.sequence.size() > 0) {
    sequences.push_back(seq);
  }

  //close the sequence file
  fclose(file);

  // remove trailing space from the sequence name
  for(int i=0; i<(int)sequences.size(); i++) {
    int len = sequences[i].label.size();
    int j = len - 1;
    while(isspace(sequences[i].label[j]) && j>=0) { j--; }
    j++;
    if(j < len) {
      sequences[i].label.erase(j);
    }
  }

  return true;
}

//----------------------------------------------------------------------
//clean_gaps is used to clean the gaps in the sequences (of course). 
//Mechanism:
// 1. for each gene in all sequence, check its value. If it is 
//    one of four types ('-', '.', '?', 'X'), then store the column
//    position of of this gene.
// 2. for all gap columns, clean the whole column in all sequences

void
clean_gaps(vector<sequence_t> &sequences, vector<int> &removed) {
  removed.clear();
  
  if(sequences.empty()) return;
  
  int len = sequences[0].sequence.size();
  int num = sequences.size();
  for(int i=0; i<len; i++) {
    for(int j=0; j<num; j++) {
      if(sequences[j].sequence[i] == '.' ||
	 sequences[j].sequence[i] == '-' ||
	 sequences[j].sequence[i] == '?' ||
	 sequences[j].sequence[i] == 'X') {
	removed.push_back(i); //memorize the column position of gene gap
	break;
      }
    }
  }

  //erase the whole column out of the gene sequences
  vector<int>::reverse_iterator k;
  for(k = removed.rbegin(); k != removed.rend(); k++) {
    for(int j=0; j<num; j++) {
      sequences[j].sequence.erase(*k, 1);
    }
  }
}



//by zyw
//clean gaps and save the position of the kept aa
//modified from the above one
void
clean_gaps(vector<sequence_t> &sequences, vector<int> &removed, vector<int> &kept) {
	removed.clear();
	kept.clear();

	if(sequences.empty()) return;
	
	int len = sequences[0].sequence.size();
	int num = sequences.size();
	for(int i=0; i<len; i++) {
		for(int j=0; j<num; j++) {
			if(sequences[j].sequence[i] == '.' ||
				sequences[j].sequence[i] == '-' ||
				sequences[j].sequence[i] == '?' ||
				sequences[j].sequence[i] == 'X') {
				removed.push_back(i); //memorize the column position of gene gap
				break;
			}
			if(j==num-1) {
				kept.push_back(i);
			}
		}
	}

	//erase the whole column out of the gene sequences
	vector<int>::reverse_iterator k;
	for(k = removed.rbegin(); k != removed.rend(); k++) {
		for(int j=0; j<num; j++) {
			sequences[j].sequence.erase(*k, 1);
		}
	}

}





//----------------------------------------------------------------------
//This method is to find out the maximal removed gap position + 1
int
pos_w_gaps(const vector<int> &removed, int pos) {
  vector<int>::const_iterator i;
  for(i = removed.begin(); i != removed.end(); i++) {
    if(*i <= pos) pos++;
  }
  
  return pos;
}

//----------------------------------------------------------------------
//This method is pretty similar to the above one. The only difference is,
//we don't need to return the gap position information back to the caller method,
//which means we don't need removed parameter.
void
clean_gaps(vector<sequence_t> &sequences) {
  if(sequences.empty()) return;
  
  int len = sequences[0].sequence.size();
  int num = sequences.size();
  for(int i=len-1; i>=0; i--) {
    for(int j=0; j<num; j++) {
      if(sequences[j].sequence[i] == '.' ||
	 sequences[j].sequence[i] == '-' ||
	 sequences[j].sequence[i] == '?' ||
	 sequences[j].sequence[i] == 'X') {
	for(j=0; j<num; j++) {
	  sequences[j].sequence.erase(i, 1);
	}
	break;
      }
    }
  }
}

//----------------------------------------------------------------------
//This method can be used to cut the sequences using the specified parameter
//(start, len). The results will be stored in the variable - cut_seqs
bool
cut_sequences(const vector<sequence_t> &seqs, int start, int len,
	      vector<sequence_t> &cut_seqs) {
  cut_seqs.clear();//clean the variable first

  vector<sequence_t>::const_iterator i;
  //for each sequence, we copy a subsequence (cut) using the specified
  //parameter (start, len) and store it into cut_seqs;
  for(i=seqs.begin(); i!=seqs.end(); i++) {
    cut_seqs.push_back(sequence_t(i->label,
				  i->sequence.substr(start, len)));
  }
  
  return true;
}



//------------------------------------------------------------------------
//print sequence (for debug purpose only).
#if 0
void
print_seqs(const vector<sequence_t> &seqs) {
  vector<sequence_t>::const_iterator i;
  for(i=seqs.begin(); i!=seqs.end(); i++) {
    cout << '>' << i->label << '\n' << i->sequence << '\n';
  }
}

#endif
