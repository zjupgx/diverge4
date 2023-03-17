#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <vector>

#define appname "genePlei"


#define version "V1.0"

class summary_t {
public:
  std::string name;
  std::vector<double> values;
};

class result_t {
public:
  int pos;
  std::vector<double> values;
};

#endif
