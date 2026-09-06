#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <vector>

#define appname "genePlei"



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
