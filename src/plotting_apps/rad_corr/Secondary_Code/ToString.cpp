#include <sstream>

using namespace std;

TString ToString(int num){
  std::ostringstream start;
  start << num;
  string start1 = start.str();
  return start1;
};
