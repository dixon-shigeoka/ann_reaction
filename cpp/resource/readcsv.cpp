#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <iostream>

//using namespace std;

int main()
{
  // read the statistic data from csv files //

  std::ifstream ifs("statistics.csv");
  std::string line;
  double statistics[2][10];
  const std::string delim = ",";

  int row = 0;
  int col = 0;
  std::string field;

  while ( std::getline(ifs, line) ) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::vector<std::string> strvec;
    while (getline(stream,field,',')){
      printf("%s\n",field.c_str());
      result.push_back(field);
      strvec = result;
    }
    for (int i=0; i<strvec.size();i++){
      statistics[col][i] = std::stod(strvec.at(i));
      printf("%2.19lf\n",std::stod(strvec.at(i)));
    }
    col++;
  }
}
