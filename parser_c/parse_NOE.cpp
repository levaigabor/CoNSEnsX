#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;
using std::getline;

#include <vector>
using std::vector;

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;



int main()
{
  // create a file-reading object
  ifstream fin;
  fin.open("data (Kopie).txt"); // open a file
  if (!fin.good())
    return 1; // exit if file not found

  // read each line of the file
  while (!fin.eof())
  {
    // read an entire line into memory
    string line;
    getline(fin,line);

    if (line.empty() || line[0] == '#')
    {
        continue;
    }

    boost::iterator_range<string::iterator> r(line.begin(), line.end());
    vector<iterator_range<string::const_iterator> > result;

    algorithm::split(result, r, is_any_of(" \n\t"), algorithm::token_compress_on);

    cout << "SIZE: " << result.size() << endl;

    for (auto i : result)
        cout << "value: '" << i << "'\n";

  }
}

// #include <boost/python.hpp>

// BOOST_PYTHON_MODULE(parse_NOE)
// {
//     using namespace boost::python;
//     def("printPrimes", printPrimes);
// }

