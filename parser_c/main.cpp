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



const string NOE_key   = "_Gen_dist_constraint_list.Constraint_type";
const string NOE_value = "NOE";

const string ID        = "_Gen_dist_constraint.ID";
const string Seq_ID_1  = "_Gen_dist_constraint.Seq_ID_1";
const string Seq_ID_2  = "_Gen_dist_constraint.Seq_ID_2";
const string Comp_ID_1 = "_Gen_dist_constraint.Comp_ID_1";
const string Comp_ID_2 = "_Gen_dist_constraint.Comp_ID_2";
const string Atom_ID_1 = "_Gen_dist_constraint.Atom_ID_1";
const string Atom_ID_2 = "_Gen_dist_constraint.Atom_ID_2";
const string Dist_val  = "_Gen_dist_constraint.Distance_upper_bound_val";

const string loopstart = "loop_";
const string stop      = "stop_";


int main()
{
  // create a file-reading object
  ifstream fin;
  fin.open("data.txt"); // open a file
  if (!fin.good())
    return 1; // exit if file not found

  bool NOE_found = false;
  bool in_loop   = false;
  bool loop_start = false;




  // read each line of the file
  while (!fin.eof())
  {
    // read an entire line into memory
    string line;
    getline(fin,line);

    // skip comment and empty lines
    if (line.empty() || line[0] == '#')
        continue;


    // store tokens in vector
    boost::iterator_range<string::iterator> r(line.begin(), line.end());
    vector<iterator_range<string::const_iterator> > result;

    algorithm::split(result, r, is_any_of(" \n\t"),
                     algorithm::token_compress_on);

    // find NOE data information
    if (std::find(result.begin(), result.end(), NOE_value) != result.end()) {
        cout << "NOE DATA FOUND\n";
        NOE_found = true;
    }


    loop_start = std::find(result.begin(),
                           result.end(),
                           loopstart) != result.end();

    if (loop_start)
        in_loop = true;

    if (NOE_found && in_loop)
    {
        cout << "Found loop in NOE!" << endl;
    }

    cout << line << endl;


    // cout << "SIZE: " << result.size() << endl;

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

