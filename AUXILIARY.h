#ifndef AUX
#define AUX

#include<string>
#include<map>
#include<vector>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<iostream>
#include<math.h>
#include<typeinfo>
#include<cassert>
#include<sys/stat.h>
#include <string_view>
#include <charconv>
//#include<filesystem> //not available at BNL

//using namespace std; // bad style.

/*
Changed on 10.02.2022



*/

// VOID TYPE ###################################################################################################################################

void checkOFile(std::ofstream & ofstr, const std::string &adress); // check if ofstream is good.

void checkIFile(std::ifstream & ifstr, const std::string &adress); // check if ifstream is good.

void line_to_vec(std::string line, std::vector<std::string> &stringVector, const char &delimiter); // take line like a	b	c	d and put it into vector = [a,b,c,d] which is still string as you can see.

void line_to_vec(std::string line, std::vector<std::string> &stringVector, const std::string &delimiter); //same priciple but take string instead of char.

//void readQ2Grid(std::vector<double> & Q2binLowVec, std::vector<double> & Q2binHighVec, std::vector<double> & Q2binVec, const std::string &adress, const char &delimiter); //expects sorted file of Q2 bins

void getline(std::ifstream &inc, std::string& line, int linecounter); // get a line into the string line argument, but the function is a return void type. ( like the "original")

void getline(std::ifstream &, std::vector<std::string> &, const char); //getline of stream, but already split it(by char) to vector and store it there.

void skiplines(std::ifstream&, const unsigned int); // skip a few lines.

void fileopen(std::ifstream&, const std::string &file,const unsigned int); // skip a few lines.

void progress2terminal(const unsigned int ctr, const unsigned int total, const double ratio);// print to terminal the progress, if you have more than certain % of total events progressed

void cut_file(const std::string &src, const std::string &target, unsigned int firstline, unsigned int lastline); // truncate top and bottom of a file and create newfile

void transpone_file(const std::string &file); // transpones file, switch lines and rows

void vec_cpy_sd(std::vector<std::string> & s_vec, std::vector<double> & d_vec); //copy a string vector to a double vector

void getline_asVec_withEntry_atPosition(std::ifstream & f, std::vector<std::string> & vec, const char delimiterSymbol, std::string & entry, const unsigned int i);//this function should return the line (as a vector<string>) which contains the entry "entry" at the coloumn i. It will use just getline, thus, if the ifstream is new opened, it will take the first line it finds that matches, but if in more general application this should also return simply the next line that matches the entry-position condition that is found by the ifstream.

void true_false_input(bool & value, std::string& question); // std::cin to ask whether the something is true or not. 

void extract2_1(std::string & file, std::string &newfile, std::vector<int> &positions, std::vector<double> & values, const char delimiter); //check function body, no idea right now

void getBin(const int line, std::ifstream &in, double &varBin, double &varLow, double &varHigh); // adapted to specific csv conversion.

void round_double_vec(std::vector<double>& , const int n); // rounds double to n digits behind comma.

void error(std::string text); // removed the errorcode ->, int errorcode); // prints text and exits

void warning(std::string text); // prints text and continues

void shell_text(std::string text); // prints text and continues

void save_values_to_file(std::ofstream &f, std::vector<double> &v); // converts double - vector to string separated with TAB (in original version) and prints this line + std::endl to ofstream.

void append_values_to_file(const std::string & filename, std::vector<double> &v); // appends (by calling "save_values_to_file") the vector to the file

void file_append_1to2(const std::string& file1, const std::string& file2); //appends file 1 to file 2

void cprint(std::string, char colour); // print with colour that is mapped by char.

void print_dif(std::vector<double> &ov, std::vector<double> &nv); // print the difference of two vectors

// CHAR TYPE ###################################################################################################################################

char get_Parameter_char(std::string file, std::string parametercode); // get parameter as char, see: std::string get_Parameter

// BOOL TYPE ###################################################################################################################################

bool check_File_exists(const std::string & file); // true if file exists, false if not.

bool getline_as_vec(std::ifstream &ifs, std::vector<std::string> &line_vec, const char delimiter); //to use the above function in while-loop like while(getline(...)).

bool IsPathExist(const std::string &directory); // true if dir exists, false if not. Function body copied from https://www.codegrepper.com/code-examples/cpp/how+to+find+out+if+a+directory+exists+in+cpp

bool get_Parameter_bool(std::string file, std::string parametercode); // get parameter as char, see: std::string get_Parameter

// INT TYPE ###################################################################################################################################

int get_Parameter_int(std::string file, std::string parametercode); // get parameter as integer, see: std::string get_Parameter

int modsys(std::string cmd); // allows for "modsys("mkdir " + dirname);" for example. function calls "system(cmd.c_str())" and returns. thats all. ( returns status in integer)

// UNSIGNED INT TYPE ###################################################################################################################################

// LONG UNSIGNED INT TYPE ###################################################################################################################################

long unsigned int final_line_ctr(const std::string&); // counter of final line in document

long unsigned int final_nonempty_line_ctr(const std::string&); // final line that is not empty. calls final_line_ctr and moves upwards (but "moving upwards" is calling a function that reads the line above. but currently, to read line(i) it reads from beginning to line i. thus, for a long document, this is very inefficient. Most likely there are much better solutions. For small documents and preambel or finalizing functions this is irrelevant.

long unsigned int first_line_without_signature(const std::string& file, const std::string& signature); // first line that does not start on a certain string.

long unsigned int polyadisch_inverse(const unsigned int basis, const std::vector<unsigned int>& representation); // inverse operation : create number from rep.



// STD::STRING TYPE ###################################################################################################################################

std::string vec_to_line(std::vector<std::string>&, const char delimiter); // inverse line_to_vec (but use std::string as return type, which is more direct

std::string vec_to_line(const std::vector<double> &vec, const char delimiter); // like for a vector<string> but casts to double to string with std::to_string() automatically.

std::string vec_to_line(const std::vector<int> &vec, const char delimiter); // cast int to string, call vec_to_line(string, char)...

std::string vec_to_line_precision(std::vector<double> &vec, const char delimiter, int precision); // like vec_to_line but here it enforces a given precision for the conversion..

std::string vec_to_line_scientific(std::vector<double> &vec, const char delimiter, int precision); // like vec_to_line but here it enforces a given precision for the conversion and use scientific notation

//std::string scientific(double d, int dec); // convert a double to a string that has abc_e_+-_n format.

//std::string getline(std::ifstream &inc, std::string& line, int linecounter); //get a line in a file with the specific linenumber (linecounter)

std::string remove_chars(std::string source, std::string chars); // remove all charactrers appearing in #chars from #source. return the remaints of #source. Not fully tested so may contain a lot of bugs.         

std::string get_Parameter(std::string file, std::string parametercode); //get a paramter from a file

std::string get_Parameter_L(const std::string& line, std::string parametercode); // takes line and not a file

std::string TrueFalse(bool); // return "True" or "False";

// DOUBLE TYPE ###################################################################################################################################

double average1(double x0, double x1, std::vector<double> grid, std::vector<double> vals); // function to compute weighted average over a grid with values assigned to it. all one dimensional here

double computeBjorkenx(const double nu, const double Q2, const double M); //function i needed to compute xB from the variables given (nu, Q2 and Mass).

double get_Parameter_double(std::string file, std::string parametercode); // get parameter as double, see: std::string get_Parameter

// STD::VECTOR<INT> TYPE ###################################################################################################################################

std::vector<int> v_stoi(const std::vector<std::string> & svec); //stoi but for vectors

// STD::VECTOR<UNSIGNED INT> TYPE ###################################################################################################################################

std::vector<unsigned int> polyadisch(const unsigned int basis, const unsigned int number); // create vector containing polyadische Zerlegung of number to basis, whith the respective names

// STD::VECTOR<DOUBLE> TYPE ###################################################################################################################################

std::vector<double> v_stod(const std::vector<std::string> & svec); //stod but for vectors

std::vector<double> matrix_col(std::vector<std::vector<double>> matrix, const unsigned int col); // get the coloumn "col" from a std::vector<std::vector<double>> type matrix. to get the row, just access the position like matrix[i]

std::vector<double> logGrid(double minVal, double maxVal, double basis, int steps, char mode); // create log grid with range min-maxVal, have exponential basis "basis" and have steps, total (mode = 't') or relative (mode = 'r') to the basis step.

std::vector<double> linGrid(double minVal, double maxVal, int steps); // create linear grid with range min-maxVal and number of steps. 1 step will be one interval

std::vector<double> linGrid(double minVal, double maxVal, double steps_width); // create linear grid with range min-maxVal and defined step_width

std::vector<double> grid_around_vector(std::vector<double>&); // create grid that fits the values of the vector given in it. not very accurate.

std::vector<double> empty_vector(unsigned int size); // only works for double since template failed?? Should not there exist a constructor for this?

std::vector<double> get_Parameter_bin(std::string file, std::string parametercode, std::string separator); // get a pair of parameters separated by some string. this is bound to double values.

std::vector<double> difference(std::vector<double> &v1, std::vector<double>& v2); // diff = v1 - v2

std::vector<double> readbins(std::string file); // create vector, read bins from file, using existing methods, returns as vector. reduces lines of code in effective code.

// STD::VECTOR<STD::STRING> TYPE ###################################################################################################################################

std::vector<std::string> stov(const std::string&, const std::string& delimiter); // string to vec

std::vector<std::string> get_Parameter_bin_L(const std::string& line, std::string parametercode, std::string delimiter); //same as get_Parameter_bin but takes line instead of file as first argument

std::vector<std::string> line_to_vec(const std::string line, const char delimiter); // a line to vec method

std::vector<std::string> line_to_vec(const std::string line, const char delimiter, char mode); // a line to vec method. allows for multiple modes, including: 's' single char as delimiter -> refers to previous method | 'u' undefined amount of the same char type as delimiter

std::vector<std::string> file_to_vec(const std::string& file); //read a whole file and the lines will become entries of the vector. this probably only makes sense for small files..

// STD::VECTOR<STD::VECTOR<TYPE> TYPE ###################################################################################################################################

std::vector<std::vector<double>> file_to_matrix(const std::string &file, const char delimiter); // read a whole file ( THAT MUST ONLY CONSIST OF NUMBERS!! ) to a matrix: vec<vec<double>>. I updated it so now it is necessary to pass a delimiter ( before it was by just TAB automatically). I can add a function that only takes the file as an argument and automatically calls this function with TAB as second argument.

std::vector<std::vector<std::string>> file_to_stringmatrix(const std::string &file); // read a whole file to a matrix: vec<vec<string>>.
// TEMPLATES ###################################################################################################################################


// TYPENAME RETURN VALUE TYPE TEMPLATES

template <typename T>
T sum_vector(const std::vector<T>&v){
	T s = 0;
	for(T x : v){
		s += x;
	
	}
	return s;
}


template <typename T>
T find_minimum(std::vector<T> &vec)
{
	if(vec.size() == 0){
		//std::cout << "Error: the vector in which the minimum should be found is empty." << std::endl;
		//std::exit(1);
		error("Error: the vector in which the minimum should be found is empty.");
		std::exit(3);
	}
	T val = vec[0];
	for(size_t i = 1; i < vec.size(); ++i)
	{
		if(val > vec[i]) val = vec[i];
	}
	return val;
}

template <typename T>
T find_maximum(std::vector<T> &vec)
{
	if(vec.size() == 0){
		//std::cout << "Error: the vector in which the maximum should be found is empty." << std::endl;
		//std::exit(1);
		error("Error: the vector in which the maximum should be found is empty.");
		std::exit(3);
	}
	T val = vec[0];
	for(size_t i = 1; i < vec.size(); ++i)
	{
		if(val < vec[i]) val = vec[i];
	}
	return val;
} 



// VOID TYPE TEMPLATES ###################################################################################################################################


template <typename T>
void vec_to_file(std::vector<T> v, std::string file ){

	std::ofstream f;
	f.open(file);
	for(unsigned int i = 0;  i < v.size(); ++i){
		f << v[i] << std::endl;
	}
	f.close();
	return;
}

template <typename T>
void print(const std::vector<T>& v, std::ostream& os, const char delimiter){
	for(unsigned int i = 0; i < v.size()-1; ++i){
		os << v[i] << delimiter;
	}
	os << v[v.size()-1];
	return;
}


template <typename T>
void print_double_type_vector_scientific(const std::vector<T>& v, std::ostream& os, const char delimiter, const int nGEZ){
	
	
	//std::vector<std::string> sv = v_to_string(v);
	std::stringstream s;
	s << std::scientific << std::setprecision(nGEZ);
	for(unsigned int i = 0; i < v.size()-1; ++i){
		s << v[i] << delimiter;
		//os << v[i] << delimiter;
	}
	s << v[v.size()-1];
	std::string line = s.str();
	os << line;
	return;
}

template <typename T>
void flip(std::vector<T> &vec)
{
//	std::cout << "in flip ( auxlib: inverse order of a vector)" << std::endl;
	if(vec.size() == 1) // inverse of 1 element vector is trivial
		return;
		
	std::vector<T> cpy = vec; // a copy
//	std::cout << " vector copied " << std::endl;
	for(unsigned int i = 0; i < vec.size(); ++i)
	{
//		std::cout << i << std::endl;
		vec[i] = cpy[vec.size() - i - 1]; // inverts order of vector;
//		std::cout << vec[i] << std::endl;
	} 
	return;
}

template <typename T>
void print1dvec(std::vector<T> &vec)
{
	if(vec.size() == 0){
	
		//std::cout << "the vector you asked to print is empty" << endl;
		warning("the vector you asked to print is empty");
		return;
		
	}
	
	//old:
	/*for(unsigned int i = 0; i < vec.size(); ++i)
	{
		std::cout << vec[i] << std::endl;
	}
	*/
	//new: in a line
	std::cout << vec[0];
	for(unsigned int i = 1; i < vec.size(); ++i)
	{
		std::cout << "	" << vec[i]; // tab inserted
	}
	std::cout << std::endl;
	return;
}

template <typename T>
void print1dvec(const std::vector<T> &vec)
{
	std::vector<T> auxvec = vec;
	print1dvec(auxvec);
	return;
}

template <typename T>
void write_vec_to_stream(std::ofstream &f, std::vector<T> &vec, const char delimiter)
{
	if(vec.size() == 0){
	
		//std::cout << "asked to write an empty vector to ofstream" << endl;
		warning("asked to write an empty vector to ofstream");
		return;
		
	}
	
	f << vec[0];
	for(unsigned int i = 1; i < vec.size(); ++i){
			f << delimiter;
			f << vec[i];
		}
		f << std::endl;
	return;
}

template <typename T>
void read_1dFile_to_vec(const std::string &file, std::vector<T> &vec) // read a 1d file (so its just one number per line) into a vector, and since its a template i can use it for double, int, (so far thats all i have)
{
	std::ifstream f;
	f.open(file);
	checkIFile(f, file); // check if file is good.
	
	//safety checks:
	//checkIFile(f,file);
	if(vec.size() != 0)
	{
//		std::cout << "Warning: the vector given was not empty" << std::endl;
		warning("Warning: the vector given was not empty");
	}
	std::string line = "emptyline";
	
	if(typeid(T) == typeid(int))
	{
		while(getline(f,line))
		{
			vec.push_back(stoi(line));
		}
	}
	
	if(typeid(T) == typeid(double))
	{
		while(getline(f,line))
		{
			vec.push_back(stod(line));
		}
	}
	
	f.close();
	return;
}

// BOOL TYPE TEMPLATES

template <typename T>
bool onList(const T &item, std::vector<T> &list) // check if some value (item) is on a list (list)
{
	for(unsigned int i = 0; i < list.size(); ++i)
	{
		if(list[i] == item)
			return true;
	}
	return false;
}

template <typename T>
bool check_sorted(std::vector<T> &vec) // check if the vector values are increasing
{
	for(unsigned int i = 0; i < vec.size()-1; ++i)
	{
		if(vec[i] >= vec[i+1])
		{
			return false;
		}
	}
	return true; 
}

template <typename T>
bool isOrdered(std::vector<T> &vec, const char c){
	
	if(vec.size() == 0){
		//std::cout << "Warning: the vector to be checked for ordering is empty" << std::endl;
		warning("Warning: the vector to be checked for ordering is empty");
	}
	
	size_t start = 0, step = 0; // irrelevant default values.
	
	switch(c){
		case '<': //increasing values
			start = 1;
			step = +1;
			break;		
		case '>': //declining values
			start = vec.size()-1;
			step = -1;
			break;
		default:
			//std::cout << "Error: the symbol\"" << c << "\" is not suited to check the ordering of a vector. Use \"<\" for increasing or \">\" for decreasing values." << std::endl;
			//std::exit(1);
			error("Error: the symbol\"" + std::to_string(c) + "\" is not suited to check the ordering of a vector. Use \"<\" for increasing or \">\" for decreasing values.");
			std::exit(4);
	}
	
	T val0 = vec[start];
	for(size_t i = start + step; i < vec.size() && i >= 0; i+=step){ //two conditions to include increasing and declining values.
		
		T val1 = vec[i];
		if(val0 > val1) return false;
			
	}
	return true;
}

// UNSIGNED INT TYPE TEMPLATES ###################################################################################################################################

template <typename T>
T getMaxValue(const std::vector<T> &myvec){
	T maxval = myvec[0];
	
	for(unsigned int i = 1; i < myvec.size(); ++i){
	
		if(maxval < myvec[i])
			maxval = myvec[i];
	}
	return maxval;
}

template <typename T>
T getMinValue(const std::vector<T> &myvec){
	T minval = myvec[0];
	
	for(unsigned int i = 1; i < myvec.size(); ++i){
	
		if(minval < myvec[i])
			minval = myvec[i];
	}
	return minval;
}

unsigned int getIndex(std::vector<std::string> &vec, std::string entry); // this is here because its just a hard copy of the below function, the template, for strings, because c++ complains if it should convert a string into a string via "std::to_string(s)".

template <typename I>
unsigned int getIndex(std::vector<I> &vec, I entry)
{
	int position = 0; // definition here is pointless but not defined rivet compilation (due to flags) returns an error.
	bool already_found = false;
	for(unsigned int i = 0; i < vec.size(); ++i)
	{
		if(vec[i] == entry)
		{
			if(!already_found)
			{
				position = i;
				already_found = true;
				break; // if this is used, the warning is blocked. if i comment it out i access the warning, this i updated on 13.04.2023
			}
			else
			{
			
				//cout << "the vector seems to contain at least two times the entry \"" << entry << "\" at the positions " << position << " and " << i << ".\nAbort." << endl;
				//std::exit(2); better leave it as a warning but keep running:
				//std::cout << "the vector seems to contain at least two times the entry \"" << entry << "\" at the positions " << position << " and " << i << ".\nThe entry: " << entry << std::endl << "The vetcor" <<  std::endl;
				warning("the vector seems to contain at least two times the entry \"" + std::to_string(entry) + "\" at the positions " + std::to_string(position) + " and " + std::to_string(i) + ".\nThe entry: " + std::to_string(entry) + "\nThe vetcor");
				print1dvec(vec);
				
			}
		}	
	}
	if(already_found)
	return position;
	else
	{
		error("the vector does not seem to contain the entry \"" + std::to_string(entry) + "\". ");
//		std::cout << "the vector does not seem to contain the entry \"" << entry << "\". " << endl;
		std::cout << "the vector is: " << std::endl;
		std::cout << vec << std::endl;
		//error("abort");
		std::exit(3);
	}
}


template <typename T>
T find_maximum_index(std::vector<T> &vec)
{
	T entry = find_maximum(vec);
	unsigned int i = getIndex(vec, entry);
	return i;
}

template <typename I>
unsigned int getBinIndex(std::vector<I> &vec, I entry)
{
	//assume i am no idiot and the vector binning contains this entry
	
	for(unsigned int i = 0; i < vec.size()-1; ++i)
	{
		if(vec[i] <= entry && vec[i+1] > entry )
		{
			return i;
		}
	}
	//if the entry is exactly the very very last border, in this case return the last bin!
	if(entry == vec[vec.size()-1])
		return vec.size()-1;
	
	error("the vector binning does not seem to contain the entry \"" + std::to_string(entry) + "\". ");
	//std::cout << "the vector binning does not seem to contain the entry \"" << entry << "\". " << std::endl;
	std::cout << "the vector is: " << std::endl;
	std::cout << vec << std::endl;
	//error("abort");
	std::exit(3);
}

// STD::OSTREAM TYPE TEMPLATES ###################################################################################################################################

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
	os << '[';
	for(unsigned int i = 0; i < v.size(); ++i){
		os << ' '<< v[i];
	}
	os << ' ' << ']';
	return os;
}


// STD::VECTOR<STD::STRING> TYPE TEMPLATES ###################################################################################################################################

template <typename T>
std::vector<std::string> v_to_string(const std::vector<T>& v){
	
	std::vector<std::string> v_s;
	
	for(unsigned i = 0; i < v.size(); ++i){
		v_s.push_back(std::to_string(v[i]));
	}
	return v_s;
}


// STD::VECTOR<TEMPLATE> TYPE TEMPLATES ###################################################################################################################################
template <typename T>
std::vector<T> operator+(const std::vector<T> &v1, const std::vector<T>& v2){
	assert(v1.size() == v2.size());
	
	std::vector<T> result(v1.size(), 0); // will work for integer and double.
	//result.reserve(v1.size());
	for(size_t i = 0; i < result.size(); ++i){
		result[i] = v1[i] + v2[i];
	}
	return result;
} 

template <typename T>
std::vector<T> get_vec_entries_from_file(std::string adress)
{
	std::vector<T> vec;
	
	std::ifstream f;
	f.open(adress);
	
	checkIFile(f, adress);
	
	if(vec.size() != 0)
	{
		//std::cout << "Warning: the vector given was not empty" << std::endl << "Code 57192" << std::endl;
		warning("Warning: the vector given was not empty");
	}
	std::string line = "emptyline";
	
	if(typeid(T) == typeid(int))
	{
		while(getline(f,line))
		{
			vec.push_back(stoi(line));
		}
	}
	
	if(typeid(T) == typeid(double))
	{
		while(getline(f,line))
		{
			vec.push_back(stod(line));
		}
	}
	
	f.close();
	
	
	return vec;
}

template <typename T>
std::vector<T> order_vector(const std::vector<T> &vec, const char c){
	
	if(vec.size() == 0){
		//std::cout << "Warning: the vector to be ordered is empty" << std::endl;
		warning("Warning: the vector to be ordered is empty");
	}
	
	bool flip_vector = false;
	
	switch(c){
		case '<': //increasing values
			flip_vector = false;
			break;
		case '>': //declining values
			flip_vector = true;
			break;
		default:
			//std::cout << "Error: the symbol\"" << c << "\" is not suited to check the ordering of a vector. Use \"<\" for increasing or \">\" for decreasing values." << std::endl;
			
			error("Error: the symbol\"" + std::to_string(c) + "\" is not suited to check the ordering of a vector. Use \"<\" for increasing or \">\" for decreasing values.");
			std::exit(4);
	}
	
	// expect that this vector has been checked whether is ordered or not, and there is a reason why this function is being called.
	
	
	std::vector<T> vec_out = vec; // create copy.
	std::vector<T> toy = vec; //	
	// just order it increasingly. if you want a decreasing vector, flip the vector.
	
	for(size_t i = 0; i < vec.size(); ++i){
		T val = find_minimum(toy);
		
		toy.erase( toy.begin() + getIndex(toy, val) );
		
		vec_out[i] = val;
		
	}
	
	//after this a vector should be orderd increasingly
	
	if(flip_vector) 
		flip(vec_out);
	
	return vec_out;
		
}

std::vector<int> load_distribution_fromfile_int(std::string f); //because template fails..
std::vector<double> load_distribution_fromfile_double(std::string f); //because template fails..
/*
template <typename T>
std::vector<T> load_distribution_fromfile(std::string f){
	std::vector<T> r;
	std::vector<std::string> file = file_to_vec(f);
	for( std::string s : file){
		const std::string_view str(s);
//	for(auto s : file){
		T v;
		auto [ptr, ec] {std::from_chars(str.data(), str.data() + str.size(), v)};
		//auto v = std::from_chars(
		if(ec != std::errc()){
			error("Problem with converting the string " + s + " into a numeric in load_distribution_fromfile(std::string )");
		}
		else
			r.push_back(v);
	}
	return r;
}*/

// OTHER SUCH AS MATRIX OPERATIONS

template <typename T>
std::vector<std::vector<T>> transpos_matrix(const std::vector<std::vector<T>>& matrix){



	std::vector<std::vector<T>> transposed_matrix; // copy (does it work)?
//	for(unsigned int j = 0; matrix[0].size(); ++j) // lines
	for(unsigned int i = 0; i < matrix[0].size(); ++i){ // columns
//		std::cout << i << std::endl;
		std::vector<T> v;
		for(unsigned int j = 0; j < matrix.size(); ++j){ // rows
			/*std::cout << " " << j << std::endl;
			std::cout << matrix[j][i] << std::endl;*/
			v.push_back(matrix[j][i]);

			
		}
		transposed_matrix.push_back(v);
	}
	return transposed_matrix;
}


//THIS GUY actually must appear below the template function "write_vec_to_stream"...
template <typename T>
void matrix_to_file(const std::vector<std::vector<T>>& matrix, const std::string & file, const char delimiter){ // demand delimiter as argument as well as well.

	std::ofstream f(file);
	checkOFile(f, file);
	
//	const char delimiter = '	';
	
	const unsigned int nrows = matrix.size();
	
	for(unsigned int i = 0; i < nrows; ++i){ // rows
		std::vector<T> row = matrix[i];
		write_vec_to_stream(f, row, delimiter);
	}	
	f.close();

	return;
}

/* THE FOLLOWING FUNCTION DID NOT WORK. 
template <typename T>
T read_1dFile_to_vec_direct(const std::string &file) // read a 1d file (so its just one number per line) into a vector, and since its a template i can use it for double, int, (so far thats all i have)
{
	std::ifstream f;
	f.open(file);
	checkIFile(f, file); // check if file is good.
	
	string line = "";
	
	if(typeid(T) == typeid(std::vector<int>))
	{
		std::vector<int> vec;
		while(getline(f,line))
		{
			vec.push_back(stoi(line));
		}
		f.close();
		return vec;
	}
	
	if(typeid(T) == typeid(std::vector<double>))
	{
		std::vector<double> vec;
		while(getline(f,line))
		{
			vec.push_back(stod(line));
		}
		f.close();
		return vec;
	}
	
	
}
*/



#endif
