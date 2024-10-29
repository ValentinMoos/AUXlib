#include"Aux.h"

void shell_text(std::string text){
	std::cout << "\033[37m" << text << "\033[0m" << std::endl;
	return;
}

void warning(std::string text){
	std::cout << "\033[33m" << text << "\033[0m" << std::endl;
	return;
}

void error(std::string text){
	std::cout << "\033[1m\033[31m" << text << "\033[0m" << std::endl;
	//exit(errorcode); // better do it by really call exit in that function, because c++ understands that it can end a function with "exit" and therefore does not need a return value in this case.
	
	
	/*
	errorcodes:
	0 would be no error, all good
	1 unknown file, string, something that was expected to be found was not
	2 conversion error?
	3 something not found in a vector
	4 user input error. user is me
	*/
}

void line_to_vec(std::string line, std::vector<std::string> &stringVector, const char &delimiter)
{	
	// clear vector:
	stringVector.clear();
	
	
	

	static const size_t noPos = -1; //copied from forum
	#ifdef DEBUG
	std::cout << "In function :line_to_vec" << std::endl;
	#endif
	
	if(line.find(delimiter) == noPos)
	{
		// did not find the symbol, hence i expect that the array is only one symbol
		//std::cout << "Warning: The string \"" << line << "\" should be transformed to a vector with the elements separated by \"" << delimiter << "\"\nSince no \"" << delimiter << "\" was found expect that the line is only one entry." << std::endl;
		std::string warning_msg = "Warning: The string \"" + line + "\" should be transformed to a vector with the elements separated by \"" + delimiter + "\"\nSince no \"" + delimiter + "\" was found expect that the line is only one entry.";
		warning(warning_msg);
		stringVector.push_back(line);
		return;
	}
	//if there are more that one event (thus the delimiter appears)
	for(int i=0; line.find(delimiter) != noPos; ++i) // condition: there exists a delimiter symbol 
	{
		#ifdef DEBUG
		std::cout << "In function :line_to_vec: 	for loop: i=" << i << std::endl;
		#endif
		size_t pos;
		pos = line.find(delimiter);
		stringVector.push_back(line.substr(0,pos));
		line = line.substr(pos+1,line.length()-(pos+1));
		//last symbol:
		if(line.find(delimiter) == noPos)
		{
			stringVector.push_back(line.substr(0,line.length()));
		}
		#ifdef DEBUG
		std::cout << "stringVector[" << i << "]	" << stringVector[i] << "	current line string:	" << line << std::endl;
		std::cout << "pos: " << pos << std::endl;
		std::cout << "next pos: " << line.find(delimiter) << std::endl;
		#endif
		
	}
	return;
}

std::vector<std::string> line_to_vec(const std::string line, const char delimiter){
	std::vector<std::string> sv;
	std::string copy = line;
	
	if(line.find(delimiter) == std::string::npos){
		warning("Warning: The string \"" + line + "\" should be transformed to a vector with the elements separated by \"" + delimiter + "\"\nSince no \"" + delimiter + "\" was found expect that the line is only one entry.");
		sv.push_back(line);
		return sv;
	}
	
	for(unsigned i = 0; copy.find(delimiter) != std::string::npos; ++i){
		
		unsigned nextPos = copy.find(delimiter);

		std::string rest = copy.substr(nextPos+1);

		std::string entry = copy.substr(0,nextPos);

		sv.push_back(entry);
		copy = rest;
		
	}
	//now only one entry left:
	if(copy.find(delimiter) != std::string::npos){
		error("After loop finished condition still satisfied in \"std::vector<std::string> line_to_vec(const std::string line, const char delimiter)\"");
		exit(3);
	}

	sv.push_back(copy);
	return sv;
}

void line_to_vec(std::string line, std::vector<std::string> &stringVector, const std::string &delimiter){	
	// clear vector:
	//stringVector.clear();
	
	if(line.find(delimiter) == std::string::npos)
	{
		// did not find the symbol, hence i expect that the array is only one symbol
		//std::cout << "Warning: The string \"" << line << "\" should be transformed to a vector with the elements separated by \"" << delimiter << "\"\nSince no \"" << delimiter << "\" was found expect that the line is only one entry." << std::endl;
		std::string warning_msg = "Warning: The string \"" + line + "\" should be transformed to a vector with the elements separated by \"" + delimiter + "\"\nSince no \"" + delimiter + "\" was found expect that the line is only one entry.";
		warning(warning_msg);
		stringVector.push_back(line);
		return;
	}
	//if there are more that one event (thus the delimiter appears)
	for(int i=0; line.find(delimiter) != std::string::npos; ++i) // condition: there exists a delimiter symbol 
	{
		size_t pos = line.find(delimiter);
		stringVector.push_back(line.substr(0,pos));
		line = line.substr(pos+delimiter.size(),line.length()-(pos+delimiter.size()));
		//last symbol:
		if(line.find(delimiter) == std::string::npos)
		{
			stringVector.push_back(line.substr(0,line.length()));
		}
	}
	return;
}
void fileopen(std::ifstream& f, const std::string &file,const unsigned int i){

	f.open(file);
	checkIFile(f,file);
	skiplines(f, i);
	return;
}

std::string remove_chars(std::string source, std::string chars){
	std::string s = source; // copy
	//std::cout << "s=" << s << std::endl;
	for(unsigned int i = 0; i < chars.size(); ++i){
		char c = chars[i];
		//std::cout << "ELEMENT->" << c << "<-" << std::endl;
		//std::cout << "s.size() = " << s.size() << std::endl;
		for(unsigned int j = 0; j < s.size(); ++j){
			//std::cout << s[j] << std::endl;
			if(s[j] == c){
				//std::cout << s << std::endl;
				s.erase(s.begin()+j);
				//std::cout << s << std::endl;
				j--;
			}
		}
		//std::cout << s << std::endl;	
	
	}
	//std::cout << s << std::endl;
	return s;
}

std::vector<std::string> line_to_vec(const std::string line, const char delimiter, char mode){
	switch(mode){

		case 's':{ 
			return line_to_vec(line, delimiter);
			break;
		}
		case 'u':{ // reduce problem to solveable.
			std::string s = line;
			//std::cout << "the string>>>" << s << "<<<" << std::endl;
			bool previous_is_char = false;
			for(size_t i = 0; i < s.size(); ++ i){
			//	std::cout << i << "_ _" << s.size() << "***" << s[i] << "***" << std::endl << s << std::endl;
				if(s[i] == delimiter){
			//		std::cout << " A " << std::endl;
					if(previous_is_char){
			//			std::cout << " B " << std::endl;
						s.erase(s.begin() + i);
			//			std::cout << " C " << std::endl;
						--i; // otherwise i would skip next	
					}					
					previous_is_char = true;

				}
				else
					previous_is_char = false;
			}
			//std::cout << " D " << std::endl;
			//check first and last element
			if(s.front() == delimiter){
			//	std::cout << " z " << std::endl;
				s.erase(s.begin() + 0);
			//	std::cout << " z2 " << std::endl;
			}
			if(s.back() == delimiter){
			//	std::cout << " w " << std::endl;
				s.erase(s.size()-1,1);
			//	std::cout << " w2 " << std::endl;
			}
			//std::cout << "the string>>>" << s << "<<<" << std::endl;
			//apply known function:
			const std::string sconst = s;
			return line_to_vec(sconst, delimiter);
		}
		default:{
			error("undefined mode in \"std::vector<std::string> line_to_vec(const std::string line, const char delimiter, const char mode)\". Abort.");
			exit(4);
		}

	}
}

bool IsPathExist(const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

/*void readQ2Grid(vector<double> & Q2binLowVec, vector<double> & Q2binHighVec, vector<double> & Q2binVec, const string &adress, const char &delimiter) //expects sorted file of Q2 bins
{
	#ifdef DEBUG
	std::cout << "in readQ2Grid" << std::endl;
	#endif 
	vector<double> Q2; 
	vector<double> Q2min;
	vector<double> Q2max;
	
	map<string,int> posInGrid;
	posInGrid["Q2"] = 3;
	posInGrid["Q2min"] = 0;
	posInGrid["Q2max"] = 2;
	posInGrid["xB"] = 7;
	posInGrid["xBmin"] = 4;
	posInGrid["xBmax"] = 5;
	std::ifstream f;
	string line;
	vector<string> stringVec;
	
	int posQ2 = posInGrid["Q2"]; //=map(Q2)
	int posQ2min = posInGrid["Q2min"];//=map(Q2min)
	int posQ2max = posInGrid["Q2max"]; //=map(Q2max)
	
	bool newValues = false;
	
	f.open(adress);
	if(not f.good()){
		std::cerr << "Failed to open \"" << adress << "\" for writing." << endl;
	}
	for(int i=0; getline(f,line);)
	{
		 
		
		line_to_vec(line, stringVec, delimiter);
		// expect all bins to be meaningful, no strange overlap whatsoever
		
		//the following looks complicated. My idea is to be work around the problem "program might want to check Q2[-1, 0] or similar while it does not exist.
		if(i==0)
		{
			newValues = true;
		}
		else
		{
			if(Q2[i-1]!=stod(stringVec[posQ2]) || Q2min[i-1]!=stod(stringVec[posQ2min]) || Q2max[i-1]!=stod(stringVec[posQ2max]))
			{	
				
				newValues = true;
			}
			else{ newValues = false;}
		}
		
		// if value is new..
		if(newValues)
		{
			Q2.push_back(stod(stringVec[posQ2]));
			Q2min.push_back(stod(stringVec[posQ2min]));
			Q2max.push_back(stod(stringVec[posQ2max]));
			
			++i;
			
		}
		
	}
	
	
	//write new valus into vectors:
	Q2binLowVec =  Q2min;
	Q2binHighVec = Q2max;
	Q2binVec = Q2;
	return;
}*/

void getline(std::ifstream &inc, std::string& line, int linecounter)
{
	if(not inc.good()){
		std::cerr << "Failed to read from stream." << std::endl;
	}
	std::string buffer;
	inc.seekg(std::ios::beg); // return to start
	for(int i = 1; i<=linecounter; ++i)
	{
		getline(inc,buffer);
//		std::cout << buffer << std::endl;
	}
	line = buffer;
	return;
}
/*
string getline(std::ifstream &inc, string& line, int linecounter)
{
	if(not inc.good()){
		std::cerr << "Failed to read from stream." << endl;
	}
	string buffer;
	inc.seekg(std::ios::beg); // return to start
	for(int i = 1; i<=linecounter; ++i)
	{
		getline(inc,buffer);
	}
	line = buffer;
	return line;
}
*/

std::vector<double> v_stod(const std::vector<std::string> & svec){
	std::vector<double> dvec;
	
	for(unsigned int i = 0; i < svec.size(); ++i){
		dvec.push_back(std::stod(svec[i]));
	}
	return dvec;
}

std::vector<int> v_stoi(const std::vector<std::string> & svec){
	std::vector<int> ivec;
	
	for(unsigned int i = 0; i < svec.size(); ++i){
		ivec.push_back(stoi(svec[i]));
	}
	return ivec;
}

void true_false_input(bool& value, std::string& question){
	//std::cout << question << endl;
	shell_text(question);
	shell_text("y is true, n is false");
//	std::cout << "y is true, n is false" << endl;
	char tempchar;
	std::cin >> tempchar;

	switch(tempchar){
		case 'y':
			value = true;
			shell_text("continue with \"true\"");
			break;
		case 'n':
			value = false;
			shell_text("continue with \"false\"");
			break;
		default:
			error("input was " + std::to_string(tempchar) + "\nplease try again.");
			std::exit(4);
			//true_false_input(value, question);
	}
}

void getBin(const int line, std::ifstream &in, double &varBin, double &varLow, double &varHigh)
{
	std::string str;
	getline(in, str, line);
	//cout << str << endl;
	size_t begin = 0, end = 0; // because of an warning->error (yes it was an error) initialized the values for begin and end.
	
	size_t pos = str.find(",,,", 0);
	for(size_t i = pos+3; i < str.size(); ++i)
	{
		
		if(str[i] != ' ') //then expect a[i] is beginn of first number
		{
			begin = i;
			break;
		}
	}
	
	for(size_t i = begin+1; i < str.size(); ++i)
	{
		
		if(str[i] == ' ') //then expect a[i] is end of first number
		{
			end = i-1;
			break;
		}
	}
	
	std::string a = str.substr(begin,end-begin+1);
	varBin = std::stod(a);
	
	pos = str.find('=', 0);
	
	for(size_t i = pos+1; i < str.size(); ++i)
	{
		//cout << str[i] << endl;
		if(str[i] != ' ') //then expect a[i] is beginn of first number
		{
			begin = i;
			break;
		}
	}
	//cout << "begin is " << begin << endl;
	for(size_t i = begin+1; i < str.size(); ++i)
	{
		//cout << str[i] << endl;
		if(str[i] == ' ') //then expect a[i] is end of first number
		{
			end = i-1;
			break;
		}
	}
	//cout << "end is " << end << endl;
	
	a = str.substr(begin,end-begin+1);
	varLow = std::stod(a);
	//cout << varMin << endl;
	
	pos = str.find("TO", 0);
	
	for(size_t i = pos+2; i < str.size(); ++i)
	{
		//cout << str[i] << endl;
		if(str[i] != ' ') //then expect a[i] is beginn of first number
		{
			begin = i;
			break;
		}
	}
	//cout << "begin is " << begin << endl;
	for(size_t i = begin+1; i < str.size(); ++i)
	{
		//cout << str[i] << endl;
		if(str[i] == ')') //then expect a[i] is end of first number
		{
			end = i-1;
			break;
		}
	}
	//cout << "end is " << end << endl;
	
	a = str.substr(begin,end-begin+1);
	//cout << "#" << a << "#" << endl;
	varHigh = std::stod(a);
	//cout << varMax << endl;
	
}

std::string TrueFalse(bool b){
	if(b) return "True";
	else return "False";
}

void round_double_vec(std::vector<double>& v, const int n){
	for(size_t i = 0; i < v.size(); ++i){
		double factor = std::pow(10.,n);
		v[i] = std::ceil(v[i] * factor) / factor; 
	}	
	return;
}

void vec_cpy_sd(std::vector<std::string> & s_vec, std::vector<double> & d_vec)
{
	d_vec.clear();
	for(unsigned int i = 0; i < s_vec.size(); ++i)
	{
		d_vec.push_back(std::stod(s_vec[i]));
	}
	return;
}

void checkOFile(std::ofstream & ofstr, const std:: string &adress){
	if(!ofstr.is_open()){
		//std::cout << "Could not open the ofstream with adress " << adress << " aborting process." << std::endl;
		//exit(1);
		error("Could not open the ofstream with adress " + adress + " aborting process.");
		std::exit(1);
	}
	return;
}


void checkIFile(std::ifstream & ifstr, const std::string &adress)
{
	if(!ifstr.is_open())
	{

		error("Could not open the ifstream with adress " + adress + " aborting process.");
		std::exit(1);
	}
	return;
}

void getline(std::ifstream & f, std::vector<std::string> & vec, const char delimiter)
{
	std::string s;
	getline(f,s);
	line_to_vec(s, vec, delimiter);
	return;
}
double computeBjorkenx(const double nu, const double Q2, const double M)
{
	double x = Q2/(2.*M*nu);
	return x;
}

double average1(double x0, double x1, std::vector<double> grid, std::vector<double> vals){
	if(x0 < grid.front() || x1 > grid.back()){
		//std::stringstream s;
		//s << "average1: x0-x1 = " << x0 << "-" << x1 << " while grid is " << grid;
		//error(s.str());
		std::string s;
		const char space = ' ';
		s += "average1: x0-x1 = " + std::to_string(x0) + "-" + std::to_string(x1) + " while grid is " + vec_to_line(grid,space);
		error(s);		
		exit(3);
	}
	
	if(grid.size() != vals.size() + 1){

		//s << "amount of values should be equal to amount of bins. \n binning:	" << grid << "\n values:	" << vals;
		//error(s.str());
		std::string s;
		const char space = ' ';
		s += "amount of values should be equal to amount of bins. \n binning:	" + vec_to_line(grid,space) + "\n values:	" + vec_to_line(vals,space);
		error(s);
		exit(3);
	}
	if(!isOrdered(grid,'<')){
		//std::stringstream s;
		//s << "average1: grid not ordered:" << grid;
		//error(s.str());
		std::string s;
		const char space = ' ';
		s += "average1: grid not ordered:" + vec_to_line(grid,space);
		error(s);
		exit(3);
	}
	
	unsigned iMin = getBinIndex(grid,x0);
	unsigned iMax = getBinIndex(grid,x1);
	//std::cout << iMin << "-" << iMax << std::endl;
	double weightedSum = 0;
	double totalVolume = 0;
	//first
	double measure = grid[iMin + 1] - x0;
	weightedSum += vals[iMin] * measure;
	totalVolume +=  measure;
	//std::cout << iMin << " " << measure << " " << vals[iMin] << std::endl;
	//all values in between.
	for(unsigned i = iMin+1; i < iMax; ++i){
		double measure = grid[i+1]-grid[i];
		weightedSum += vals[i] * measure;
		totalVolume += measure;
		//std::cout << i << " " << measure << " " << vals[i] << std::endl;
	}
	//last
	measure = x1 - grid[iMax];
	weightedSum += vals[iMax] * measure;
	totalVolume += measure;
	//std::cout << iMax << " " << measure << " " << vals[iMax] << " " << weightedSum << " "<< totalVolume << std::endl;

	//std::cout << weightedSum << " " << totalVolume << std::endl;
	//normalise
	weightedSum *= 1./(totalVolume);
	
	return weightedSum;	
}

long unsigned int first_line_without_signature(const std::string& file, const std::string& signature){
	std::ifstream f;
	f.open(file);
	checkIFile(f,file);
	std::string line;
	
	long unsigned int i = 0;
	while(std::getline(f,line)){
		++i;
		
		bool has_signature(line.substr(0,signature.size()) == signature);
		
		if(!has_signature){ // this is the first line without this signature.
			//std::cout << line << std::endl;
			return i;
		}
	}
	//std::cout << "End of file " << file << ". End of signature " << signature << " was not found." << std::endl;
	warning("End of file " + file + ". End of signature " + signature + " was not found.");
	return 0;
}

long unsigned int final_line_ctr(const std::string& file){
	long unsigned int i = 0;
	std::ifstream f;
	f.open(file);
	checkIFile(f,file);
	std::string s;
	while(getline(f,s)){
		//std::cout << s << std::endl;
		++i;
	}
	f.close();
	return i;
}
long unsigned int final_nonempty_line_ctr(const std::string& file){
	const long unsigned int lastline = final_line_ctr(file); // this already checks for the input file.
	long unsigned int i = lastline;
	
	std::ifstream f;
	f.open(file);
	std::string line;
	getline(f, line, i);
	
	while(line == ""){ // while line does not contain info
		//std::cout << i << ":::::::" << line << std::endl;
		--i;
		getline(f, line, i);
	}
	return i;
}

void skiplines(std::ifstream& f, const unsigned int n){ // skip a few lines.
	for(unsigned int i = 0; i < n; ++i){
		std::string s;
		getline(f, s);
	}
	return;
}

bool check_File_exists(const std::string & file)
{
	std::ifstream f;
	f.open(file);
	if(f.is_open())
		return true;
	else
		return false;
}

std::string vec_to_line(std::vector<std::string> &vec, const char delimiter){
	// not to have spaceing after the last element.
	std::string line =""; // emptystring.
	for(size_t i = 0; i < vec.size()-1; ++i){
		line += vec[i] + delimiter;
	}
	line += vec.back(); // add last element by hand
	
	return line;
}

std::string vec_to_line(const std::vector<double> &vec, const char delimiter){
	// not to have spaceing after the last element.
	std::string line =""; // emptystring.
	for(size_t i = 0; i < vec.size()-1; ++i){
		line += std::to_string(vec[i]) + delimiter;
	}
	line += std::to_string(vec.back()); // add last element by hand
	
	return line;
}

std::string vec_to_line(const std::vector<int> &vec, const char delimiter){
	std::vector<std::string> sv;
	for(auto x : vec){
		sv.push_back(std::to_string(x));
	}
	std::string line = vec_to_line(sv, delimiter);
	return line;
}

std::string vec_to_line_precision(std::vector<double> &vec, const char delimiter, int precision){
	// not to have spaceing after the last element.
	std::string line =""; // emptystring.
	std::stringstream stream;
	for(size_t i = 0; i < vec.size()-1; ++i){
		stream << std::setprecision(precision) << vec[i] << delimiter;
	}
	stream << std::setprecision(precision) << vec.back(); // add last element by hand
	line = stream.str();
	
	return line;
}

std::string vec_to_line_scientific(std::vector<double> &vec, const char delimiter, int precision){
	// not to have spaceing after the last element.
	std::string line =""; // emptystring.
	std::stringstream stream;
	for(size_t i = 0; i < vec.size()-1; ++i){
		stream << std::setprecision(precision) << std::scientific << vec[i] << delimiter;
	}
	stream << std::setprecision(precision) << vec.back(); // add last element by hand
	line = stream.str();
	
	return line;
}

std::vector<std::string> file_to_vec(const std::string& file){ //read a whole file and the lines will become entries of the vector. this probably only makes sense for small files..
	std::ifstream f;
	f.open(file);
	checkIFile(f,file);
	
	std::vector<std::string> v;
	std::string s;
	while(getline(f,s)){
		v.push_back(s);
	}
	return v;
}

void save_values_to_file(std::ofstream &f, std::vector<double> &v){
	//i have to decide with which char i want to separate my values. use tab.
	const char tab = '	';
	std::string s = vec_to_line(v, tab);
	//std::string s = vec_to_line_precision(v, tab, 50); // hard coded..😐️
	f << s << std::endl;
	return;
	
}

void append_values_to_file(const std::string & filename, std::vector<double> &v){
	std::ofstream f;
	f.open(filename, std::ios_base::app);
	save_values_to_file(f, v);
	f.close();
}

int modsys(std::string cmd){
	int i = system(cmd.c_str());
	return i;
}

void file_append_1to2(const std::string& file1, const std::string& file2){
	std::ofstream ofs;
	std::ifstream ifs;
	
	ifs.open(file1);
	ofs.open(file2,std::ios_base::app);
	checkIFile(ifs,file1);
	checkOFile(ofs,file2);
	std::string s;
	while(getline(ifs,s)){
		ofs << s << std::endl;
	}
	ifs.close();
	ofs.close();
}

std::string get_Parameter(std::string file, std::string parametercode){
	
	//std::cout << "trying to find \"" << parametercode << "\" in file \"" << file << "\"..." << std::endl;
	shell_text( "trying to find \"" + parametercode + "\" in file \"" + file + "\"...");
	bool parameter_found = false;
	std::ifstream ifs;
	ifs.open(file);
	checkIFile(ifs, file);
	
	//int foundline = -1;
	int linectr = 0;
	
	std::string s;
	
	//here i define the symbols i will use:
	//const char separator = '=';
	const char comment = '!';
	
	
	std::string line;
	while(getline(ifs,line)){
		++linectr;
		//std::cout << line << std::endl;
		if(line == "") continue; // emptyline
		//std::cout << "not empty" << std::endl;
		if(line[0] == comment) continue; // a comment line
		//std::cout << "not comment" << std::endl;
		//if(line.find(separator) == string::npos) continue;
		//std::cout << "not without seperator" << std::endl;
		std::string cmd = line.substr(0, parametercode.size());
		//std::cout << "the command of this line is: \"" << cmd << "\"" << std::endl;
		if(cmd != parametercode) continue;
		
		//to check:
		//std::cout << "|" << cmd << "|" << std::endl;
		std::string parameter = line.substr(cmd.size() + 0);
		//std::cout << "|" << parameter << "|" << std::endl;
		if(!parameter_found){
			parameter_found = true;
			//foundline = linectr;
			s = parameter;
		}
		else{
			warning("multiple definition of the parameter \"" + parametercode + "\": \"" + parameter + "\". Already set and used value is \"" + s + "\".");
		}
	}
	ifs.close();
	//std::cout << "left loop over input file" << std::endl;
	
	if(!parameter_found){
		//std::cout << "The parameter \"" << parametercode << "\" was not found in file \"" << file << "\"" << std::endl;
		//exit(1);
		error("The parameter \"" + parametercode + "\" was not found in file \"" + file + "\"");
		std::exit(1);
	}
	
	//std::cout << "get_Parameter() routine finito... with value: " << s << std::endl;
	shell_text("get_Parameter() routine finito... with value: " + s);
	return s;
}

void print_dif_element(double a, double b){
	double dif = b - a;
	if(dif == 0){
	    //continue;// ignore those without change..
	    std::cout << "\033[0m";
	   }
	else if(dif > 0)
		std::cout << "\033[1;34m";
	else if(dif < 0)
		std::cout << "\033[1;33m";
	
	std::cout << b;
	std::cout << "\033[0m";
}

std::vector<double> difference(std::vector<double> &v1, std::vector<double>& v2){
	if(v1.size() != v2.size()){
		error("Called to compute difference for two vectors not of the same size:");
		std::cout << v1 << std::endl << v2 << std::endl;
		exit(3);
	}
	std::vector<double> diff;
	for(unsigned int i = 0; i < v1.size(); ++i){
		diff.push_back(v1[i]-v2[i]);
	}
	return diff;
}

std::vector<double> readbins(std::string file){
	std::vector<double> v;
	read_1dFile_to_vec(file, v);
	return v;
}

void print_dif(std::vector<double> &ov, std::vector<double> &nv){
	
	
	
	
	for(unsigned int i = 0; i < ov.size(); ++i){
//		double o, n, dif; //old, new, diff
		double o, n; //old, new
		o = ov[i];
		n = nv[i];
		if(i != 0){
			std::cout << '	';
		}
		print_dif_element(o,n);
	}
	/*std::cout << std::endl;
	for(unsigned int i = 0; i < ov.size(); ++i){
		double o, n, dif; //old, new, diff
		o = ov[i];
		n = nv[i];
		if(i != 0){
			std::cout << '	';
		}
		std::cout << dif;
		
		
	}*/
	
	//old version
	/*
	std::cout << "1st	2nd	dif" << std::endl;	
	
	
	for(unsigned int i = 0; i < ov.size(); ++i){
		double o, n, dif; //old, new, diff
		o = ov[i];
		n = nv[i];
		dif = n - o;
		
		if(dif == 0){
		    continue;// ignore those without change..
		    std::cout << "\033[0m";
		   }
		else if(dif > 0)
			std::cout << "\033[1;34m";
		else if(dif < 0)
			std::cout << "\033[1;33m";
		
		std::cout << o << "	" << n << "	" << dif << std::endl;
	}
	
	
	*/
    std::cout << "\033[0m";
    std::cout << std::endl;
    return;
}

/*std::string get_Parameter(std::string file, std::string parametercode){
	
	std::cout << "trying to find \"" << parametercode << "\" in file \"" << file << "\"..." << std::endl;	
	bool parameter_found = false;
	std::ifstream ifs;
	ifs.open(file);
	checkIFile(ifs, file);
	
	std::string s;
	
	//here i define the symbols i will use:
	const char separator = '=';
	const char comment = '!';
	
	
	std::string line;
	while(getline(ifs,line)){
		//std::cout << line << std::endl;
		if(line == "") continue; // emptyline
		//std::cout << "not empty" << std::endl;
		if(line[0] == comment) continue; // a comment line
		//std::cout << "not comment" << std::endl;
		if(line.find(separator) == string::npos) continue;
		//std::cout << "not without seperator" << std::endl;
		std::string cmd = line.substr(0, line.find(separator));
		//std::cout << "the command of this line is: \"" << cmd << "\"" << std::endl;
		if(cmd != parametercode) continue;
		
		//to check:
		std::cout << "|" << cmd << "|" << std::endl;
		std::string parameter = line.substr(line.find(separator) + 1);
		std::cout << "|" << parameter << "|" << std::endl;
		
		parameter_found = true;
		s = parameter;
		break;
	}
	
	//std::cout << "left loop over input file" << std::endl;
	
	if(!parameter_found){
		std::cout << "The parameter \"" << parametercode << "\" was not found in file \"" << file << "\"" << std::endl;
		exit(1);
	}
	
	
	ifs.close();
	//std::cout << "get_Parameter() routine finito..." << std::endl;
	return s;
}*/

double get_Parameter_double(std::string file, std::string parametercode){
	//std::cout << "find a double" << std::endl;
	std::string s = get_Parameter(file, parametercode);
	//std::cout << s << std::endl;
	//double x = stod(s);
	//std::cout << x << std::endl;
	return std::stod(s);
}

int get_Parameter_int(std::string file, std::string parametercode){
	std::string s = get_Parameter(file, parametercode);
//	std::cout << s << std::endl;
	return stoi(s);
}

char get_Parameter_char(std::string file, std::string parametercode){
	std::string s = get_Parameter(file, parametercode);
	return s[0];
}

bool get_Parameter_bool(std::string file, std::string parametercode){
	int i = get_Parameter_int(file, parametercode);
	
	bool b = bool(i);
	//std::cout << "value of parameter \"" << parametercode << "\":" << b << i << std::endl;
	return b;
}

bool getline_as_vec(std::ifstream &ifs, std::vector<std::string> &line_vec, const char delimiter){
	std::string line;
	if(getline(ifs,line)){
	
		
		line_to_vec(line, line_vec, delimiter);
		return true;
	}
	else{
		return false;
	}

}

void getline_asVec_withEntry_atPosition(std::ifstream & f, std::vector<std::string> & vec, const char delimiterSymbol, std::string & entry, const unsigned int i)
{
	std::string line;
	std::vector<std::string> line_vec;
	while(getline(f, line))
	{
		line_to_vec(line, line_vec, delimiterSymbol);
		if(line_vec[i] == entry)
		{//got it
			vec = line_vec;
			return;
		}
		
	}
	//std::cout << "From: getline_asVec_withEntry_atPosition: \nreached the end of the document( out of scope of \"while-getline\") and did not find a line that matches the condition:\ncondition = the entry \"" << entry << "\" at position " << i << std::endl;
	warning("From: getline_asVec_withEntry_atPosition: \nreached the end of the document( out of scope of \"while-getline\") and did not find a line that matches the condition:\ncondition = the entry \"" + entry + "\" at position " + std::to_string(i));
	return;
}

void progress2terminal(const unsigned int ctr, const unsigned int total, const double ratio){// print to terminal the progress, if you have more than certain % of total events progressed
	unsigned int x = (ratio*double(total));
	if( ctr % x == 0 )
		std::cout << (double(ctr)/  total) * 100. << "%" <<  " done: " << ctr << " of " << total << " steps done." <<  std::endl;
		std::string txt = std::to_string( (double(ctr)/  total) * 100.) + "%" + " done: " + std::to_string(ctr) + " of " + std::to_string(total) + " steps done.";
		shell_text(txt);
	return;
}
void extract2_1(std::string & file, std::string &newfile, std::vector<int> &positions, std::vector<double> & values, const char delimiter)
{
	if(!check_File_exists(file))
	{
		//std::cout << "file \"" << file << "\" does not exist" << std::endl;
		warning("file \"" + file + "\" does not exist");
		return;
	}
	std::ifstream in; 
	in.open(file);
	std::ofstream off;
	off.open(newfile);
	
	std::string line;
	std::vector<std::string> line_vec;
	while(getline(in, line))
	{
		line_to_vec(line, line_vec, delimiter);
		bool line_matches = true;
		for(unsigned int i = 0; i < positions.size(); ++i)
		{
			if(std::stod(line_vec[positions[i]]) != values[i])
			line_matches = false;
		}

		if(line_matches)
		off << line << std::endl;		
	}
	
	in.close();
	off.close();
}

/*std::string scientific(double d, int dec){
	if(d == 0) {
		std::string s("0",dec+1);
		s.insert(1,".");
		s.append("e+00");
		return s;
	}
	int errorcode = 0;
	std::string s,as;
	as = std::to_string(d);
	s = as;
	std::cout << "###################" << s << std::endl;
	bool problem = false;
	char sgn;
	if(s[0] == '-') sgn = '-';
	else sgn = '+';
	
	std::cout << "sign is: "<< sgn << std::endl;
	
	if(sgn == '-'){ 
		s.erase(0);
	}	
	std::cout << " s now is " << s << std::endl;
	
	size_t dot = s.find('.');
	if(dot == string::npos) {errorcode = 1; problem = true;}
	s.erase(dot,1);
	std::cout << s << std::endl;
	bool positive_10_exp;
	size_t dot_shift;
	std::cout << "PPPP" << s << std::endl;
	if(s[0] == '0'){
		std::cout << s << std::endl;
		positive_10_exp = false;
		for(unsigned i = 1; i < s.size(); ++i){
			if(s[i] != '0'){
				dot_shift = i;
				std::cout << "***" << i << " xxx " << s[i] << std::endl;
				break;
			}
		}
		std::cout << s << std::endl;
		errorcode = 2; problem = true;
	}
	else{
		positive_10_exp = true;
	}
	
	std::string exponent = "e";
	if(positive_10_exp){
		exponent += "+" + std::to_string(dot - 1);
	}
	else{
		exponent += "-" + std::to_string(dot_shift + 1);
	}
	
	s.insert(1,".");
	
	s = s.substr(0, 1 + 1 + dec); // first number, dot, and then decimals
	//now add potenz.
	s += exponent;
	
	
	if(problem) {
		std::cout << "need help with double " << d << " and #decimals " << dec << ". I cannot do this on my own." << std::endl;
		std::cout << errorcode << std::endl;
		//exit(1);
	}
	return s;
}*/

std::vector<double> get_Parameter_bin(std::string file, std::string parametercode, std::string separator){ // get a pair of parameters separated by some string. this is bound to double values.
	std::string s = get_Parameter(file, parametercode);
	//std::cout << s << std::endl;
	std::vector<std::string> v = stov(s, separator);
	
	std::vector<double> vd = v_stod(v);
	
	return vd;
}

std::string get_Parameter_L(const std::string& line, std::string parametercode){
	if(line.substr(0, parametercode.size()) == parametercode){
		return line.substr(parametercode.size());
	}
	else{warning(line + " does not start with \"" + parametercode + "\", returning empty string as parameter."); return "";}
}

std::vector<std::string> get_Parameter_bin_L(const std::string& line, std::string parametercode, std::string delimiter){
	std::string s = get_Parameter_L(line, parametercode);
	std::vector<std::string> sv;
	line_to_vec(s, sv, delimiter);
	return sv;
}

std::vector<std::string> stov(const std::string& s, const std::string& delimiter){
	std::vector<std::string> sv;
	line_to_vec(s, sv, delimiter);
	return sv;
}

std::vector<unsigned int> polyadisch(const unsigned int basis, const unsigned int number){
	int max_expo = 1;
	while(pow(basis,max_expo) < number){
		 ++max_expo;
	}
	//std::cout << "number is " <<  number << " basis is " << basis << " first overshooting expo is" << max_expo << std::endl;
	//now j is highest expo of our candidate i to basis "basis"
	//this is now what i want to have
	if(pow(basis,max_expo) < number){
		//std::cout << "problem in checking highest exponent" << std::endl; exit(1);
		error("problem in checking highest exponent");
		std::exit(2);
	}
	
	std::vector<unsigned int> nums;
	unsigned int rest = number;
	for(int j = max_expo-1; j >= 0; --j){
		int b = 0; // i ONLY define this variable here because a different programm complained about this during compilation. 
		//std::cout << "rest is " << rest << std::endl;
		for(unsigned int x = 0; x < basis; ++x){
			unsigned int val = (x+1) * int(pow(basis,j));
			//cout << val << " " << j << endl;
			if(val > rest){
				b = x;
				//cout << "break with coeff: " << b << endl;
				break;
			}
		}
		rest -= b * int(pow(basis,j));
		
		nums.push_back(b);
		
	}
	//print1dvec(nums);
	return nums;
}

std::vector<std::vector<double>> file_to_matrix(const std::string &file, const char delimiter){ // read a whole file ( THAT MUST ONLY CONSIST OF NUMBERS!! ) to a matrix: vec<vec<double>>.
	std::vector<std::string> vec_file = file_to_vec(file);
	std::vector<std::vector<double>> matrix; // acces matrix elements through matrix[i][j], therefore a good convention is: i are lines, j are coloumns.
	
	//const char delimiter = '	'; // TAB as default.. # this is updated after creating the original code. Thus it is highly likely that there are complications with older version. Just edit code by also passing delimiter when calling the function
	
	for(auto l : vec_file){ // iterate over lines of file (already in vector stored)
		
		std::vector<std::string> vs = line_to_vec(l, delimiter, 's'); // s is single symbol. Maybe will be adapted some day.
		std::vector<double> vd = v_stod(vs); // convert to double type vector. here it may fail if there is somewhere a non convertable character..
		
		matrix.push_back(vd);
		
	}
	
	return matrix;
	
}

std::vector<double> matrix_col(std::vector<std::vector<double>> matrix, const unsigned int col){
	std::vector<double> vec;
	unsigned int l = matrix.size();
	for(unsigned int i = 0; i < l; ++i){
		vec.push_back(matrix[i][col]);
	}
	return vec;
}


std::vector<std::vector<std::string>> file_to_stringmatrix(const std::string &file){ // read a whole file ( THAT MUST ONLY CONSIST OF NUMBERS!! ) to a matrix: vec<vec<double>>.
	std::vector<std::string> vec_file = file_to_vec(file);
	std::vector<std::vector<std::string>> matrix; // acces matrix elements through matrix[i][j], therefore a good convention is: i are lines, j are coloumns.
	
	const char delimiter = '	'; // TAB as default..
	
	for(auto l : vec_file){ // iterate over lines of file (already in vector stored)
		
		std::vector<std::string> vs = line_to_vec(l, delimiter, 's'); // s is single symbol. Maybe will be adapted some day.
		
		matrix.push_back(vs);
		
	}
	
	return matrix;
	
}



void cut_file(const std::string &src, const std::string &target, unsigned int firstline, unsigned int lastline){

	std::ofstream ofs;
	std::ifstream ifs;
	
	ifs.open(src);
	checkIFile(ifs,src);
	
	ofs.open(target);
	checkOFile(ofs,target);
	
	for(unsigned int i = firstline; i <= lastline; ++i){
		std::string line;
		getline(ifs, line, i);
		ofs << line << std::endl;	
		
	}
	
}

void cprint(std::string s, char colour){ // print with colour that is mapped by char.
	std::string colourout;
	std::string back_to_normal = "\033[0m";
	switch(colour){
		case 'k':
		colourout = "\033[1;30mbold";
		break;
		
		default:
		colourout = "\033[1;35mbold";
		break;
	}
	
	std::cout << colourout << s << back_to_normal;
	return; 
}

long unsigned int polyadisch_inverse(const unsigned int basis, const std::vector<unsigned int>& representation){
	long unsigned int value = 0;
	for(unsigned int i = 0; i < representation.size(); ++i){
		value += representation[representation.size() - (i + 1)] * pow(basis,i);
		//cout << value << endl;
	}
	return value;
}

std::vector<double> logGrid(double minVal, double maxVal, double basis, int steps, char mode) // create log grid with range min-maxVal, have exponential basis "basis" and have steps, total (mode = 't') or relative (mode = 'r') to the basis step.
{
	std::vector<double> grid;
	switch(mode){
	case 't':{
		for(int i = 0; i <= steps; ++i)
		{	
			grid.push_back(pow(basis,log10(minVal)/log10(basis)+double(i)/steps * log10(maxVal/minVal)/log10(basis)));
		}
		break;}
		
	case 'r':{
		// compute #total steps resulting and call the function again with that argument. genius.
		int decadesCovered = int(log10(maxVal/minVal)/log10(basis)); 
		int s = decadesCovered  * steps ; // <- thats #decades-covered  x stepsPerPower = steps in total
		
		grid = logGrid(minVal, maxVal, basis, s, 't');
		break;}
		
	default:{
		std::cout << "std::vector<double> logGrid(double minVal, double maxVal, double basis, int steps, char mode) -- Did not recognize the char \"mode\" provided: " << mode << std::endl << "	expected options for \"mode\" are: \n't'	for total number of steps provided in \"steps\"\n'r' for relative number of steps per power provided in \"steps\"\n" << std::endl;
		warning("std::vector<double> logGrid(double minVal, double maxVal, double basis, int steps, char mode) -- Did not recognize the char \"mode\" provided: " + std::to_string(mode) + "\n	expected options for \"mode\" are: \n't'	for total number of steps provided in \"steps\"\n'r' for relative number of steps per power provided in \"steps\"\n");
		break;}
	}
	/*
	
for(int i = 0 ; i < sizeOfArrays; ++i)
	{
		range[i]=pow(10.,log10(lowBorder)+double(i)/binsPerDecade);
		unnormalizedvalue[i] = 0.;
	}
*/
	return grid;
}

std::vector<double> grid_around_vector(std::vector<double>& v){
	std::vector<double> grid(v.size() + 1, 0.);
	for(unsigned i = 0; i < v.size()-1; ++i){
		grid[i+1] = (v[i]+v[i+1])/2.;
	}
	grid[0] = 0.5 * v.front();
	grid[v.size()] = 2.*v.back();
	return grid;
}

std::vector<double> linGrid(double minVal, double maxVal, int steps) // create linear grid with range min-maxVal and number of steps. 1 steps will be one interval
{
	std::vector<double> grid;
	double delta = (maxVal-minVal)/steps;
	//grid.push_back(minVal);
	for(int i = 0; i <= steps; ++i)
	{
		grid.push_back(minVal + i*delta);
	}
	return grid;
	
}

std::vector<double> linGrid(double minVal, double maxVal, double steps_width){
	std::vector<double> grid;
	if(minVal > maxVal) error("linear Grid cannot be created with the max value smaller than the min value: (double minVal, double maxVal, double steps_width) = (" + std::to_string(minVal) + ", " + std::to_string(maxVal) + ", " + std::to_string(steps_width) + ")");
	//grid.push_back(minVal);
	double v = minVal;
	while( v < maxVal + steps_width/2. ){
		grid.push_back(v);
		v += steps_width;
	}
	return grid;
	
}

std::vector<int> load_distribution_fromfile_int(std::string f){
	std::vector<int> r;
	std::vector<std::string> file = file_to_vec(f);
	for( std::string s : file){
		
		int v = std::stoi(s);
		r.push_back(v);
	}
	
	return r;
}

std::vector<double> load_distribution_fromfile_double(std::string f){
	std::vector<double> r;
	std::vector<std::string> file = file_to_vec(f);
	for( std::string s : file){
		
		double v = std::stod(s);
		r.push_back(v);
	}
	
	return r;
}

std::vector<double> empty_vector(unsigned int size) // incredibly stupid because there is a standard function for this.
{	
	std::vector<double> v;
	for(unsigned int i = 0; i < size; ++i){
		v.push_back(0.); // 0 should be then ( i think ) automatically interpreted as the correct data type.
	}
	return v;
}

unsigned int getIndex(std::vector<std::string> &vec, std::string entry)
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
			}
			else
			{
			
				//cout << "the vector seems to contain at least two times the entry \"" << entry << "\" at the positions " << position << " and " << i << ".\nAbort." << endl;
				//std::exit(2); better leave it as a warning but keep running:
				//std::cout << "the vector seems to contain at least two times the entry \"" << entry << "\" at the positions " << position << " and " << i << ".\nThe entry: " << entry << std::endl << "The vetcor" <<  std::endl;
				warning("the vector seems to contain at least two times the entry \"" + entry + "\" at the positions " + std::to_string(position) + " and " + std::to_string(i) + ".\nThe entry: " + entry + "\nThe vetcor");
				print1dvec(vec);
				
			}
		}	
	}
	if(already_found)
	return position;
	else
	{
		error("the vector does not seem to contain the entry \"" + entry + "\". ");
//		std::cout << "the vector does not seem to contain the entry \"" << entry << "\". " << endl;
		std::cout << "the vector is: " << std::endl;
		std::cout << vec << std::endl;
		//error("abort");
		std::exit(3);
	}
}
