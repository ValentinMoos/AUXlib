#ifndef INPUTFILEREADER
#define INPUTFILEREADER

#include "Aux.h"
//#include <vector>

/*idea of this class:
	have an object that reads the variables you want from a certain file. Since the file will be often the same one for multiple readouts, it makes sense to save it in a class object.


*/
class InFiR{

	private:
		//member
		//const char delimiter = '	'; //delimiter, TAB
		
		std::string src_file;
		
		char default_data_type = 's'; //default: as string, available types are string s, char c, integer i, double d, float f, boolean b
		
		std::string commentblock = "_!"; // if some line should be ignored fully
		std::string initials = ""; // the first few chars in a line that are irrelevant
		std::string separation = ""; // seperation between name of parameter and value of parameter
		
		std::vector<std::string> parameter_codes;
		std::vector<std::string> parameter_values;
		std::vector<std::string> non_accountables;
		
		std::vector<std::string> bin_vars;
		std::vector<std::vector<std::string>> bin_vals;
		
		size_t spsz = 0;                                    // dimension = # of variables, size of "values"
		//function	
		
		unsigned int lines_to_skip_01 = 0;
		
		std::string sep_var_bin = "";
		std::string sep_bin_itself = "";
		
		std::string current_parameter = "";
		void get_parameter();
		
	public:
		
		InFiR(); //default constructor
		InFiR(const std::string sourcefile);               // declare sourcefile only
		InFiR(const std::string sourcefile, const std::string _commentblock);
		InFiR(const std::string sourcefile, const std::string _commentblock, const std::string _separation);
		
		void read_file(); // read the file, if no arg is given read the src_file, if an arg is given read the file added.
		void read_file(std::string);
		
		
		
		//void set_source(const std::string&);
		void set_source(const std::string);
		
		//void set_default_returntype(const char);
		
		void set_separation(const std::string);
		
		void set_initials(const std::string);
		//void set_initials(std::string);
		
		void set_commentblock(const string);
		//void set_commentblock(const string&);
		
		void set_separator_bin(const std::string s);
		
		void set_separator_variable_bin(const std::string s);
		
		void info();                                       //print info(member data) to terminal
		
		void t_print(); // print whole info to terminal
		
		std::string read_param_as_string(std::string);
	/*	
		template <typename T>
T read(const std::string keyword){
	std::string s = read_param_as_string(keyword);
	T val;
	//bool type_matched = false;
	if(typeid(val) == typeid(double)){
		val = stod(s);
	}
	else if(typeid(val) == typeid(int)){
		val = stoi(s);
	}
	else if(typeid(val) == typeid(char)){
		val = s[0];
	}
	else if(typeid(val) == typeid(bool)){
		val = true ? (s == "True", s == "true", s == "yes", s == "wahr", s == "Wahr", s == "1") : false;
	}
	else if(typeid(T) == typeid(std::string)){
		val = s;
	}
	else if(typeid(T) == typeid(float)){
		val = float(stod(s));
	}
	else
		std::cout << "Warning::from Inputfile::Reader.h::the type could not be matched." << std::endl;
	
	return val;
	
}*/
	double read_double(std::string);
	int read_int(std::string);
	float read_float(std::string);
	char read_char(std::string);
	bool read_bool(std::string); // true is one of the following strings: {True, true, yes, wahr, Wahr, 1}. Everything else is considered to be false.
		
};



namespace InfRead{ }
#endif
