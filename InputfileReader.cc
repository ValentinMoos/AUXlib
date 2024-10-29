#include"InputfileReader.h"
//Constructors
InFiR::InFiR(){} //uses default variables which are already set in header

InFiR::InFiR(const std::string sourcefile){
	
	set_source(sourcefile);
	
}

InFiR::InFiR(const std::string sourcefile, const std::string _commentblock) : InFiR(sourcefile){
	
	set_commentblock(_commentblock);
}

InFiR::InFiR(const std::string sourcefile, const std::string _commentblock, const std::string _separation) : InFiR(sourcefile, _commentblock){
	
	set_separation(_separation);
}

/*InFiR::InFiR(const std::string s1, const std::string s2, const std::string s3 sourcefile, const std::string _commentblock, const std::string _separation) : InFiR(sourcefile, _commentblock){
	
	set_separation(_separation);
}*/



// set - functions
/*void InFiR::set_source(const std::string s){
	const std::string &s_ref = s;
	set_source(s_ref);
	return;
}*/
void InFiR::set_source(const std::string s){
	if(!check_File_exists(s)){
		std::cout << "the file " << s << " does not seem to be existent." << std::endl;
		exit(1);
	}
	src_file = s;
	
	return;
}

void InFiR::set_separator_bin(const std::string s){
	sep_bin_itself = s;
	
	return;
}

void InFiR::set_separator_variable_bin(const std::string s){
	sep_var_bin = s;
	
	return;
}

/*void InFiR::set_default_returntype(const char returntype){
	
	if(returntype == 'c' || returntype == 'i' || returntype == 'd' || returntype == 's' || returntype == 'f' || returntype == 'b'){
		default_data_type = returntype;
	}
	else{
		std::cout << "the returntype \"" << returntype << "\" is unknown." << std::endl;
		exit(1);
	}
	return;
	
}*/

void InFiR::set_commentblock(const std::string _commentblock){
	commentblock = _commentblock;
	
	return;
}

void InFiR::set_initials(const std::string s){
	initials = s;
	
	return;
}

/*void InFiR::set_initials(std::string s){

	const std::string aux_s = s;
	
	set_initials(aux_s);
	return;
}*/

void InFiR::set_separation(const std::string _separation){
	separation = _separation;

	spsz = separation.size();	
	return;
}



//readers

/*InFiR::get_parameter(){
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
		if(line.substr(0, cbs) == commentblock) continue; // a comment line
		//std::cout << "not comment" << std::endl;
		if(line.find(separator) == string::npos){ std::cout << "This should not occurr: found a non comment line without declaration" continue;}
		//std::cout << "not without seperator" << std::endl;
		
		size_t para_name_begin = line.find(parametercode);
		if(para_name_begin == string::npos)continue; // the parameter is not to be found in this line
		
		//now we are sure that we found the parametercode in this line. Check if after the parameter there is the separation code
		
		std::string rest_of_line = line.substr(para_name_begin + parametercode.size());
		if(rest_of_line.find(separation))
		
		if()
		
		
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

double InFiR::read_double(std::string keyword){
	std::string s = read_param_as_string(keyword);
	
	return stod(s);
}

int InFiR::read_int(std::string keyword){
	std::string s = read_param_as_string(keyword);
	
	return stoi(s);
}

char InFiR::read_char(std::string keyword){
	std::string s = read_param_as_string(keyword);
	
	return s[0];
}

float InFiR::read_float(std::string keyword){
	std::string s = read_param_as_string(keyword);
	
	return stof(s);
}

bool InFiR::read_bool(std::string keyword){
	std::string s = read_param_as_string(keyword);
	bool b = true ? (s == "True", s == "true", s == "yes", s == "wahr", s == "Wahr", s == "1") : false;
	return b;
}


void InFiR::read_file(std::string file){
	//std::string src_backup = src_file;
	
	src_file = file; // make a copy because this function should leave the src_file parameter intact.
	read_file();
	//set_source(src_backup);	// restore src_file parameter.
	return;
}



void InFiR::read_file(){

	std::ifstream f;
	
	f.open(src_file);
	checkIFile(f, src_file);
	std::string line;
	
	std::cout << "Reading from file " << src_file << std::endl;
	
	unsigned int info_lines_ctr = 0;	
	
	// heart code of this class
	
	std::vector<std::string> info_depot;
	
	//start at top of line and read until standard beginning of line ends.
	while(getline(f,line)){
		info_lines_ctr++; // raise counter
		//std::cout << line << std::endl;
		if(line == "") continue; // emptyline, skip it
		//std::cout << "not empty" << std::endl;
		std::string comblock_substring = line.substr(0, commentblock.size());
		//std::cout << comblock_substring << "|" << commentblock << std::endl;
		if(comblock_substring == commentblock) continue; // a comment line
		//std::cout << "not comment" << std::endl;
		
		if(! (line.substr(0,initials.size()) == initials)){
			std::cout << "understand this line as the end of header-info." << std::endl;
			break; // if the initials are not found at the beginning, and its not a comment and its not empty we break this.
		}
		
		//std::cout << "add this line" << std::endl;
		
		info_depot.push_back(line); // if we are here: the line is not empty, its not a comment, and its contains the initials at the beginning.
		/*
		if(line.find(separator) == string::npos){ std::cout << "This should not occurr: found a non comment line without declaration" continue;}
		//std::cout << "not without seperator" << std::endl;
		
		size_t para_name_begin = line.find(parametercode);
		if(para_name_begin == string::npos)continue; // the parameter is not to be found in this line
		
		//now we are sure that we found the parametercode in this line. Check if after the parameter there is the separation code
		
		std::string rest_of_line = line.substr(para_name_begin + parametercode.size());
		if(rest_of_line.find(separation))
		
		if()
		
		
		//std::cout << "the command of this line is: \"" << cmd << "\"" << std::endl;
		if(cmd != parametercode) continue;
		
		//to check:
		std::cout << "|" << cmd << "|" << std::endl;
		std::string parameter = line.substr(line.find(separator) + 1);
		std::cout << "|" << parameter << "|" << std::endl;
		
		parameter_found = true;
		s = parameter;
		break;*/
	}
	f.close();
	//maybe read a certain amount of extra lines?
	//if so: also increase counter
	
	
	//tell me how many lines there are comments and info lines until the data is given
	
	std::cout << "lines added: " << info_depot.size() << std::endl;
	
	for(size_t i = 0; i < info_depot.size(); ++i){
		std::cout << info_depot[i] << std::endl;
		std::string s = info_depot[i].substr(initials.size());
		std::cout << s << std::endl;
		std::string keyword, value;
		size_t sepos = s.find(separation);
		std::cout << sepos << std::endl;
		if(sepos == string::npos){
			std::cout << "non accountable" << std::endl;
			non_accountables.push_back(s);
			std::cout << s << std::endl;
			continue;
		}
		keyword = s.substr(0,sepos);
		value   = s.substr(sepos + separation.size()); // till end, without the separation itself
		std::cout << "keyword: " << keyword << " | value: " << value << std::endl;
		parameter_codes.push_back(keyword);
		parameter_values.push_back(value);
		
	}
	info_depot = non_accountables; // everthing is now either in "keyword + value" or in "non accountable.
	non_accountables.clear();
	//now for bins

	for(size_t i = 0; i < info_depot.size(); ++i){
		// does it have the signature?
		std::string s = info_depot[i];
		std::cout << s << std::endl;
		size_t p1 = s.find(sep_var_bin);
		size_t p2 = s.find(sep_bin_itself);
		
		if(p1 == string::npos || p2 == string::npos){
			non_accountables.push_back(s);
			//std::cout << s << std::endl;
			continue;
		}
		else{
			std::cout << p1 << "|" << p2 << std::endl;
			std::string var = s.substr(0,p1);
			std::string bin1 = s.substr(p1 + sep_var_bin.size(), p2 - (p1 + sep_var_bin.size()) );
			std::string bin2 = s.substr(p2 + sep_bin_itself.size());
			
			std::cout << bin1 << std::endl << bin2 << std::endl;
			
			std::vector<std::string> bin {bin1, bin2};
			
			bin_vars.push_back(var);
			bin_vals.push_back(bin);
			
		}
	}
	
	
	return;
}

std::string InFiR::read_param_as_string(std::string param){
	for(size_t i=0; i < parameter_codes.size(); ++i){
		if(parameter_codes[i] == param){
			return parameter_values[i];		
		}
	}
	std::string warning = "The parameter \"" + param + "\" was not found.";
	std::cout << warning << std::endl;
	return warning;
}

// other member functions
void InFiR::info(){
	std::cout << "InFiR has the data:" << std::endl << "src_file	 = 	" << src_file << std::endl << "commentblock	=	" << commentblock << std::endl << "default data return type	=	" << default_data_type << std::endl;
	return;
}

void InFiR::t_print(){
	std::cout << "The following information is stored: " << std::endl << "parameter	" << "|-|" << "	value" << std::endl;
	
	for(size_t i=0; i < parameter_codes.size(); ++i){
		std::cout << parameter_codes[i] << "	" << "|-|" << "	" << parameter_values << std::endl;
	}


	std::cout << "The following information is stored for bins" << std::endl << "variable	" << "|-|" << "	min_val	" << "-" << "	max_val" << std::endl;
	
	for(size_t i=0; i < bin_vars.size(); ++i){
		std::cout << bin_vars[i] << "	" << "|-|" << "	" << bin_vals[0][i] << "	-	" << bin_vals[1][i] << std::endl;
	}
	
	std::cout << "The following lines were saved but not stored as a parametervalue: " << std::endl;
	
	for(size_t i=0; i < non_accountables.size(); ++i){
		std::cout << non_accountables[i] << std::endl;
	}
	
	return;
}



