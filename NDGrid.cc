#include"NDGrid.h"
//CONSTRUCTORS
NDGrid::NDGrid(){} //default does not do anything

NDGrid::NDGrid(std::vector<std::vector<double>> bins){
	//initialize all the members. 
	
	//values = bins; // the whole grid info.
	
	//reSetGrid(); // the rest is done here
	
	newGrid(bins); // (re) sets the bins - does no matter if its the first initialization or it overwrites all info.
}

NDGrid::NDGrid(std::string file){
	std::vector<std::vector<double>> m = file_to_matrix(file,'	'); // default: use tab.
	std::vector<std::vector<double>> grid;
	const unsigned int lines = m.size(); // the lines
	const unsigned int cols = m[0].size();
	
	unsigned int step = 1;
	const unsigned int vars = (cols-1) / 2; // = dimension of NDGrid
	//std::cout << "nvariables " << vars << std::endl;
	
	for(unsigned int i = 0; i < vars; ++i){
		//std::cout << i << std::endl;
		std::vector<double> bin_i {m[0][2*vars - 2*(i+1)], m[0][2*vars - 2*(i+1)+1]}; // get first two values in binning for this variable.
		
		unsigned int j = 1;
		bool variable_done = false;
		while(!variable_done){
			unsigned int line = j*step;
			if(line >= lines){
				//its done. this means i will assume i did the last bin and can break;
				break;
			}
			double element = m[line][2*vars - 2*(i+1)+1];
			//std::cout << element << std::endl;
			//std::cout << bin_i << std::endl;
			if(element == bin_i[1]){
				//std::cout << " variable done." << std::endl;
				variable_done = true;
				//std::cout << bin_i << std::endl;
			}
			else{
				//std::cout << " not done yet." << std::endl;
				bin_i.push_back(element);
				++j;
			}
			//std::cout << bin_i << std::endl;
		
		}
		grid.push_back(bin_i);
		step *= bin_i.size();
		
		//std::cout << " while done " << std::endl;
	}
	//std::cout << "all vars done" << std::endl;
	//invert grid ( correct order.. )
	//std::cout << grid << std::endl;
	flip(grid);
	//std::cout << grid << std::endl;
	//its done.
	*this = NDGrid(grid, file);
	
}

NDGrid::NDGrid(std::vector<std::vector<double>> bins, std::string file){
	//invoke standard constructor:
	*this = NDGrid(bins);
	//read values in:
	fill_from_file(file);
}

NDGrid::NDGrid(std::vector<double> bins1d){
	//invoke standard constructor:
	std::vector<std::vector<double>> bins {bins1d};
	*this = NDGrid(bins);
}
// part of constructor

void NDGrid::reSetGrid(){
	//assumes values have been set. Now it can set lengths, aux_lengths, and the (empty!) grid.

	dim = values.size(); // dimension of grid
	lengths = std::vector<unsigned int>(dim,0);
	aux_lengths = lengths;
	
	int size_of_grid = 1; // whole length of 1d projection of bins
	for(unsigned i = 0; i < dim; ++i){
		lengths[i] = values[i].size(); //save the lengths of individual bins
		size_of_grid *= (lengths[i]-1); //whole grid size, the #bins in a grid = length-1
		
	}
	grid = std::vector<int>(size_of_grid, 0);
	
	
	//create aux sizes to easier compute the entry-coordinates
	int aux_size = 1;
	for(unsigned i = 0; i < dim; ++i){
		// invert counter and i am to stupid to do it directly
		int j = dim -i -1;
		aux_lengths[j] = aux_size;
		aux_size *= (values[j].size()-1);
	}
}

//MEMBER FUNCTIONS

void NDGrid::newGrid(std::vector<std::vector<double>> newgrid){

	values = newgrid;
	
	reSetGrid();
	return;
}

void NDGrid::newGrid(std::vector<double> newgrid){
	std::vector<std::vector<double>> tmp;
	tmp.push_back(newgrid);
	newGrid(tmp);
	return;
}

unsigned int NDGrid::grid_size(){
	return grid.size();
}

std::vector<std::vector<double>> NDGrid::copy_grid(){
	return values;
}
int NDGrid::access_position(std::vector<unsigned int> positions){
	return grid[index(positions)];
}

void NDGrid::fill_from_file(std::string file){

	//setup stream
	std::ifstream ifs;
	ifs.open(file);
	checkIFile(ifs, file);
	//cheap solution: if the file is created properly, then line "i" is index "i" in my grid. Thus simply proceed as follows.
	std::vector<std::string> line_vec;
	
	for(unsigned int i = 0; i < this->grid.size(); ++i){
		
		getline(ifs, line_vec, delimiter);
		
		
		this->grid[i] = stoi(line_vec[line_vec.size() - 1]);	//last element of line is counter. All previous elements are binning info, but in this cheap solution i expect that the binning matches.
		
	}
	ifs.close();
}

void NDGrid::add_file(std::string f){
	//assume file f has SAME structure as already existing NDGrid. If not this function does not make sense, and will most likely return an error message.. if you are lucky.
	NDGrid tmp(f);
	*this += tmp;
	return;
}

std::vector<unsigned int> NDGrid::compute_indices(std::vector<double> vals){
	
	std::vector<unsigned int> positions(dim, 0);
	
	for(unsigned int i = 0; i < dim; ++i){
		const std::vector<double> bins = values[i];
		//print1dvec(bins);
		//const double val = vals[i];
		
		
		positions[i] = getBinIndex(values[i], vals[i]); // find the position of value i (= vals[i], input of this function) in the binning vector i (=values[i])
		
		
	}
	//print1dvec(positions);
	return positions;
}

bool NDGrid::pointCovered(std::vector<double> &vals){

	for(unsigned int i = 0; i < vals.size(); ++i){
		double val = vals[i];
		if(val < values[i].front() or val > values[i].back())
			return false;
	}
	return true;
}

void NDGrid::fill(std::vector<double> vals){
	if(pointCovered(vals))
		grid[index(vals)] += 1;
}

void NDGrid::fill(double val){
	std::vector<double> vals {val};
	fill(vals);
}

void NDGrid::raise_counter(std::vector<unsigned int> positions){
	grid[index(positions)] += 1;
}

void NDGrid::set_counter(std::vector<unsigned int> positions, const int counter){
	grid[index(positions)] = counter;
}

unsigned int NDGrid::index(std::vector<unsigned int> positions){
	unsigned int pos = 0;
	for(unsigned int i = 0; i < dim; ++i){
		pos += positions[i] * aux_lengths[i];
	}
	return pos;
}

unsigned int NDGrid::index(std::vector<double> vals){

	return index(compute_indices(vals));
}

void NDGrid::print_singleLine(std::ofstream &f){
	for(unsigned int i = 0; i < grid.size()-1; ++i){
		f << grid[i] << delimiter;
	}
	f << grid[grid.size()-1] << std::endl; // last element and end line
	return;
}

/*void NDGrid::read_from_file(std::string& adress){
	std::ifstream f;
	for(unsigned int i = 0; i < grid.size()-1; ++i){
		f << grid[i] << delimiter;
	}
	f << grid[grid.size()-1] << std::endl; // last element and end line
	return;
}*/ // I will think of implementing this once i think it is relevant. Right now i think this is highly unneccessary.

void NDGrid::print(std::ofstream &f){
	for(unsigned int i = 0; i < grid.size(); ++i){
		//invert index i:
		
		std::vector<unsigned int> index = multiindex(i);
		// now print binning
		for(unsigned int j = 0; j < dim; ++j){
			f << values[j][index[j]] << delimiter << values[j][index[j]+1] << delimiter;
		}
		//print value
		f << grid[i] << std::endl;
	}
	return;
}
/* REMOVED FROM AGENDA
void NDGrid::print_sidis_format(std::ofstream &f){ // this means the very last grid is pid. this is my convention
	for(unsigned int i = 0; i < grid.size(); i+=lengths[dim-1]){
		//invert index i:
		
		std::vector<unsigned int> index = multiindex(i);
		// now print binning
		for(unsigned int j = 0; j < dim-1; ++j){
			f << values[j][index[j]] << delimiter << values[j][index[j]+1] << delimiter;
		}
		//print values for all pids
		for(unsigned int l = 0; l < lengths[dim-1]; ++l){
			f << grid[i+l];
			if(l != lengths[dim-1]-1){
				f << delimiter;}
		}
		f << endl;
	}
	return;
}*/

std::vector<unsigned int> NDGrid::multiindex(unsigned int index){
	std::vector<unsigned int> multiindex (dim, 0);
	for(unsigned int j = 0; j < dim; ++j){
		multiindex[j] =  index / aux_lengths[j];
		
		index = index % aux_lengths[j];
	}
	return multiindex;
}
/*REMOVED FROM AGENDA
void NDGrid::print_sidis_format(std::string &adress){
	std::ofstream f;
	f.open(adress);
	print_sidis_format(f);
	f.close();
	return;
}*/

void NDGrid::print(std::string &adress){
	std::ofstream f;
	f.open(adress);
	print(f);
	f.close();
	return;
}
//the following again will be extremely complicated while also useless if one just does it as 4 d vector explicitly...
/*NDGrid_plus NDGrid::createM(NDGrid &dis, std::vector<bool> positions){ // take this grid and divide by DIS grid, output will be then multiplicity grid
	assert(positions.size() == dim);
	// the bool are : true if i have to divide by it and false if not, but if not divide by it i should normalize by the binsize there. This is multiplicity generation
	
	std::vector<double> mlt (grid.size(),0.);
	for(unsigned int i = 0; i < grid.size(); ++i){
		std::vector<unsigned int> index = multiindex(i);
		//create multiindex for dis
		std::vector<unsigned int> dis_ndx (dis.grid.size(), 0);
		mlt[i] = double(grid[i])/dis.grid()
	}
	
}*/

void NDGrid::print_normalized(std::string &normalized_file,const int nEvents){

	std::vector<double> ngrd (grid.size(),0.);
	
	for(unsigned int i = 0; i < grid.size(); ++i){
		double binsize = 1.;
		
		std::vector<unsigned int> indx = multiindex(i);
		
		for(unsigned int j = 0; j < dim; ++j){
			double d_var = values[j][indx[j]+1] - values[j][indx[j]];
			binsize *= d_var; 
		}
		
		ngrd[i] = double(grid[i]) / (nEvents * binsize);
	}
	
	std::ofstream f;
	f.open(normalized_file);
	
	for(unsigned int i = 0; i < grid.size(); ++i){
		//invert index i:
		std::vector<unsigned int> index = multiindex(i);
		// now print binning
		for(unsigned int j = 0; j < dim; ++j){
			f << values[j][index[j]] << delimiter << values[j][index[j]+1] << delimiter;
		}
		//print norm. value with also counter info
		//f << ngrd[i] << delimiter << grid[i] << endl;
		//print norm. value without counter info
		f << ngrd[i] << std::endl;
	}
	f.close();
	return;
}


/*void ::NDGrid::info(std::string &file){
	std::ofstream f;
	f.open(file); 
}*/
//OPERATORS

NDGrid NDGrid::operator=(const NDGrid& other){ //copy everything.
	if(this == & other) return *this;
	
	this->grid = other.grid;
	this->values = other.values;
	this->lengths = other.lengths;
	this->aux_lengths = other.aux_lengths;
	this->dim = other.dim;
	return *this;
}

NDGrid* NDGrid::operator+=(const NDGrid & other){
	//check if both have the same binning
	assert( this->dim == other.dim);
	for(unsigned int i = 0; i < this->dim; ++i){
		assert( this->values[i].size() == other.values[i].size());
		for(unsigned int j = 0; j < this->values[i].size(); ++j){
		
			assert( this->values[i][j] == other.values[i][j]);
		}
	}
	std::vector<int> a = other.grid;
	
	for(unsigned int i = 0; i < this->grid.size(); ++i){ // I USE THIS TO COMPILE WITH RIVET
		this->grid[i] += other.grid[i];
	}
	
	//this->grid = this->grid + a; //this version gave me trouble for compilation with RIVET
	
	return this;
}

NDGrid NDGrid::integrate(const unsigned int variable_order){
	//create new grid:
	
	//probably this should not be done if i integrate a one dimensional grid and would obtain just a number.
//	std::cout << "integrating position " << variable_order << " in " << dim << " dimensional grid." << std::endl;
	std::vector<std::vector<double>> grid = this->copy_grid(); // copy the grid from *this*. Isolate the target variable and then delete it from the grid:
//	std::cout << grid << std::endl;
	const std::vector<double> target_var = grid[variable_order];
	grid.erase(grid.begin() + variable_order);
	
	//now assign new values.
	NDGrid new_obj(grid);
	
	for(unsigned int i = 0; i < new_obj.grid_size(); ++i){
//		std::cout << "i = " << i << std::endl;
		std::vector<unsigned int> new_pos = new_obj.multiindex(i);
		int new_val = 0;
		std::vector<unsigned int> old_pos = new_pos;
		old_pos.insert(begin(old_pos) + variable_order, 0);
		
		for(unsigned int j = 0; j < target_var.size()-1; ++j){
//			std::cout << "j = " << j << std::endl;
			old_pos[variable_order] = j;
			
			new_val += this->access_position(old_pos);
			
			
						
			/*double measure = (target_var[j+1] - target_var[j]) ? normalize_yes_no : 1.; makes no sense in integer-valued counter class.
			//warning
			if(measure < 0){
				std::cout << "Warning: negative measure in integration. The vector:" << std::endl;
				print1dvec(target_var);
			}*/
			
			
		}
		
		new_obj.set_counter(new_pos, new_val);
		
	}
	
	return new_obj;
	
}

/* auf eis gelegt
NDGrid NDGrid::integrate(const std::vector<unsigned int> &variable_orders){

	//order: decreasing, integrate over last bin first..
	//if(!isOrdered(variable_orders,'>'){
	//	orderVector(variable_orders,'>');
	//}
	
	NDGrid temp1(), temp2();
	
	temp1 = *this; // copy
	
	for(unsigned int i = 0; i => variable_orders.size(); ++i){
		std::cout << i < std::endl;
		temp2 = temp1.integrate(variable_orders[i]);
		temp1 = temp2;
	}
	
}*/

namespace NDGridAF{
	//NON MEMBER FUNCTION THAT HOWERVER IS BOUND TO THIS CLASS
	/*void multiplicity_to_file(std::string file, NDGrid &sidis, NDGrid &dis, std::vector<bool> &vars, const char delimiter){
		multiplicity_to_file(*file, sidis, dis, vars, delimiter);
		return;
	}*/

	void print_normed(std::string &file, NDGrid &hist, std::vector<bool> &vars, double globalNormV, const char delimiter){ // essentially i copied the code of the function "multiplicity_to_file"...
		std::vector<std::vector<double>> grid = hist.copy_grid();
		std::vector<double> normvals (hist.grid_size(), 0.);
		for(unsigned int i = 0; i < hist.grid_size(); ++i){
			std::vector<unsigned int> indices = hist.multiindex(i);
			double val = double(hist.access_position(indices));
			double measure = 1;
			for(unsigned int j = 0; j < vars.size(); ++j){
				if(vars[j]){
				//means this variable has to be normalized;
					std::vector<double> v = grid[j];
					measure *= (v[indices[j]+1] - v[indices[j]]); // measure of this bin of the j-th variable.
					}
			}
			measure *= globalNormV;
			
			normvals[i] = val / measure;
		}
		std::ofstream f;
		f.open(file);
		for(unsigned int i = 0; i < hist.grid_size(); ++i){
			std::vector<unsigned int> indices = hist.multiindex(i);
			for(unsigned int j = 0; j < grid.size(); ++j){
				f << grid[j][indices[j]] << delimiter << grid[j][indices[j]+1] << delimiter; // lets trust this 
			}
			f << normvals[i];
			f << std::endl;
		}

		return;
	}

	void multiplicity_to_file(std::string &file, NDGrid &sidis, NDGrid &dis, std::vector<bool> &vars, const char delimiter){
		std::vector<std::vector<double>> grid = sidis.copy_grid();
		std::vector<double> mult (sidis.grid_size(), 0.); // = sidis grid with zeros. now we compute the values
		const unsigned int nPID = grid[grid.size()-1].size()-1;  // last grid dimension is for pid. the size of this vector -1 (because its also binned!) then is the amount of particle types i care about in this analysis. Certainly not the most intuitive formula but correct.
		
		for(unsigned int i = 0; i < sidis.grid_size(); ++i){
			std::vector<unsigned int> indices = sidis.multiindex(i);
			double val = double(sidis.access_position(indices));
			double measure = 1;
			for(unsigned int j = 0; j < vars.size(); ++j){
				if(vars[j]){
				//means this variable has to be normalized;
					std::vector<double> v = grid[j];
					measure *= (v[indices[j]+1] - v[indices[j]]); // measure of this bin of the j-th variable.
					}
			}
			std::vector<unsigned int> dis_indices = {indices[0], indices[1]}; // this can be adapted if dis should be somehow different organized, or whatever. However right now, this is what i need.
			int dis_number = dis.access_position(dis_indices);
			mult[i] = val / (measure * dis_number);
		}
		
		std::ofstream f;
		f.open(file);
		checkOFile(f,file);
		
		for(unsigned int i = 0; i < sidis.grid_size(); i+=nPID){
			// move forward in the stepsize of amount of pids, because i would like to store them all in one file but next to each other. this will make it easier to plot with my routines...
			std::vector<unsigned int> indices = sidis.multiindex(i);
			for(unsigned int j = 0; j < grid.size()-1; ++j){
			//CAUTION: gridsize is the size of the LOCAL vector. Thus the gridsize here equals the dimension of the sidis grid!
			//CAUTION: here is -1 because the last grid is to be understood as the pid of the particle. thus we dont need to output any grid here (even thought i use a "grid" for it to store it.
				f << grid[j][indices[j]] << delimiter << grid[j][indices[j]+1] << delimiter;
			}// now the binning is made clear
			for(unsigned int ipid = 0; ipid < nPID; ++ipid){
				f << mult[i+ipid];
				if(ipid != nPID-1) f << delimiter; // such that the delimiter is only printed if its not the last element of the line.
			}
			f << std::endl;
		}
		
		
		f.close();

		return;
	}

	//void NDG_to_oldform(string &infile, string& outfile, std::vector<int> &pidList){ use just amoung of pids, thats the only info that matters here.
	void NDG_to_oldform(std::string &infile, std::string& outfile, const unsigned int amountOfPidBins){
		
		const char delimiter = '	'; // tab
		
		std::ifstream ifs; 
		std::ofstream ofs;
		
		ifs.open(infile);
		checkIFile(ifs,infile);
		
		ofs.open(outfile);
		checkOFile(ofs,outfile);
		
		std::string line;
		
		while(getline(ifs, line)){ // loop over document, and already read the first line here.

			std::vector<std::string> line_vec, aux_vec;
			line_to_vec(line, line_vec, delimiter);
			
			//remove pid binning in vector. The very last entry of my vector is the counter - its is the important information. The two entries in front of the last entry however are the pid binning and i want to remove it and collect it in one line.
			
			line_vec.erase(line_vec.end()-3);//the pid binning slots are removed #1
			line_vec.erase(line_vec.end()-2);//#2
		
			// now add the other pid lines, but only the counter info
			for(unsigned int i = 0; i < amountOfPidBins - 1; ++i){// i already have the first entry, so size -1 elements remaining!
			
				getline(ifs, aux_vec, delimiter);
				
				line_vec.push_back(aux_vec[ aux_vec.size()-1 ] ); // add the last element which is the counter.
			
				
			}
		
			// now write the line, then move on to the next kinematic bin
			write_vec_to_stream(ofs, line_vec, delimiter); 
		}
		
		ofs.close();
		ifs.close();
		
		return;
	}

	void isolate_bin(const std::vector<double> &values, const std::vector<unsigned int> &positions, const std::string &sourcefile, const std::string &outfile){
		//safetychecks
		if(values.size() != positions.size()){
			//std::cout << "a problem occurred in \"isolate_bin\": different vector lengths.\nFirst vector(values):" << endl;
			error("a problem occurred in \"isolate_bin\": different vector lengths.\nFirst vector(values):");
			print1dvec(values);
			std::cout << "Second vector(positions):" << std::endl;
			print1dvec(positions);
			
			
			exit(1);
		}

		std::ifstream ifs; 
		std::ofstream ofs;
		
		ifs.open(sourcefile);
		checkIFile(ifs,sourcefile);
		
		ofs.open(outfile);
		checkOFile(ofs,outfile);
		
		std::string line;
		
		const char delimiter = '	';
		
		while(getline(ifs, line)){
			//std::cout << line << std::endl;
			bool printline = true;
			
			
			std::vector<std::string> line_vec;
			line_to_vec(line, line_vec, delimiter);
			std::vector<double> line_vec_d = v_stod(line_vec);
			
			//print1dvec(line_vec_d);
			
			//std::cout << line << std::endl;
			
			
			for(unsigned int i = 0; i < positions.size(); ++i){
			
				unsigned int j = positions[i];
				//std::cout << line_vec_d[j] << "	" << values[i] << std::endl;
				if(line_vec_d[j] != values[i]){
					//std::cout << "doesntfit" << std::endl;
					printline = false;
					//std::cout << line << std::endl;
					break;
				}
				
			}
			
			if(printline){
				//std::cout << line << std::endl;
				ofs << line << std::endl;
			}
		}
		
		return;
	}
	
	void isolate_bin_reduce(const std::vector<double> &values, const std::vector<unsigned int> &positions, const std::string &sourcefile, const std::string &outfile){
		//safetychecks
		if(values.size() != positions.size()){
			//std::cout << "a problem occurred in \"isolate_bin\": different vector lengths.\nFirst vector(values):" << std::endl;
			error("a problem occurred in \"isolate_bin\": different vector lengths.\nFirst vector(values):");
			print1dvec(values);
			std::cout << "Second vector(positions):" << std::endl;
			print1dvec(positions);
			
			
			exit(1);
		}
		std::vector<unsigned int> ordered_positions = order_vector(positions, '<'); // this returns a vector with increasing values. 
		
		
		std::ifstream ifs; 
		std::ofstream ofs;
		
		ifs.open(sourcefile);
		checkIFile(ifs,sourcefile);
		
		ofs.open(outfile);
		checkOFile(ofs,outfile);
		
		std::string line;
		
		const char delimiter = '	';
		
		while(getline(ifs, line)){
			
			bool printline = true;
			
			
			std::vector<std::string> line_vec;
			line_to_vec(line, line_vec, delimiter);
			std::vector<double> line_vec_d = v_stod(line_vec);
			
			//print1dvec(line_vec_d);
			
			//cout << line << endl;
			
			
			for(unsigned int i = 0; i < positions.size(); ++i){
			
				unsigned int j = positions[i];
				//cout << line_vec_d[j] << "	" << values[i] << endl;
				if(line_vec_d[j] != values[i]){
					//cout << "doesntfit" << endl;
					printline = false;
					break;
				}
				
			}
			
			if(printline){
				for(unsigned int i = ordered_positions.size()-1; i >= 0; --i){ //start from vector end.
					line_vec.erase(line_vec.begin() + ordered_positions[i] + 1); 
					line_vec.erase(line_vec.begin() + ordered_positions[i]);
				}
				std::string outline = vec_to_line(line_vec, delimiter);
				ofs << outline << std::endl;
			}
		}
		
		return;
	}

}
