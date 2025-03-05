#ifndef NDGRID
#define NDGRID

#include "Aux.h"
#include <vector>

/*idea of this class:
	get a grid, which is always ONE vector of multiple vectors, defining a more dimensional grid (in general).
	the class then will project this grid always to a 1d counter grid (therefore int - values).
	the grid is accessed via functions that understand the binning.
	
	Example: have two dim binning v1 = [0,2,4,6], v2 = [0.3,0.7,1.0]
	this class should then generate a grid with 3 x 2 = 6 entries.
	then take two exemplary values i have to add to my histogram: 
	a = [3, 0.8] -> 2nd bin in v1 and 2nd bin in v2. Thus put it into bin number #v1binnumber-1 x #sizeofsubbins + #v2binnumber-1 x #sizeofsubbins(=1 for v2 because its the "last" grid)  = 1 x 2  +  1 x 1 = 3 
	#binnumber-1 is because in c++ 0 is the first number so first bin is actually bin 0
	b = [1, 0.4] -> 1st bin in v1 and 1st bin in v2. Thus put it into bin number #v1binnumber-1 x #sizeofsubbins + #v2binnumber-1 x #sizeofsubbins = 0 + 0 = 0 = very first bin
	c = [5, 0.9] -> 3rd bin in v1 and 2nd bin in v2. Thus put it into bin number #v1binnumber-1 x #sizeofsubbins + #v2binnumber-1 x #sizeofsubbins = 4 + 1 = 5 = very last bin


*/
class NDGrid{

	private:
		//member
		const char delimiter = '	'; //delimiter, TAB
		std::vector<int> grid;  //the grid
		
		std::vector<std::vector<double>> values; // = the input bins
		
		std::vector<unsigned int> lengths;                               // lengths of the single 1-d vectors
		std::vector<unsigned int> aux_lengths;                           // lengths for binning, biggest to lowest (=1) stepsize
	
		unsigned int dim;                                                 // dimension = # of variables, size of "values"
		//function	
		
		std::vector<unsigned int> compute_indices(std::vector<double> vals);     // 
		
		//void add(unsigned int i); // add +1 to the position i of the vector
		void reSetGrid();

	public:
		NDGrid(std::vector<std::vector<double>> bins);              // create with grid and initialize the members
		NDGrid(std::vector<double> bins1d);					 // create 1d grid by calling general case. Simlifies life, does not change the mechanism.
		NDGrid();                                                    // default to create with no member-info
		NDGrid(std::vector<std::vector<double>> bins, std::string file); // read from file. This should first create the grid with the standard procedure ( IS NOT default! (i already forgot what i meant with "not default". this calls the NDGrid(std::vector<std::vector<double>> ) constructor, as one can check.) ) and then read the values from the file.
		NDGrid(std::string file); // read from a file, which has the file structure of an output of a NDGrid. So it can understand the binning and then read in the values.
		
		NDGrid operator=(const NDGrid& other);                       // copy object	
		NDGrid* operator+=(const NDGrid & other);					 // add grid of other to this ( sum entries )
		NDGrid integrate(const unsigned int variable_order);   //integrate over the variable with order "variable_order"
//		NDGrid integrate(const std::vector<unsigned int> &variable_orders);   //call above multiple times..
		
		void fill(std::vector<double> vals);                         // raise counter in bin depending on values
		void fill(double vals);                                      // raise counter in bin depending on value -- for convenience in 1d case. SL, DNCTM!!
		void raise_counter(std::vector<unsigned int> positions);     // raise counter in bin depending on indices   
		void fill_from_file(std::string);                            // read from a file all the counter values. The file must be in standard format of this grid - the standard output format, which is NOT oldform!
		void add_file(std::string);
		void set_counter(std::vector<unsigned int> positions, const int counter); // set a counter to a certain value (counter) instead of raising it by 1.
		void newGrid(std::vector<std::vector<double>> newgrid);	     //(re)define the binning. Clears all entries.
		void newGrid(std::vector<double> newgrid);                   // 1d case
		
		unsigned int index(std::vector<double> vals);                // get index of values
		unsigned int index(std::vector<unsigned int> positions);     // get index of indices
		std::vector<unsigned int> multiindex(unsigned int index);    // inversion of "get index of indices"
		
		int access_position(std::vector<unsigned int> positions);    // return value of certain position
		
		
		unsigned int grid_size();                                    // return grid.size()
		bool pointCovered(std::vector<double> &vals);                 // tells you if a point is covered by the range in which the histogram can fill values.
		
		void print_singleLine(std::ofstream &f); // print to stream (in one line)
		// DO NOT USE THIS FUNCTION - IT PRODUCES A CORRUMPT FILE void print_sidis_format(std::string &adress); // print into ofstream created linked to adress
		// DO NOT USE THIS FUNCTION - IT PRODUCES A CORRUMPT FILE void print_sidis_format(std::ofstream &f); // print to stream all kinematic binning and then all pids with the respective counter. special function for my sidis purpose.
		void print(std::ofstream &f); // print to stream fully listed all bins.
		void print(std::string &adress); //print into ofstream created linked to adress
		void print_normalized(std::string &normalized_file, const int nEvents);
		std::vector<std::vector<double>> copy_grid(); //returns "values"
		//void info(std::string &file); // print some information about this NDGrid into a file. Not urget or necessary atm.
};

namespace NDGridAF{ // auxiliary functions that I do not want to be member function because that is not necessary here.
/*
void multiplicity_to_file(std::string file, NDGrid &sidis, NDGrid &dis, std::vector<bool> &vars, const char delimiter); */// read sidis and dis grid, and compute and then print the multiplicities.

void print_normed(std::string &file, NDGrid &hist, std::vector<bool> &vars, double globalNormV, const char delimiter); // prints to file the hist, the vars that are true will be normalized to bin size and also every  bin will be divided by the last argument(globalNormV).

void print_normedInverse(std::string &file, NDGrid &hist, std::vector<bool> &vars, double globalNormV, const char delimiter); // prints to file the hist, the vars that are true will be normalized to bin size INVERSE (multiplied by size) and also every  bin will be divided by the last argument(globalNormV).

void multiplicity_to_file(std::string &file, NDGrid &sidis, NDGrid &dis, std::vector<bool> &vars, const char delimiter); // same function but takes string as reference. the above function calls this one.

//void NDG_to_oldform(string &infile, string& outfile, std::vector<int> &pidList); // convert to old format.

void NDG_to_oldform(std::string &infile, std::string& outfile, const unsigned int amount_of_pids_to_bin); // convert to old format.

void isolate_bin(const std::vector<double> &values, const std::vector<unsigned int> &positions, const std::string &sourcefile, const std::string &outfile); //project onto variables fixed by first two arguments. for example with values = {7,4} and positions {0,3} the function will look throught thefile "sourcefile" and save all lines in "outfile" that have a 7 on the first(index 0) and a 4 on the fourth(index 3) position (assumes a data file with delimiter TAB atm.).
}
#endif
