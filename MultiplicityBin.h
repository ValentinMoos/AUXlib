#ifndef MULTIPLICITYBIN
#define MULTIPLICITYBIN


/*
Changed on 08.06.2021
features: 
	Designed for logarithmic SIDIS analysis/plots/histograms. CAUTION if use for other purpose.
	Use constructors to initialize values. To avoid user( that is, my own) errors
	
*/

using namespace std;
#include"Aux.h"


class MultiplicityBin{
	private:
		//variables:
		const char delimiter = '	';
		int PID;
		
		//grid
		std::vector<double> Q2values; //Q2 values, from which i obtain the bining in Q2. Analog for the other 3 variables
		std::vector<double> xvalues;
		std::vector<double> zvalues;
		std::vector<double> ptvalues;
		
		//4d grid
		std::vector<std::vector<std::vector<std::vector<double>>>> counter_grid;         // this is the bin counter. It is double instead of integer so i could copy it to obtain the other bins directly
		std::vector<std::vector<std::vector<std::vector<double>>>> multiplicity_grid;    // multiplicity grid in all 4 variables. It is obtained by taking the SIDIS event counter (already normalized by bin size) and divide any bin (Q2,x,z,pt) by the amount of (also by bin sized normalized -) DIS events in the same (Q2, x) region.
		std::vector<std::vector<std::vector<std::vector<double>>>> intermediate_grid;    //  this is for intermediate "steps", but right now only for the step from couter -> counter / bin size.
		std::vector<std::vector<double>> two_dim_grid;  // placeholder for outputs in 2 variables (so its always a projection (integration) from 4 d to 2 d
		std::vector<double> one_dim_grid;               // placeholder for outputs in 1 variable  (so its always a projection (integration) from 4 d to 1 d
		
		void reset_grid(vector<double>& vec, unsigned int length); 	                                 // reset the placeholder, such that it then has the desired length and all entries are 0
		void reset_grid(vector<vector<double>>& vec, unsigned int length1, unsigned int length2);    // reset the placeholder, such that it then has the desired length(2d) and all entries are 0
		
		void Integrate2(string A, string B, vector<vector<double>> &twodgrid);            // integrate 2 variables out of 4 -> 2d grid remaining
		void Integrate3(string A, string B, string C, vector<double> &onedgrid);          // integrate 3 variables out of 4 -> 1d grid remaining
		void Integrate3_counter(string A, string B, string C, vector<double> &onedgrid);  // integrate 3 variables in the counter. I think this was for testing only.
		
	public:
		MultiplicityBin(const string &Q2binAdress, const string &xbinAdress, const string &zbinAdress, const string &ptbinAdress, const int particleId); //constructor, reads grids and creates an empty grid of the needed size in 4 dimension
		void addValue(const double Q2, const double x, const double z, const double pt, const int particleId); // adds the value to the right bin if the particle ID is correct for this Hist
		void normalize_bins(); //divides by bin size
		void create_multiplicities(std::vector<std::vector<double>> &DIS_grid);    // divide the bins by the respective DIS counter(normalized)
		void create_multiplicities2(std::vector<std::vector<double>> &DIS_grid);    // divide the bins by the respective DIS counter(unormalized, true counter)
		
		void print_toStream_twodimGrid(ofstream &f, string var1, string var2);        //prints 2d Hist in the variables given to a steam f
		void print_toStream_onedimGrid(ofstream &f, string var);                      //prints 1d Hist in the variable  given to a steam f
		
		void printToFile(std::string &, std::vector<string> &);                       //prints to the file (first argument) the Histogram projection in the variables of the 2nd argument [var1,var2] (it is checked if one or two variables are passed on (or more, but then it returns and error message)
		
		void print_Hist(std::string variables, std::string adress);
		void printToFile_two_variable_mult(std::string file, std::vector<double> &vec1, std::vector<double> &vec2);
		void printToStream_two_variable_mult(std::ofstream &f, std::vector<double> &vec1, std::vector<double> &vec2);
		
		void printToFile_one_variable_mult(std::string file, std::vector<double> &vec1);
		void printToStream_one_variable_mult(std::ofstream &f, std::vector<double> &vec1);
		
		void printAllToDirectory(const std::string & s); //prints "all" (most) Histograms to a file. Since the projection in the experiment depends on the binning, most files produced wont make too much sense.
		
		void print_toTerminal_twodimGrid(string var1, string var2);
		void print_toTerminal_onedimGrid(string var);
		void print_toTerminal_onedimGrid_counter(string var);
		
		void printFullGridTerminal(const char Mode); // print the full grid, either c -> counter, m -> multiplicity, i -> intermediate
		void saveGrid(const char Mode, std::string &adress); // save the grid with options  either c -> counter, m -> multiplicity, i -> intermediate (best will be to save the counter so its transparent whats going on.)
		
		void read_grid(const char mode, std::string &adress); // reads grid from adress
		
		int acceptedEvents; // counter how many events are accepted in the Histogram. Since i had trouble putting this to private for a still undiscovered reason this is public.
		//int access_acceptedEvents();
		void printBin(std::string); // to show that the bins are correct.
		void print_all_bins_Terminal(); // print all bins to show that the bins are correct.
		//void print(std::string, std::string &location); // 
		void printToFile_two_variable_mult_onefixed(std::string directoryAdress, std::vector<double> &vec_fixed, unsigned int fixPos, std::vector<double> &vec_running, bool correctOrder); // prints to file the 2 d grid which one of the grid axis fixed at the position. The boolean is needed such that i know whether to acess the grid[i,j] or [j,i].
		void printToFile_two_variable_mult_allCombis(std::string adress, std::vector<double> &vec_fixedVal, std::vector<double> & vec_runningVal, bool correctOrder); // evokes other function to print to file the 2 d grid. This will cause a file for each of the values of the first vector (vec_fixedVal). the boolean is to access the grid in correct order.
};

#endif
