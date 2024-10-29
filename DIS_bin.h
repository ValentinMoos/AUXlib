#ifndef DIS_BIN
#define DIS_BIN


#include"Aux.h"

class DIS_bin{
	private:
		//variables:
		const char delimiter = '	';
		bool normalized = false;
		//int PID;
		
		std::vector<double> Q2values;
		std::vector<double> xvalues;
		
		std::vector<std::vector<double>> counter_grid;
		std::vector<std::vector<double>> normalized_grid;
		std::vector<std::vector<double>> two_dim_grid;
		
		void reset_grid(vector<vector<double>>& vec);
		
		
	public:
		
		
		DIS_bin(); // default constructor to initialize it, it is only here to allow the initialization of this object while not using it!
		
		DIS_bin(const string &Q2binAdress, const string &xbinAdress); //constructor, reads grids
		void addValue(const double Q2, const double x);
		void normalize_bins();
		void print_toStream_Grid(ofstream &f);
		
		std::vector<std::vector<double>>& accessGrid(const char); //return the grid, but via this is can probably also be converted, so probably it would be good to set it as constant at some point...
		
		void printTo_Terminal(const char mode);
		void saveGrid(const char mode, const std::string &adress);
		void read_grid(const char mode, std::string &adress); //read DIS grid from adress
		
		DIS_bin& operator=(const DIS_bin& other);  //makes it possible to assign an other DIS_bin to this one.
};

#endif
