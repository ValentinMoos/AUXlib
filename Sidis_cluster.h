#ifndef SIDIS_CLUSTER
#define SIDIS_CLUSTER

#include"Aux.h"


#include"DIS_bin.h"
#include"MultiplicityBin.h"

class Sidis_cluster{

	private:
		std::vector<int> PID_list;
		
		
		std::vector<MultiplicityBin> cluster;
		
		DIS_bin DIS;
	
		std::string find_bin(const std::string &directory, std::string variable); // find the bins in the given directory automatically to simplify the construction of an object.
	
	public:
		Sidis_cluster(std::vector<int> & pid_list, const std::string & Q2grid, const std::string & xgrid, const std::string & zgrid, const std::string & ptgrid);//give all binning files seperately as references
		Sidis_cluster(std::vector<int> & pid_list, const std::string directory);  // give directory where I expect all names (Q2, x, z, pt) (or + "_assumed")

		void addDIS(const double Q2, const double x); // add event to DIS
		void addSIDIS(const double Q2, const double x, const double z, const double pt, const int pid); // add event to SIDIS grids
		void normalize_to_multiplicity(); // normalize all sidis grids by dividing them by the dis grids;
		void normalize_to_multiplicity2(); // normalize by counter grid of DIS and only grid sizes in z and pt
		void normalize_grids_by_size(); // normalize all grids by dividing them by the bin size. This affects all sidis grids and the dis grid;
		
		void print_Hists(std::string, std::string location); // print all hists (all pid) with that variable combination to the directory "location"
		void print_Hists2(std::string, std::string location); // print all hists (all pid) with that variable combination to the directory "location"
		
		void printFullGridTerminal(const char Mode); // print the full grid, either c -> counter, m -> multiplicity, i -> intermediate
		void saveGrids(const char Mode, const std::string & adress); // save grids with option either c -> counter, m -> multiplicity, i -> intermediate (best will be to save just the counter first.)  The adress will be a directory, and in there will be stored: the DIS counter grid, and all the counter grids for SIDIS for each particle ID.
		
		void print_bins_Terminal(); //prints all bins to the Terminal for checking, so it loops over all "member" MultiplicityBins and calls them to print their info to terminal
		
		void read_grid(const char mode, const std::string &directory, const std::string &DIS_grid, const std::string &SIDIS_names); // read in values from stored files instead of counting by running programm
		
};

#endif
