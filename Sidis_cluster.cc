#include"Sidis_cluster.h"

Sidis_cluster::Sidis_cluster(std::vector<int> &pid_list, const std::string & Q2grid, const std::string & xgrid, const std::string & zgrid, const std::string & ptgrid)
{
	PID_list = pid_list;
	//create DIS grid
	DIS_bin DIS(Q2grid, xgrid);
	//DIS = DIS_bin(Q2grid, xgrid);
	
	// create SIDIS grids
	
	for(unsigned int i = 0; i < pid_list.size(); ++i)
	{
		cluster.push_back(MultiplicityBin(Q2grid, xgrid, zgrid, ptgrid, pid_list[i]));
	}
	
}


Sidis_cluster::Sidis_cluster(std::vector<int> & pid_list, const std::string directory)
{
	//find the adresses
	std::string Q2grid = find_bin(directory,"Q2");
	std::string  xgrid = find_bin(directory,"xB");
	std::string  zgrid = find_bin(directory,"z");
	std::string ptgrid = find_bin(directory,"pT");;
	
	PID_list = pid_list;
	//create DIS grid
	DIS_bin DIS_constructed(Q2grid, xgrid);
	DIS = DIS_constructed;
	//DIS = DIS_bin(Q2grid, xgrid);
	
	// create SIDIS grids
	
	for(unsigned int i = 0; i < pid_list.size(); ++i)
	{
		cluster.push_back(MultiplicityBin(Q2grid, xgrid, zgrid, ptgrid, pid_list[i]));
	}
	
	//Sidis_cluster(pid_list, Q2grid, xgrid, zgrid, ptgrid);
}

void Sidis_cluster::read_grid(const char mode, const std::string &directory, const std::string &DIS_grid, const std::string &SIDIS_names)
{
		
	// read DIS grid
	std::string adress = directory + DIS_grid;
	DIS.read_grid(mode, adress);
	
	// read SIDIS grid
	for(unsigned int i = 0; i < PID_list.size(); ++i)
	{
		adress = directory + to_string(PID_list[i]) + SIDIS_names;
		cluster[i].read_grid(mode, adress);
	}
}

/*
Sidis_cluster::Sidis_cluster(std::vector<int> & pid_list, const std::string directory)
	: Sidis_cluster{
		pid_list,
		find_bin(directory,"Q2"),
		find_bin(directory,"xB"),
		find_bin(directory,"z"),
		find_bin(directory,"pT"),
	}
{}
*/

std::string Sidis_cluster::find_bin(const std::string &directory, std::string variable)
{
	std::string adress = directory + variable;
	if(check_File_exists(adress))
		return adress;
	else if(check_File_exists(adress + "_assumed"))
		return adress + "_assumed";
	else{ cout << "The direcotry " << directory << " does not seem to contain the needed binning for the variable " << variable << endl; exit(46);}
}

void Sidis_cluster::addDIS(const double Q2, const double x)
{
	DIS.addValue(Q2,x);
	return;
}


void Sidis_cluster::addSIDIS(const double Q2, const double x, const double z, const double pt, const int pid)
{
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		if(pid == PID_list[i])
		{
			cluster[i].addValue(Q2, x, z, pt, pid);
		}
	}
}
void Sidis_cluster::saveGrids(const char Mode, const std::string & adress)
{
	char DIS_Mode;
	switch(Mode)
	{
		case 'c':
		DIS_Mode = Mode;
		break;
	
		case 'n':
		DIS_Mode = Mode;
		break;
				
		case 'i':
		DIS_Mode = 'c';
		break;
	
		default:
		cout << "Sidis_cluster::saveGrids unknown option: " << Mode << endl;
		exit(2);
	}
	std::string DIS_adress = adress + "DISgrid";
	DIS.saveGrid(DIS_Mode, DIS_adress);
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		std::string cluster_adress = adress + to_string(PID_list[i]) + "_SIDISgrid";
		cluster[i].saveGrid(Mode, cluster_adress);
	}
	
}


void Sidis_cluster::normalize_grids_by_size()
{
	//DIS
	DIS.normalize_bins();	
	
	//Multiplicities
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		cluster[i].normalize_bins();
	}
}

void Sidis_cluster::normalize_to_multiplicity2()
{
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		cluster[i].create_multiplicities2(DIS.accessGrid('c'));
	}
}

void Sidis_cluster::normalize_to_multiplicity()
{
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		cluster[i].create_multiplicities(DIS.accessGrid('n'));
	}
}



void Sidis_cluster::print_Hists(std::string variables, std::string location)
{
	// create vector from variables string;
	const char delimiter = '-';
	std::vector<std::string> variables_vec;
	splitLineToArray(variables, variables_vec, delimiter);
	
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		std::string filename = location + to_string(PID_list[i]) + variables + "_mult.dat";
		cluster[i].printToFile(filename, variables_vec);
	}
}


void Sidis_cluster::print_Hists2(std::string variables, std::string location)
{
	// create vector from variables string;
	
	for(unsigned int i = 0; i < cluster.size(); ++i)
	{
		std::string filename = location + to_string(PID_list[i]) + variables + "_v2_mult.dat";
		cluster[i].print_Hist(variables, filename);
	}
}

void Sidis_cluster::printFullGridTerminal(const char Mode)
{
	cout << "##########################################################################" << endl;
	cout << "DIS" << endl;
	DIS.printTo_Terminal(Mode);
	cout << "##########################################################################" << endl; 
	for(unsigned int i = 0; i< PID_list.size(); ++i)
	{
		cout << "##########################################################################" << endl << PID_list[i] << " has the grid: (option " << Mode << "):" << endl;
		cluster[i].printFullGridTerminal(Mode);
		cout << "##########################################################################" << endl;
	}
}


void Sidis_cluster::print_bins_Terminal()
{
	//DIS bins:
	DIS.printTo_Terminal('c');	
	
	//Multiplicity bins:
	for(unsigned int i=0; i<PID_list.size(); ++i)
	{
		cluster[i].print_all_bins_Terminal();
	}
}
