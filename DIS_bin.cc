#include"DIS_bin.h"

DIS_bin::DIS_bin(const string &Q2binAdress, const string &xbinAdress)
{
	// read the bin files and with it generate the vectors containing the information.
	string line = "emptyLine";
	
	read_1dFile_to_vec(Q2binAdress, Q2values);
	read_1dFile_to_vec(xbinAdress, xvalues);
	
	
	
	//Now all grids are known
	//I define the order in grids as Q2-x
	//Create Q2-x grid
	std::vector<double> testvector;
	for(unsigned int l=0; l<xvalues.size(); ++l)
	{
		testvector.push_back(0.);
	}
	std::vector<std::vector<double>> testvector2;
	for(unsigned int l=0; l<Q2values.size(); ++l)
	{
		counter_grid.push_back(testvector);
	}
	// conveniently create other grids as copies of the first one
	normalized_grid = counter_grid;
}

DIS_bin::DIS_bin()
{
	// only placeholder.
}

DIS_bin& DIS_bin::operator=(const DIS_bin& other)
{	
	//self assignment
	if( this == &other)
		return *this;
	//serious assignment of other DIS
	normalized = other.normalized;
	Q2values   = other.Q2values;
	xvalues    = other.xvalues;
	counter_grid = other.counter_grid;
	normalized_grid = other.normalized_grid;
	return *this;
}

void DIS_bin::normalize_bins()
{
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			//cout << i << "-" << j << "-" << k << "-" << l << endl;
			double dQ2 = Q2values[i+1]-Q2values[i];
			double dx  =  xvalues[j+1]- xvalues[j];
			if(dQ2 <= 0. || dx <= 0.)
			{
				cout << i << "-" << j << "	Please check the binning." << endl; 
			}
			normalized_grid[i][j] = counter_grid[i][j]/(dQ2*dx);
		}
	}
	normalized = true;
	return;
	
	
}

void DIS_bin::reset_grid(vector<vector<double>>& vec)
{
	for(decltype(vec.size()) i = 0; i < vec.size(); ++i)
	{
		for(decltype(vec[i].size()) j = 0; j < vec.size(); ++j)
		{
			vec[i][j]=0.;
		}
	}
}

void DIS_bin::addValue(const double Q2, const double x)
{
	int Q2pos=-1;
	for(unsigned int i=0; i<Q2values.size()-1; ++i)
	{
		if(Q2>Q2values[i] && Q2<=Q2values[i+1])
		{
			Q2pos=i;
		}
	}
	
	int xpos=-1;
	for(unsigned int i=0; i<xvalues.size()-1; ++i)
	{
		if(x>xvalues[i] && x<=xvalues[i+1])
		{
			xpos=i;
		}
	}
	
	if(Q2pos==-1 || xpos == -1)
	{
		//cannot fit the value in the grid
		return;
	}
	else
	{
		counter_grid[Q2pos][xpos]+=1;	
	}
	return;
	
}

void DIS_bin::print_toStream_Grid(ofstream &f)
{
	for(unsigned int i=0; i < Q2values.size()-1; ++i)
	{
		for(unsigned int j=0;j < xvalues.size()-1; ++j)
		{
			f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << xvalues[j] << delimiter << xvalues[j+1] << delimiter << counter_grid[i][j] << endl;
		}
	}
}

std::vector<std::vector<double>>& DIS_bin::accessGrid(const char c)
{
	switch(c){
	case 'c':
		return counter_grid;
		break;
	case 'n':
		if(normalized)
			return normalized_grid;
		else
		{
			cout << "The bins are not normalized yet." << endl;
			exit(202);
		}
	default :
		cout << "Unknown option in DIS_bin::accessGrid with argument " << c << endl;
		exit(202);
	}
	
}

void DIS_bin::printTo_Terminal(const char mode)
{
	std::vector<std::vector<double>> grid;
	switch(mode)
	{
		case 'c':
		grid = counter_grid;
		break;
		
		case 'n':
		grid = normalized_grid;
		break;
		
		case 'm':
		grid = normalized_grid;
		break;
		
		case 'i':
		grid = counter_grid;
		break;
		
		default:
		cout << "Unknown option in DIS_bin::printTo_Terminal with argument " << mode << endl;
		exit(2); 
	}
	for(unsigned int i=0; i < Q2values.size()-1; ++i)
	{
		for(unsigned int j=0;j < xvalues.size()-1; ++j)
		{
			cout << Q2values[i] << delimiter << Q2values[i+1] << delimiter << xvalues[j] << delimiter << xvalues[j+1] << delimiter << grid[i][j] << endl;
		}
	}
}

void DIS_bin::read_grid(const char mode, std::string &adress)
{
	auto grid = counter_grid;
		
	
	ifstream f;
	f.open(adress);
	checkIFile(f, adress);
	
	std::string line;
	
	while(getline(f,line))
	{
		//cout << line << endl;
		std::vector<std::string> line_vec;
		splitLineToArray(line, line_vec, delimiter);
		
		unsigned int i = getIndex(Q2values, stod(line_vec[0])); // Q2 positions are 0 and 1
		unsigned int j = getIndex( xvalues, stod(line_vec[2])); //  x positions are 2 and 3
		
		grid[i][j] = stod(line_vec[4]);                         // data position is 4
		
	}
	f.close();
	switch(mode)
	{
		case 'c':
		counter_grid = grid;
		break;
		
		case 'n':
		normalized_grid = grid;
		break;
		
		default:
		cout << "Sidis_cluster::read_grid	unknown option: mode \"" << mode << "\"" << endl;
		exit(2);
	}
}


void DIS_bin::saveGrid(const char mode, const std::string &adress)
{
	std::vector<std::vector<double>> grid;
	switch(mode)
	{
		case 'c':
		grid = counter_grid;
		break;
		
		case 'n':
		grid = normalized_grid;
		break;
		
		default:
		cout << "DIS_bin::saveGrid unknown option " << mode << endl;
		exit(2);
	}
	ofstream f;
	f.open(adress);
	for(unsigned int i=0; i < Q2values.size()-1; ++i)
	{
		for(unsigned int j=0;j < xvalues.size()-1; ++j)
		{
			f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << xvalues[j] << delimiter << xvalues[j+1] << delimiter << grid[i][j] << endl;
		}
	}
	f.close();
	
}

