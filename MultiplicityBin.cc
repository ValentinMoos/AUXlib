#include"MultiplicityBin.h"

MultiplicityBin::MultiplicityBin(const string &Q2binAdress, const string &xbinAdress, const string &zbinAdress, const string &ptbinAdress, const int particleID)
{
	// read the bin files and with it generate the vectors containing the information.
	string line = "emptyLine";
	
	
	//Q2
	
	read_1dFile_to_vec(Q2binAdress, Q2values);
	read_1dFile_to_vec(xbinAdress, xvalues);
	read_1dFile_to_vec(zbinAdress, zvalues);
	read_1dFile_to_vec(ptbinAdress, ptvalues);
	
	//Now all grids are known
	//I define the order in grids as Q2-x-z-pt
	//Create Q2-x-z-pt grid
	std::vector<double> testvector;
	for(decltype(ptvalues.size()) l=0; l<ptvalues.size(); ++l)
	{
		testvector.push_back(0.);
	}
	std::vector<std::vector<double>> testvector2;
	for(decltype(zvalues.size()) l=0; l<zvalues.size(); ++l)
	{
		testvector2.push_back(testvector);
	}
	std::vector<std::vector<std::vector<double>>> testvector3;
	for(decltype(xvalues.size()) l=0; l<xvalues.size(); ++l)
	{
		testvector3.push_back(testvector2);
	}
	
	for(decltype(Q2values.size()) l=0; l<Q2values.size(); ++l)
	{
		counter_grid.push_back(testvector3);
	}
	multiplicity_grid = counter_grid; // conveniently create other grids as copies of the first one
	intermediate_grid = counter_grid;
	//Define the necessary variables
	acceptedEvents =0;
	PID = particleID;
	
	if(!(check_sorted(Q2values) && check_sorted(xvalues) && check_sorted(zvalues) && check_sorted(ptvalues)))
	{
		cout << "Not all binning values are sorted. Check this." << endl;
		exit(1);
	}
	
}

/*
void MultiplicityBin::checkIFile(ifstream & ifstr, const string &adress)
{
	if(!ifstr.is_open())
	{
		cout << "Could not open the file with adress " << adress << " aborting process." << endl;
		exit(1);
	}
	return;
}
*/
/*
int MultiplicityBin::access_acceptedEvents()
{
	return acceptedEvents;
}
*/

void MultiplicityBin::normalize_bins()
{
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			for(decltype(zvalues.size()) k=0; k<zvalues.size()-1; ++k)
			{
				for(decltype(ptvalues.size()) l=0; l<ptvalues.size()-1; ++l)
				{
					//cout << i << "-" << j << "-" << k << "-" << l << endl;
					double dQ2 = Q2values[i+1]-Q2values[i];
					double dx  =  xvalues[j+1]- xvalues[j];
					double dz  =  zvalues[k+1]- zvalues[k];
					double dpt = ptvalues[l+1]-ptvalues[l];
					if(dQ2 <= 0. || dx <= 0. || dz <= 0. || dpt <= 0.)
					{
						cout << i << "-" << j << "-" << k << "-" << l << "	Please check the binning." << endl << "dQ2=" << dQ2 << "	dx=" << dx << "	dz=" << dz << "	dpt=" << dpt << endl; 
					}
					intermediate_grid[i][j][k][l] = counter_grid[i][j][k][l]/(dQ2*dx*dz*dpt);
				}
			}
		}
	}
}

void MultiplicityBin::printFullGridTerminal(const char Mode)
{
	std::vector<std::vector<std::vector<std::vector<double>>>> grid;
	switch(Mode)
	{
		case 'c':
		grid = counter_grid;
		break;
		
		case 'm':
		grid = multiplicity_grid;
		break;
		
		case 'i':
		grid = intermediate_grid;
		break;
		
		default:
		cout << "MultiplicityBin::printFullGridTerminal unknown option with argument " << Mode << endl;
		exit(2);
	}
	//const std::string space = " | ";
	const std::string space = "	"; //tab
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			for(decltype(zvalues.size()) k=0; k<zvalues.size()-1; ++k)
			{
				for(decltype(ptvalues.size()) l=0; l<ptvalues.size()-1; ++l)
				{
					cout << Q2values[i] << space << Q2values[i+1] << space << xvalues[j] << space << xvalues[j+1] << space << zvalues[k] << space << zvalues[k+1] << space << ptvalues[l] << space << ptvalues[l+1] << space << grid[i][j][k][l] << endl;
					//cout << i << space << j << space << k << space << l << space << grid[i][j][k][l] << endl;
				}
			}
		}
	}
}

void MultiplicityBin::read_grid(const char mode, std::string &adress)
{

	auto grid = counter_grid; // to get the structure, the values are being overwritten anyway.
	
	ifstream f;
	f.open(adress);
	checkIFile(f, adress);
	
	std::string line;
	
	while(getline(f,line))
	{
		std::vector<std::string> line_vec;
		splitLineToArray(line, line_vec, delimiter);
		
		unsigned int i = getIndex(Q2values, stod(line_vec[0])); // Q2 positions are 0 and 1
		unsigned int j = getIndex( xvalues, stod(line_vec[2])); //  x positions are 2 and 3
		unsigned int k = getIndex( zvalues, stod(line_vec[4])); //  z positions are 4 and 5
		unsigned int l = getIndex(ptvalues, stod(line_vec[6])); // pt positions are 6 and 7
		
		grid[i][j][k][l] = stod(line_vec[8]);                         // data position is 8
	}
	f.close();
	switch(mode)
	{
		case 'c':
		counter_grid = grid;
		break;
		
		case 'm':
		multiplicity_grid = grid;
		break;
		
		case 'i':
		intermediate_grid = grid;
		break;
		
		default:
		cout << "Sidis_cluster::read_grid	unknown option: mode \"" << mode << "\"" << endl;
		exit(2);
	}
}

void MultiplicityBin::saveGrid(const char Mode, std::string &adress)
{
	std::vector<std::vector<std::vector<std::vector<double>>>> grid;
	switch(Mode)
	{
		case 'c':
		grid = counter_grid;
		break;
		
		case 'm':
		grid = multiplicity_grid;
		break;
		
		case 'i':
		grid = intermediate_grid;
		break;
		
		default:
		cout << "MultiplicityBin::printFullGridTerminal unknown option with argument " << Mode << endl;
		exit(2);
	}

	ofstream f; 
	f.open(adress);	
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			for(decltype(zvalues.size()) k=0; k<zvalues.size()-1; ++k)
			{
				for(decltype(ptvalues.size()) l=0; l<ptvalues.size()-1; ++l)
				{
					f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << xvalues[j] << delimiter << xvalues[j+1] << delimiter << zvalues[k] << delimiter << zvalues[k+1] << delimiter << ptvalues[l] << delimiter << ptvalues[l+1] << delimiter << grid[i][j][k][l] << endl;
					//cout << i << space << j << space << k << space << l << space << grid[i][j][k][l] << endl;
				}
			}
		}
	}
	f.close();
}

void MultiplicityBin::addValue(const double Q2, const double x, const double z, const double pt, const int particleId)
{
	//if its the wrong particle the whole process can be skipped.
	if(particleId != PID)
	return;
	
	//get coordinates:
	int Q2pos=-1;
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		if(Q2>Q2values[i] && Q2<=Q2values[i+1])
		{
			Q2pos=i;
		}
	}
	
	int xpos=-1;
	for(decltype(xvalues.size()) i=0; i<xvalues.size()-1; ++i)
	{
		if(x>xvalues[i] && x<=xvalues[i+1])
		{
			xpos=i;
		}
	}
	
	int zpos=-1;
	for(decltype(zvalues.size()) i=0; i<zvalues.size()-1; ++i)
	{
		if(z>zvalues[i] && z<=zvalues[i+1])
		{
			zpos=i;
		}
	}
	
	int ptpos=-1;
	for(decltype(ptvalues.size()) i=0; i<ptvalues.size()-1; ++i)
	{
		if(pt>ptvalues[i] && pt<=ptvalues[i+1])
		{
			ptpos=i;
		}
	}
	if(Q2pos==-1 || xpos == -1 || zpos == -1 || ptpos== -1)
	{
		//cannot fit the value in the grid
		//cout << "coult not accept event with values: Q2=" << Q2 << "	x=" << x << "	z=" << z << "	pt=" << pt << "	Positions: " << Q2pos << xpos << zpos << ptpos << endl; 
		return;
		
	}
	else
	{
		counter_grid[Q2pos][xpos][zpos][ptpos]+=1;
		acceptedEvents++;
	}
	return;
}

void MultiplicityBin::create_multiplicities(std::vector<std::vector<double>> &DIS_grid)
{
	//normalized DIS grid is expected
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			for(decltype(zvalues.size()) k=0; k<zvalues.size()-1; ++k)
			{
				for(decltype(ptvalues.size()) l=0; l<ptvalues.size()-1; ++l)
				{
					//cout << i << "-" << j << "-" << k << "-" << l << endl;
					multiplicity_grid[i][j][k][l] = (intermediate_grid[i][j][k][l])/(DIS_grid[i][j]);
					if(DIS_grid[i][j]==0.)
					{
						cout << i << "-" << j << "-" << k << "-" << l << "	" << DIS_grid[i][j] << " entry in DIS seems to be zero, thus normalization fails." << endl; 
					}
				}
			}
		}
	}
}

void MultiplicityBin::create_multiplicities2(std::vector<std::vector<double>> &DIS_grid)
{
	//here the DIS grid is not normalized, it is just the counter
	for(decltype(Q2values.size()) i=0; i<Q2values.size()-1; ++i)
	{
		for(decltype(xvalues.size()) j=0; j<xvalues.size()-1; ++j)
		{
			if(DIS_grid[i][j]==0.)
			{
				cout << i << "-" << j << "	" << DIS_grid[i][j] << " entry in DIS seems to be zero, thus normalization fails." << endl; 
			}
			for(decltype(zvalues.size()) k=0; k<zvalues.size()-1; ++k)
			{
				double dz = zvalues[k+1] - zvalues[k];
				for(decltype(ptvalues.size()) l=0; l<ptvalues.size()-1; ++l)
				{
					//cout << i << "-" << j << "-" << k << "-" << l << endl;
					double dpT = ptvalues[l+1]-ptvalues[l];
					
					multiplicity_grid[i][j][k][l] = (counter_grid[i][j][k][l])/(DIS_grid[i][j]*dz*dpT);
				}
			}
		}
	}
}


void MultiplicityBin::printToFile_two_variable_mult_allCombis(std::string adress, std::vector<double> &vec_fixedVal, std::vector<double> & vec_runningVal, bool correctOrder)
{
	for(unsigned int i = 0; i < vec_fixedVal.size()-1; ++i)
	{
		printToFile_two_variable_mult_onefixed(adress, vec_fixedVal, i, vec_runningVal, correctOrder);
	}
}

void MultiplicityBin::printToFile_two_variable_mult_onefixed(std::string directoryAdress, std::vector<double> &vec_fixed, unsigned int fixPos, std::vector<double> &vec_running, bool correctOrder)
{
	std::string adress = directoryAdress + to_string(vec_fixed[fixPos]) + "-" + to_string(vec_fixed[fixPos+1]);
	ofstream f;
	f.open(adress);
	for(unsigned int i = 0; i < vec_running.size()-1; ++i)
	{
		if(correctOrder)
		f << vec_fixed[fixPos] << delimiter << vec_fixed[fixPos+1] << delimiter << vec_running[i] << delimiter << vec_running[i+1] << delimiter << two_dim_grid[fixPos][i] << endl;
		else
		f << vec_fixed[fixPos] << delimiter << vec_fixed[fixPos+1] << delimiter << vec_running[i] << delimiter << vec_running[i+1] << delimiter << two_dim_grid[i][fixPos] << endl;
	}
	
	f.close();
}

void MultiplicityBin::print_Hist(std::string variables, std::string adress)
{
	if(variables == "Q2-z")
	{
		cout << "print Q2-z" << endl;
		reset_grid(two_dim_grid, Q2values.size()-1, zvalues.size()-1);
		//integrate over x and pt:
		for(unsigned int i = 0; i < Q2values.size()-1;++i)
		{
			for(unsigned int j = 0; j < xvalues.size()-1; ++j)
			{
				double dx = xvalues[j+1] - xvalues[j];
				for(unsigned int k = 0; k < zvalues.size()-1; ++k)
				{
					for(unsigned int l = 0; l < ptvalues.size()-1; ++l)
					{
						double dpt = ptvalues[l+1] - ptvalues[l];
						two_dim_grid[i][k] += multiplicity_grid[i][j][k][l]*dx*dpt;	
					}
				}
			}
		}
		//now print to file
		printToFile_two_variable_mult(adress, Q2values, zvalues);
		//print all "sub possibilities too"
		std::string adress_for_fixed = adress + "Q2=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, Q2values, zvalues, true);
		adress_for_fixed = adress + "z=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, zvalues, Q2values, false);
		
	}
	if(variables == "x-z")
	{
		cout << "print x-z" << endl;
		reset_grid(two_dim_grid, xvalues.size()-1, zvalues.size()-1);
		//integrate over Q2 and pt:
		for(unsigned int i = 0; i < Q2values.size()-1;++i)
		{
			double dQ2 = Q2values[i+1] - Q2values[i];
			for(unsigned int j = 0; j < xvalues.size()-1; ++j)
			{
				for(unsigned int k = 0; k < zvalues.size()-1; ++k)
				{
					for(unsigned int l = 0; l < ptvalues.size()-1; ++l)
					{
						double dpt = ptvalues[l+1] - ptvalues[l];
						two_dim_grid[i][k] += multiplicity_grid[i][j][k][l]*dQ2*dpt;	
					}
				}
			}
		}		
		//now print to file
		printToFile_two_variable_mult(adress, xvalues, zvalues);
		//print all "sub possibilities too"
		std::string adress_for_fixed = adress + "x=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, xvalues, zvalues, true);
		adress_for_fixed = adress + "z=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, zvalues, xvalues, false);
	}
	if(variables == "z-pt")
	{
		cout << "print z-pt" << endl;
		reset_grid(two_dim_grid, zvalues.size()-1, ptvalues.size()-1);
		//integrate over Q2 and x:
		for(unsigned int i = 0; i < Q2values.size()-1;++i)
		{
			double dQ2 = Q2values[i+1] - Q2values[i];
			for(unsigned int j = 0; j < xvalues.size()-1; ++j)
			{
				double dx = xvalues[j+1] - xvalues[j];
				for(unsigned int k = 0; k < zvalues.size()-1; ++k)
				{
					for(unsigned int l = 0; l < ptvalues.size()-1; ++l)
					{
						two_dim_grid[i][k] += multiplicity_grid[i][j][k][l]*dQ2*dx;	
					}
				}
			}
		}		
		//now print to file
		printToFile_two_variable_mult(adress, zvalues, ptvalues);
		//print all "sub possibilities too"
		std::string adress_for_fixed = adress + "z=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, zvalues, ptvalues, true);
		adress_for_fixed = adress + "pt=";
		printToFile_two_variable_mult_allCombis(adress_for_fixed, ptvalues, zvalues, false);
	}
	if(variables == "z")
	{
		cout << "print z" << endl;
		reset_grid(one_dim_grid, zvalues.size()-1);
		//integrate over Q2, x and pt:
		for(unsigned int i = 0; i < Q2values.size()-1;++i)
		{
			double dQ2 = Q2values[i+1] - Q2values[i];
			for(unsigned int j = 0; j < xvalues.size()-1; ++j)
			{
				double dx = xvalues[j+1] - xvalues[j];
				for(unsigned int k = 0; k < zvalues.size()-1; ++k)
				{
					for(unsigned int l = 0; l < ptvalues.size()-1; ++l)
					{
						double dpt = ptvalues[l+1] - ptvalues[l];
						one_dim_grid[k] += multiplicity_grid[i][j][k][l]*dQ2*dx*dpt;	
					}
				}
			}
		}		
		//now print to file
		printToFile_one_variable_mult(adress, zvalues);
	}	
}
/*
void MultiplicityBin::Integrate_two(const std::string var1, const std::string var2)
{
	std::vector<double> vecA, vecB, vecC, vecD;
	std::string local_error_msg = "error in MultiplicityBin::Integrate_two with arguments: " + var1 + " and " + var2;
	//chose the correct vectors
	
	if(var1  == "Q2")
	{
		vecC = Q2values;
		if(var2 == "x")
		{
			vecD = xvalues;
			vecA = zvalues;
			vecB = ptvalues;
		}
		else
		{
			vecA = xvalues;
			if(var2 == "z")
			{
				vecD = zvalues;
				vecB = ptvalues;
			}
			else if(var2 == "pt")
			{
				vecD = ptvalues;
				vecB = zvalues;
			}
			else {cout << local_error_msg << endl; return 55;}
		}
	}
	else
	{
		vecA = Q2values;
		if(var1 == "x")
		{
			vecC = xvalues;
			if(var2 == "z")
			{
				vecD = zvalues;
				vecB = ptvalues;
			}
			else if(var2 == "pt")
			{
				vecD = ptvalues;
				vecB = zvalues;
			}
			else {cout << local_error_msg << endl; return 55;}
		}
		else if(var1 == "z" && var2 == "pt")
		{
			vecB = xvalues;
			vecC = zvalues;
			vecD = ptvalues;
		}
		else {cout << local_error_msg << endl; return 55;}
	}
	
	//now the variables are set. integrate the 4 d grid into a 2 d grid
	reset_grid(two_dim_grid, vecA.size()-1, vecB.size()-1);
	
	for(auto A = 0; A < vecA.size(); ++A)
		for(auto B = 0; B < vecB.size(); ++B)
			for(auto C = 0; C < vecC.size(); ++C)
			{
				double dC = vecC[C+1] - vecC[C];
				for(auto D = 0; D < vecD.size(); ++D)
				{
					double dD = vecD[D+1] - vecD[D];
					two_dim_grid[A][B] += multiplicity_grid[
				}
			}
	
	
	if(var1 == "Q2" && var2 == "x")
	{
		vecA = zvalues;
		vecB = ptvalues;
		vecC = Q2values;
		vecD = xvalues;
	}
	
	if(var1 == "Q2" && var2 == "z")
	{
		vecA = xvalues;
		vecB = ptvalues;
		vecC = Q2values;
		vecD = zvalues;
	}
	
	if(var1 == "Q2" && var2 == "x")
	{
		vecA = zvalues;
		vecB = ptvalues;
		vecC = Q2values;
		vecD = xvalues;
	}
	
	switch(var1)
	{
		case 'Q':
		vecA = xvalues;
		vecB = zvalues;
		vecC = ptvalues;
		vecD = Q2values;
		break;
		
		case 'x':
		vecA = Q2values;
		vecB = zvalues;
		vecC = ptvalues;
		vecD = xvalues;
		break;
		
		case 'z':
		vecA = xvalues;
		vecB = zvalues;
		vecC = ptvalues;
		vecD = Q2values;
		break;
	}
	
	
	
}*/

void MultiplicityBin::printToFile_two_variable_mult(std::string file, std::vector<double> &vec1, std::vector<double> &vec2)
{
	std::ofstream f;
	f.open(file);
	printToStream_two_variable_mult(f, vec1, vec2);
	f.close();
	return;
}

void MultiplicityBin::printToFile_one_variable_mult(std::string file, std::vector<double> &vec1)
{
	std::ofstream f;
	f.open(file);
	printToStream_one_variable_mult(f, vec1);
	f.close();
	return;
}

void MultiplicityBin::printToStream_two_variable_mult(std::ofstream &f, std::vector<double> &vec1, std::vector<double> &vec2)
{
	for(unsigned int m = 0; m < vec1.size()-1; ++m)
	{
		for(unsigned int n = 0; n < vec2.size()-1; ++n)
		{
			f << vec1[m] << delimiter << vec1[m+1] << delimiter << vec2[n] << delimiter << vec2[n+1] << delimiter << two_dim_grid[m][n] << endl;
		}
	}
	return;
}

void MultiplicityBin::printToStream_one_variable_mult(std::ofstream &f, std::vector<double> &vec1)
{
	for(unsigned int m = 0; m < vec1.size()-1; ++m)
	{
		f << vec1[m] << delimiter << vec1[m+1] << delimiter << one_dim_grid[m] << endl;
	}
	return;
}

void MultiplicityBin::printAllToDirectory(const std::string & directory)
{
	std::vector<std::string> variables;
	string file;
	//2d binned results
	
	//Q2-z
	variables.push_back("Q2");
	variables.push_back("z");
	file = directory + variables[0] +"-" + variables[1] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();	
	//Q2-pt
	variables.push_back("Q2");
	variables.push_back("pt");
	file = directory + variables[0] +"-" + variables[1] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//x-z
	variables.push_back("x");
	variables.push_back("z");
	file = directory + variables[0] +"-" + variables[1] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//x-pt
	variables.push_back("x");
	variables.push_back("pt");
	file = directory + variables[0] +"-" + variables[1] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//z-pt
	variables.push_back("z");
	variables.push_back("pt");
	file = directory + variables[0] +"-" + variables[1] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	
	//1d binned results
	//Q2
	variables.push_back("Q2");
	file = directory + variables[0] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//x
	variables.push_back("x");
	file = directory + variables[0] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//z
	variables.push_back("z");
	file = directory + variables[0] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	//pt
	variables.push_back("pt");
	file = directory + variables[0] +"_PID="+to_string(PID)+"_mult.dat";
	printToFile(file,variables);
	variables.clear();
	
	return;	
	
}

void MultiplicityBin::printToFile(std::string & file, std::vector<string> & var)
{
	cout << "Printing to File " << file << endl;
	bool printed = false;
	ofstream f;
	f.open(file);
	if(!f.is_open())
	{
		cout << "Could not open the file with adress " << file << " aborting process." << endl;
		exit(1);
	}
	if(var.size() == 1)
	{
		//cout << "Printing variable " << var[0] << endl;
		print_toStream_onedimGrid(f, var[0]);
		printed = true;
	}
	if(var.size() == 2)
	{
		//cout << "Printing variables " << var[0] << "	" << var[1] << endl;
		print_toStream_twodimGrid(f, var[0], var[1]);
		printed = true;
	}
	if(!printed)
	{
		cout << "Error in MultiplicityBin::printToFile:	does not know what to do with file " << file << " and vector of size " << var.size() << endl;
	}
}

void MultiplicityBin::print_all_bins_Terminal()
{
	cout << "##################################################################################" << endl;
	cout << "this is MultiplicityBin with PID " << PID << "\nThe binning is:" << endl;
	cout << "Q2 - binning" << endl;
	printBin("Q2");
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "xB - binning" << endl;
	printBin("x");
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "z  - binning" << endl;
	printBin("z");
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "pt - binning" << endl;
	printBin("pt");
	cout << "##################################################################################" << endl;
}

void MultiplicityBin::printBin(std::string variable)
{
	std::vector<double> vec;
	bool accepted = false;
	if(variable == "Q2")
	{
		vec = Q2values;
		accepted =true;
	}
	if(variable == "x")
	{
		vec = xvalues;
		accepted =true;
	}
	if(variable == "z")
	{
		vec = zvalues;
		accepted =true;
	}
	if(variable == "pt")
	{
		vec = ptvalues;
		accepted =true;
	}
	if(!accepted)
	{
		cout << "what to do with \"" << variable << "\" in printBin()? \n Abort. " << endl;
		return;
	}
	cout << "Print bin of variable: " << variable << endl;
	for(unsigned int i = 0; i < vec.size(); ++i)
	{
		cout << i << delimiter << vec[i] << endl;
	}
	
	cout << "therefore have the bins: " << endl;
	
	for(unsigned int i = 0; i < vec.size()-1; ++i)
	{
		cout << i << delimiter << vec[i] << "-" << vec[i+1] << endl;
	}
	
}


void MultiplicityBin::reset_grid(vector<double>& vec, unsigned int length)
{
	vec.clear();
	for(decltype(vec.size()) i = 0; i < length; ++i)
	{
		vec.push_back(0.);
	}
}

void MultiplicityBin::reset_grid(vector<vector<double>>& vec, unsigned int length1, unsigned int length2)
{
	vec.clear();
	vector<double> bufferVector;
	for(decltype(vec.size()) i = 0; i < length2; ++i)
	{
		bufferVector.push_back(0.);
	}
	
	for(decltype(vec.size()) j = 0; j < length1; ++j)
	{
		vec.push_back(bufferVector);
	}
}

void MultiplicityBin::Integrate2(string A, string B, vector<vector<double>> &twodgrid)
{
	//cout << "MultiplicityBin::Integrate2 with arguments " << A << " and " << B << endl; 
	if(A == "Q2" && B == "pt")
	{
		reset_grid(twodgrid, xvalues.size()-1, zvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dpt = ptvalues[l+1]-ptvalues[l];
						twodgrid[j][k] += multiplicity_grid[i][j][k][l]*dQ2*dpt;
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "z")
	{
		reset_grid(twodgrid, xvalues.size()-1, ptvalues.size()-1);
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dz  =  zvalues[k+1]- zvalues[k];
						twodgrid[j][l] += multiplicity_grid[i][j][k][l]*dQ2*dz;
					}
				}
			}
		}
	}
	if(A == "x" && B == "pt")
	{
		reset_grid(twodgrid, Q2values.size()-1, zvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dx  =  xvalues[j+1]- xvalues[j];
						double dpt = ptvalues[l+1]-ptvalues[l];
						twodgrid[i][k] += multiplicity_grid[i][j][k][l]*dx*dpt;
					}
				}
			}
		}
	}
	if(A == "x" && B == "z")
	{
		reset_grid(twodgrid, Q2values.size()-1, ptvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dx  =  xvalues[j+1]- xvalues[j];
						double dz  =  zvalues[k+1]- zvalues[k];
						twodgrid[i][l] += multiplicity_grid[i][j][k][l]*dx*dz;
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "x")
	{
		reset_grid(twodgrid, zvalues.size()-1, ptvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dx  =  xvalues[j+1]- xvalues[j];
						twodgrid[k][l] += multiplicity_grid[i][j][k][l]*dQ2*dx;
					}
				}
			}
		}
	}
	//cout << A << " " << B << " have been integrated." << endl;
}



void MultiplicityBin::Integrate3(string A, string B, string C, vector<double> &onedgrid)
{
	
	if(A == "Q2" && B == "x" && C == "pt")
	{
		reset_grid(onedgrid, zvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dx  =  xvalues[j+1]- xvalues[j];
						double dpt = ptvalues[l+1]-ptvalues[l];
						onedgrid[k] += multiplicity_grid[i][j][k][l]*dQ2*dx*dpt;
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "x" && C == "z")
	{
		reset_grid(onedgrid, ptvalues.size()-1);
	
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dx  =  xvalues[j+1]- xvalues[j];
						double dz  =  zvalues[k+1]- zvalues[k];
						onedgrid[l] += multiplicity_grid[i][j][k][l]*dQ2*dx*dz;
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "z" && C == "pt")
	{
		reset_grid(onedgrid, xvalues.size()-1);
		
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dQ2 = Q2values[i+1]-Q2values[i];
						double dz  =  zvalues[k+1]- zvalues[k];
						double dpt = ptvalues[l+1]-ptvalues[l];
						onedgrid[j] += multiplicity_grid[i][j][k][l]*dQ2*dz*dpt;
					}
				}
			}
		}
	}
	if(A == "x" && B == "z" && C == "pt")
	{
		reset_grid(onedgrid, Q2values.size()-1);
	
		for(unsigned int i=0; i<Q2values.size()-1; ++i)
		{
			for(unsigned int j=0; j<xvalues.size()-1; ++j)
			{
				for(unsigned int k=0; k<zvalues.size()-1; ++k)
				{
					for(unsigned int l=0; l<ptvalues.size()-1; ++l)
					{
						double dx  =  xvalues[j+1]- xvalues[j];
						double dz  =  zvalues[k+1]- zvalues[k];
						double dpt = ptvalues[l+1]-ptvalues[l];
						onedgrid[i] += multiplicity_grid[i][j][k][l]*dx*dz*dpt;
					}
				}
			}
		}
	}
	//cout << A << " " << B << " " << C << " have been integrated." << endl;
}


void MultiplicityBin::print_toStream_onedimGrid(ofstream &f, string var)
{
	if(var == "Q2")
	{
		Integrate3("x", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "x")
	{
		Integrate3("Q2", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			f << xvalues[i] << delimiter << xvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "z")
	{
		#ifdef DEBUG
		cout << "printing z now" << endl;
		#endif
		Integrate3("Q2", "x", "pt", one_dim_grid);
		for(unsigned int i=0; i < zvalues.size()-1; ++i)
		{
			f << zvalues[i] << delimiter << zvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "pt")
	{
		Integrate3("Q2", "x", "z", one_dim_grid);
		for(unsigned int i=0; i < ptvalues.size()-1; ++i)
		{
			f << ptvalues[i] << delimiter << ptvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	
}

void MultiplicityBin::print_toStream_twodimGrid(ofstream &f, string var1, string var2)
{
	if((var1 == "Q2" && var2 == "z") || (var1 == "z" && var2 == "Q2"))
	{
		Integrate2("x","pt", two_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			for(unsigned int j=0;j < zvalues.size()-1; ++j)
			{
				f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << zvalues[j] << delimiter << zvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "Q2" && var2 == "pt") || (var1 == "pt" && var2 == "Q2"))
	{
		Integrate2("x","z", two_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				f << Q2values[i] << delimiter << Q2values[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "x" && var2 == "z") || (var1 == "z" && var2 == "x"))
	{
		Integrate2("Q2","pt", two_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < zvalues.size()-1; ++j)
			{
				f << xvalues[i] << delimiter << xvalues[i+1] << delimiter << zvalues[j] << delimiter << zvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "x" && var2 == "pt") || (var1 == "pt" && var2 == "x"))
	{
		Integrate2("Q2","z", two_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				f << xvalues[i] << delimiter << xvalues[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "z" && var2 == "pt") || (var1 == "pt" && var2 == "z"))
	{
		Integrate2("Q2","x", two_dim_grid);
		for(unsigned int i=0; i < zvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				f << zvalues[i] << delimiter << zvalues[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
}

void MultiplicityBin::print_toTerminal_onedimGrid(string var)
{
	if(var == "Q2")
	{
		Integrate3("x", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			cout << Q2values[i] << delimiter << Q2values[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "x")
	{
		Integrate3("Q2", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			cout << xvalues[i] << delimiter << xvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "z")
	{
		#ifdef DEBUG
		cout << "printing z now" << endl;
		#endif
		Integrate3("Q2", "x", "pt", one_dim_grid);
		for(unsigned int i=0; i < zvalues.size()-1; ++i)
		{
			cout << zvalues[i] << delimiter << zvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "pt")
	{
		Integrate3("Q2", "x", "z", one_dim_grid);
		for(unsigned int i=0; i < ptvalues.size()-1; ++i)
		{
			cout << ptvalues[i] << delimiter << ptvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	
}

void MultiplicityBin::print_toTerminal_twodimGrid(string var1, string var2)
{
	if((var1 == "Q2" && var2 == "z") || (var1 == "z" && var2 == "Q2"))
	{
		Integrate2("x","pt", two_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			for(unsigned int j=0;j < zvalues.size()-1; ++j)
			{
				cout << Q2values[i] << delimiter << Q2values[i+1] << delimiter << zvalues[j] << delimiter << zvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "Q2" && var2 == "pt") || (var1 == "pt" && var2 == "Q2"))
	{
		Integrate2("x","z", two_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				cout << Q2values[i] << delimiter << Q2values[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "x" && var2 == "z") || (var1 == "z" && var2 == "x"))
	{
		Integrate2("Q2","pt", two_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < zvalues.size()-1; ++j)
			{
				cout << xvalues[i] << delimiter << xvalues[i+1] << delimiter << zvalues[j] << delimiter << zvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "x" && var2 == "pt") || (var1 == "pt" && var2 == "x"))
	{
		Integrate2("Q2","z", two_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				cout << xvalues[i] << delimiter << xvalues[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
	if((var1 == "z" && var2 == "pt") || (var1 == "pt" && var2 == "z"))
	{
		Integrate2("Q2","x", two_dim_grid);
		for(unsigned int i=0; i < zvalues.size()-1; ++i)
		{
			for(unsigned int j=0;j < ptvalues.size()-1; ++j)
			{
				cout << zvalues[i] << delimiter << zvalues[i+1] << delimiter << ptvalues[j] << delimiter << ptvalues[j+1] << delimiter << two_dim_grid[i][j] << endl;
			}
		}
	}
}



void MultiplicityBin::Integrate3_counter(string A, string B, string C, vector<double> &onedgrid)
{
	
	if(A == "Q2" && B == "x" && C == "pt")
	{
		reset_grid(onedgrid, zvalues.size());
		
		for(unsigned int i=0; i<Q2values.size(); ++i)
		{
			for(unsigned int j=0; j<xvalues.size(); ++j)
			{
				for(unsigned int k=0; k<zvalues.size(); ++k)
				{
					for(unsigned int l=0; l<ptvalues.size(); ++l)
					{
						onedgrid[k] += counter_grid[i][j][k][l];
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "x" && C == "z")
	{
		reset_grid(onedgrid, ptvalues.size());
	
		for(unsigned int i=0; i<Q2values.size(); ++i)
		{
			for(unsigned int j=0; j<xvalues.size(); ++j)
			{
				for(unsigned int k=0; k<zvalues.size(); ++k)
				{
					for(unsigned int l=0; l<ptvalues.size(); ++l)
					{
						onedgrid[l] += counter_grid[i][j][k][l];
					}
				}
			}
		}
	}
	if(A == "Q2" && B == "z" && C == "pt")
	{
		reset_grid(onedgrid, xvalues.size());
		
		for(unsigned int i=0; i<Q2values.size(); ++i)
		{
			for(unsigned int j=0; j<xvalues.size(); ++j)
			{
				for(unsigned int k=0; k<zvalues.size(); ++k)
				{
					for(unsigned int l=0; l<ptvalues.size(); ++l)
					{
						onedgrid[j] += counter_grid[i][j][k][l];
					}
				}
			}
		}
	}
	if(A == "x" && B == "z" && C == "pt")
	{
		reset_grid(onedgrid, Q2values.size());
	
		for(unsigned int i=0; i<Q2values.size(); ++i)
		{
			for(unsigned int j=0; j<xvalues.size(); ++j)
			{
				for(unsigned int k=0; k<zvalues.size(); ++k)
				{
					for(unsigned int l=0; l<ptvalues.size(); ++l)
					{
						onedgrid[i] += counter_grid[i][j][k][l];
					}
				}
			}
		}
	}
	//cout << A << " " << B << " " << C << " have been integrated." << endl;
}

void MultiplicityBin::print_toTerminal_onedimGrid_counter(string var)
{
	if(var == "Q2")
	{
		Integrate3_counter("x", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < Q2values.size()-1; ++i)
		{
			cout << Q2values[i] << delimiter << Q2values[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "x")
	{
		Integrate3_counter("Q2", "z", "pt", one_dim_grid);
		for(unsigned int i=0; i < xvalues.size()-1; ++i)
		{
			cout << xvalues[i] << delimiter << xvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "z")
	{
		#ifdef DEBUG
		cout << "printing z now" << endl;
		#endif
		Integrate3_counter("Q2", "x", "pt", one_dim_grid);
		for(unsigned int i=0; i < zvalues.size()-1; ++i)
		{
			cout << zvalues[i] << delimiter << zvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	if(var == "pt")
	{
		Integrate3_counter("Q2", "x", "z", one_dim_grid);
		for(unsigned int i=0; i < ptvalues.size()-1; ++i)
		{
			cout << ptvalues[i] << delimiter << ptvalues[i+1] << delimiter << one_dim_grid[i] << endl;
		}
	}
	
}
/*
void Multiplicity::print(std::string s, std::string location)
{
	//first find out what to print
		
	if(s == "Q2-z")
	{
		std::vector<std::string> variables = ["Q2","z"];
		printToFile(location, variables);	
	}
	
	return;
}
*/


