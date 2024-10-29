#include"myHist.h"

myHist::myHist(const std::string type, const std::string &gridAdress, const double &Q2, const double &Q2min, const double &Q2max)
{
	#ifdef DEBUGREADGRID
	cout << "myHist: constructor" << endl;
	#endif
	Q2bin    = Q2;
	Q2binMin = Q2min;
	Q2binMax = Q2max;
	//first determine which typ of Histogram I want: linear, logarithmic etc.
	if(type == "lin") // create linear Hist
	{
		//invent some function
		cout << "this option does not exist yet." << endl;
		return 0; // as long as there is no function
	}
	if(type == "ReducedCrossSection:logx") // create Log Hist for reduced Cross section, some Q2 bin (given by arguments) and the x Grid from the Grid (gridAdress)
	{
		readBinsPLUSasVal(gridAdress, Q2min, Q2max, Q2); //reads grids and creates it
		#ifdef DEBUGREADGRID
		cout << "myHist: constructor: log Hist has been initialized" << endl;
		cout << "The member variables are: " << endl;
		cout << "sizeOfArrays	" << sizeOfArrays << endl;
		cout << "Q2bin	" << Q2bin << endl;
		cout << "Q2binMax	" << Q2binMax << endl;
		cout << "Q2binMin	" << Q2binMin << endl;
		/*cout << "	" <<  << endl;
		cout << "	" <<  << endl;
		cout << "	" <<  << endl;
		cout << "	" <<  << endl;
		cout << "	" <<  << endl;*/
		#ifdef MANUALSTOP
		getchar(); 
		#endif
		#endif
		
	}
	checkInitialisation();
}

void myHist::checkInitialisation()
{
	//check if Histogram is ok: Bins exist,
	if(sizeOfArrays<=0 || Q2binMin >= Q2binMax || Q2bin < Q2binMin || Q2bin > Q2binMax  ) // anything that might have gone wrong..
	{
		cout << "myHist::checkInitialisation	Something went wrong. Here the member-values: " << endl;
		printMembers();
		
		getchar();
		exit(1);
	}
	return;
}


template<class TYPE>
void myHist::printArray(const int &size, TYPE *array) // print Array of arbitrary type.
{
	for(int i = 0; i < size; ++i)
	{
		cout << "entry: " << i << " : " << array[i] << endl; 
	}
	return;
}


void myHist::printMembers()
{
	cout << "sizeOfArrays: " << sizeOfArrays << endl;
	cout << "Q2bin: " << Q2bin << endl;
	cout << "Q2binMin: " << Q2binMin << endl;
	cout << "Q2binMax: " << Q2binMax << endl;
	cout << "mySigmaNorm: " << mySigmaNorm << "	(shoudl be 0 at start)" << endl;
	cout << "failedEntries: " << failedEntries << "(should be 0 at start)" << endl;
	cout << "acceptedEntries: " << acceptedEntries << "(should be 0 at start)" << endl;
	cout << "infoString " << infoString << "(ending of info data name)" << endl;
	cout << "delimiter: Z" << delimiter << "Z (delimiter(symbol between Zs (might be a space or TAB)) in Histogram Output and also for read in Grids)" << endl;
	cout << "unwantedsymbols: " << unwantedsymbols << "(list of symbols for which the grid file is checked before reading in. Should avoid user confusion of space with TAB, or \".\" with \",\")" << endl;
	
	cout << "range: " << endl;
	printArray(sizeOfArrays, range);
	
	cout << "rangeTOP: " << endl;
	printArray(sizeOfArrays, rangeTOP);
	
	cout << "unnormalizedvalue: " << endl;
	printArray(sizeOfArrays, range);
	
	cout << "normalizedvalue: " << endl;
	printArray(sizeOfArrays, range);
	
	cout << "intermediatevale: " << endl;
	printArray(sizeOfArrays, range);
	
	cout << "binAssignedValue: " << endl;
	printArray(sizeOfArrays, range);
	
	return;
}

void myHist::printCounterOfBin(int i)
{
	cout << unnormalizedvalue[i] << endl;
	return;
}
void myHist::printValueOfBin(int i)
{
	cout << normalizedvalue[i] << endl;
	return;
}

void myHist::checkSourceFileforsymbol(std::string symbols, std::string adress)
{
	std::ifstream f;	
	
	std::string line;
	int linecounter = 1;
	f.open(adress);
	while(getline(f,line))
	{
		for(long unsigned int i=0; i<line.length(); ++i)
		{
			for(long unsigned int j=0; j<symbols.length(); ++j)
			{
				if(line[i]==symbols[j])
				{
					cout << "line " << linecounter << ": unwanted symbol: " << symbols[j] << " at position: " << i << endl << line << endl << "confirm: ";
					getchar();
				}
			}
		}
	linecounter++;
	}
	return;
}

double myHist::returnLowestBin()
{
	return range[0];
}

double myHist::returnHighestBin()
{
	return rangeTOP[sizeOfArrays-1];
}

void myHist::readToBinAssiciatedValues(std::string adress)  //Has to be run AFTER bin-ranges are read in. These also read the amount of bins!
{
	binAssignedValue = new double[sizeOfArrays]; // ranges are 0- x1 -x2 - ... -xn = 1 (usually), that is why i have n ranges and also n average values, where the first average value i put by hand x1 / 2 (" x1 half")s
	std::ifstream f;
	f.open(adress);
	std::string line;
	for(int j = 0 ; j < sizeOfArrays; ++j)
	{
		getline(f,line);
		binAssignedValue[j]=stod(line);
	}
}

void myHist::printInOpenOfstreamSelective(ofstream &f,vector<std::string>  wantedVariables)
{
	for(int i=0; i<sizeOfArrays;++i)
	{
		f << printHistLineSelective(i, wantedVariables) << endl; 
	}
	return;
}



void myHist::printInOpenOfstreamHist(ofstream &f, vector<std::string> option)
{
	if(option[0]=="Full")
	{
		for(int i=0; i<sizeOfArrays;++i)
		{
			f << printHistLine(i) << endl;
		}
	}
	if(option[0]=="select")
	{
		option.erase(option.begin()); // delete first option. Rest options should be the which ones i select.
		
		for(int i=0; i<sizeOfArrays;++i)
		{
			f << printHistLineSelective(i, option) << endl; 
		}
	}
	return;
}

std::string myHist::printHistLineSelective(int &i, vector<std::string> & wantedVariables)
{
	//next line is not good style, but it works
	//may require manual adjustment
	map <std::string, double> m{{"xB",binAssignedValue[i]},{"xmin",range[i]},{"xmax",rangeTOP[i]},{"",normalizedvalue[i]},{"sigmaR",normalizedvalue[i]},{"binCounter",unnormalizedvalue[i]},{"Q2", Q2bin}};
	
	// create line in desired order:	
	std::string line = to_string(m[wantedVariables[0]]);
	for(decltype(wantedVariables.size()) j=1; j < wantedVariables.size(); ++j)
	{
		line += delimiter + to_string(m[wantedVariables[j]]);
	}
	
	#ifdef DEBUGREADGRID
	cout << "output of myHistLineSelective:" << endl;
	cout << line << endl;
	#endif
	return line;
}

void myHist::readCrossSectionNormalization(double factor)
{
	mySigmaNorm = factor;
	return;	
}
void myHist::readBins(std::string adressOfBins)
{
	std::ifstream f;
	f.open(adressOfBins);
	int i = 0;
	std::string line;
	while(getline(f,line))
	{
		++i;
	}
	f.close();
	sizeOfArrays = i;
	range = new double[sizeOfArrays];
	unnormalizedvalue = new double[sizeOfArrays];
	f.open(adressOfBins);
	for(int j = 0 ; j < sizeOfArrays; ++j)
	{
		getline(f,line);
		range[j]=stod(line);
		unnormalizedvalue[j] = 0.;
	}
}

void myHist::createBinsLog(double lowBorder, double highBorder, int binsPerDecade)
{	
	#ifdef DEBUGREADGRID
	cout << "creating log Hist with values: lower Limit: " << lowBorder << "; higher Limit " << highBorder << "; binsPerDecade " << binsPerDecade << endl;
	#endif
	sizeOfArrays =  double(binsPerDecade)*(log10(highBorder/lowBorder))+1;//// bins per decade(manual) * orders of magnitude between xlow and xhigh -1 and THEN +1 to also include xhigh!!
	#ifdef DEBUGREADGRID
	cout << "sizeOfArray " << sizeOfArrays << endl;
	#endif
	unnormalizedvalue = new double[sizeOfArrays];
	range = new double[sizeOfArrays];
	for(int i = 0 ; i < sizeOfArrays; ++i)
	{
		range[i]=pow(10.,log10(lowBorder)+double(i)/binsPerDecade);
		unnormalizedvalue[i] = 0.;
	}
	#ifdef DEBUGREADGRID
	cout << "created log Hist with values: lower Limit: " << lowBorder << "; higher Limit " << highBorder << "; binsPerDecade " << binsPerDecade << endl;
	#endif
	return;
}

void myHist::createBinsLin(double lowBorder, double highBorder, int bins) //CHEck this routine, I think it is not correct at generation. At least when i compare with in-built-Pythia Histogram, i get different results.
{
	#ifdef DEBUGREADGRID
	cout << "creating linear Hist with values: lower Limit: " << lowBorder << "; higher Limit " << highBorder << "; bins " << bins << endl;
	#endif
	
	sizeOfArrays = bins;
	//sizeOfArrays =  int(double(highBorder-lowBorder)/bins+1);//// bins per decade(manual) * orders of magnitude between xlow and xhigh -1 and THEN +1 to also include xhigh!!
	#ifdef DEBUGREADGRID
	cout << "sizeOfArray " << sizeOfArrays << endl;
	#endif
	unnormalizedvalue = new double[sizeOfArrays];
	range = new double[sizeOfArrays];
	for(int i = 0 ; i < sizeOfArrays; ++i)
	{
		range[i]=lowBorder + i*(highBorder-lowBorder)/bins;
		unnormalizedvalue[i] = 0.;
	}
	#ifdef DEBUGREADGRID
	cout << "created linear Hist with values: lower Limit: " << lowBorder << "; higher Limit " << highBorder << "; bins " << bins << endl;
	#endif
	return;
}

void myHist::cantAcceptValue(double &measured)
{
	#ifdef DEBUGREADGRID
	cout << "cannot read this value in: " << measured << ", while range is from " << range[0] << " to " << rangeTOP[sizeOfArrays-1] << endl;
	#endif
	#ifndef DEBUGREADGRID
	double a = measured; 
	a++;// so i dont have warning "measured unused variable if debug is out". 
	#endif
	failedEntries++;
}

void myHist::acceptedValue()
{
	acceptedEntries++;
}

void myHist::addValue(double measured)
{	bool eventAccepted = false;
	for(int i = 0; i < sizeOfArrays; ++i)
	{
		if((measured <= rangeTOP[i]) && (measured > range[i]))
		{
			unnormalizedvalue[i]++; //here one has to be very careful and also compare this with what Python makes out of it, especially when i do a stepfunction plot! So the value is really in that regime where it is associated to. good.
			eventAccepted = true;
		}	
	}
	
	if(!eventAccepted)//if it could not save the value (so far this condition is only whether it fits inside the x Grid)
	{
		cantAcceptValue(measured);		
	}
	else
		acceptedValue();
	
	return;
}

void myHist::addValueWeighted(double value, double measured)
{	bool x = false;
	for(int i = 0; i < sizeOfArrays-1; ++i)
	{
		if(measured<=range[i+1] && measured >range[i])
		{
			unnormalizedvalue[i+1]+=value; //here one has to be very careful and also compare this with what Python makes out of it, especially when i do a stepfunction plot! So the value is really in that regime where it is associated to. good.
			x = true;
		}	
	}
	if(measured<range[0])
		{
			unnormalizedvalue[0]++;
			x = true;
		}
	if(!x)//if it could not save the value
	{
		cantAcceptValue(measured);
	}
	else
		acceptedValue();
	return;
}
void myHist::normalizeValues() // normalize by bin size
{
	if(!valuesAre2ReducedCrossSection)
	{
		cout << "myHist::normalizeValues(): Trying to compute the normalized values, but the values have not been transfered to reduced cross section (which they should be at this step)." << endl << "any key to continue" << endl;
		getchar();
	}
	if(mySigmaNorm==0.)
	{
		cout << "myHist::normalizeValues(): SigmaNorm is 0. \n any key" << endl;
		getchar();
	}
	normalizedvalue = new double[sizeOfArrays];
	for(int i = 0; i<sizeOfArrays; ++i)
	{
		normalizedvalue[i]=mySigmaNorm*1./(rangeTOP[i]-range[i])*1./(Q2binMax-Q2binMin)*intermediatevalue[i];              //
	}
	valuesAreNormalized = true;
	return;
}


void myHist::printHistogramToAdress(std::string adress)
{
	#ifdef DEBUGREADGRID
	cout << "printHistogramToAdress" << endl;
	#endif
	normalizeValues();
	ofstream f;
	f.open(adress);
	for(int i=0; i<sizeOfArrays;++i)
	{
		f << printHistLine(i) << endl;
	}
	f.close();
	#ifdef MANUALSTOP
	cout << "myHist::printHistogramToAdress	above one should see the printed Histogram." << endl;
	getchar();
	#endif
	//createInfo(adress); not needed while there is infoclass.
}

std::string myHist::printHistLine(int &i)
{
	if(!valuesAreNormalized)
	{
		cout << "myHist::printHistLine: Values have not been normalized. Therefore it makes no sense to print them." << endl;
		getchar();
	}
	vector<std::string> allVariables = {"Q2", "xB", "xmin", "xmax", "sigmaR"};
	std::string line = printHistLineSelective(i, allVariables);
	
	//std::string line = to_string(range[i]) + delimiter + to_string(rangeTOP[i]) + delimiter + to_string(binAssignedValue[i]) + delimiter + to_string(unnormalizedvalue[i]) + delimiter + to_string(intermediatevalue[i]) + delimiter + to_string(normalizedvalue[i]);
	#ifdef DEBUGREADGRID
	cout << "myHist::printHistLine	generated outputline is " << line << endl;
	#endif
	return line;
}

void myHist::columnInfo(ofstream &f) //Adapt this according to "printHistLine"!
{
	f << "xbinlow" << delimiter << "xbinhigh" << delimiter << "xbinValue" << delimiter << "EventsInBin" << delimiter <<  "intermediateResult" << delimiter << "sigma_r";
	return;
}

void myHist::printHistogramToAdressNN(std::string adress) // for Debugging unnormalized output
{
	ofstream f;
	f.open(adress);
	for(int i=0; i<sizeOfArrays;++i)
	{
		f << range[i] << "	" << unnormalizedvalue[i] << endl; 
	}
	f.close();
}

void myHist::createInfo(std::string &adress)
{
	ofstream f;
	std::string infoAdress = adress + infoString;
	f.open(infoAdress);
	
	f << "Normalization factor	" << mySigmaNorm << endl;
	f << "Number of measurements that were not accepted:	" << failedEntries << endl;
	f << "Number of measurements that were accepted:	" << acceptedEntries << endl;
	
	f.close();
}


void myHist::printHistogramToTerminal()
{
	normalizeValues();
	for(int i=0; i<sizeOfArrays;++i)
	{
		cout << printHistLine(i) << endl;
	}
}

void myHist::fullyNormalize(const double sigmaNormR,const double sMandelstam,double& prefactor)
{
	readCrossSectionNormalization(sigmaNormR);
	crossSection2ReducedCrossSection(sMandelstam, prefactor);
	normalizeValues();
	return;
}

void myHist::crossSection2ReducedCrossSection(const double s, double& constFactor) //VERY IMPORTANT, this function determines output value by computing sigma reduced (from d^2 sigma/ dx dQ2 (x,Q2,s))
{
	if(valuesAre2ReducedCrossSection)
	{
		cout << "myHist::crossSection2ReducedCrossSection: this has already been executed." << endl;
		getchar();
		return;
	}
	for(int i=0; i< sizeOfArrays; ++i)
	{	
		//check for possible errors:
		if(s==0 || constFactor ==0)
		{
			cout << "crossSecton2ReducedCrossSection::Warning:	The used values are s=" << s << " and prefactor=" << constFactor << endl << "press any key to continue" << endl;
			getchar();
		}
	
		//const factor includes pi and alpha
		double Q2 = Q2bin;
		double x = binAssignedValue [i];
		double y = Q2/(x*s); // if using s = 2 pProton*pElectron this is exact. (because the error compensates)
		double Y_Plus = 1.+(1.-y)*(1.-y);
		intermediatevalue[i]=constFactor*Q2*Q2*x/(Y_Plus)*unnormalizedvalue[i];
	}
	valuesAre2ReducedCrossSection = true;
	return;
}

void myHist::readBinsPLUSasVal(std::string adress, const double &Q2lowlimit, const double &Q2uplimit, const double &Q2binVal)
{
	//first check if there are symbols in the file put by accident which we dont want:
	
	checkSourceFileforsymbol(unwantedsymbols, adress);
	
	//First find out how many bins I have here.(to know how big the Arrays have to be) Therefore read and check how many matches I have for the Q2 regime
	
	//Definitions:
	std::ifstream f;
	int counter = 0;
	std::string line;
	int dataperline = 8; // make this better
	std::string *entry = new std::string[dataperline];
	
		
	#ifdef DEBUGREADGRID
	cout << "adress of Grid: " << adress << " ...trying to open" << endl;
	#ifdef MANUALSTOP
	getchar();
	#endif
	#endif 
	
	f.open(adress);
	if(f.is_open())
	{
		#ifdef DEBUGREADGRID
		cout << "adress should be open now. Test entries: Lines are read out:" << endl;
		#ifdef MANUALSTOP
		getchar();
		#endif
		#endif
		
		while(getline(f,line))
		{	
			#ifdef DEBUGREADGRID
			cout << "whole line	" << line << endl;
			#endif
			for(int i=0; i < dataperline; ++i)
			{
				size_t pos;
				if((i+1)<dataperline)
				{
					pos = line.find(delimiter);
					entry[i]= line.substr(0,pos);
					line = line.substr(pos+1,line.length()-(pos+1));
				}
				else{entry[i] = line;}
				#ifdef DEBUGREADGRID
				cout << "entry[" << i << "]	" << entry[i] << "	current line std::string:	" << line << endl;
				cout << "pos: " << pos << endl;
				#endif
				
			}
			
			if((Q2lowlimit==stod(entry[0])) && (Q2uplimit==stod(entry[2])) && (Q2binVal==stod(entry[3]))) // entry[1] should be a "-" (Q2min	-	Q2max	Q2...) 
			{
				++counter;
				#ifdef DEBUGREADGRID
				cout << counter << " data with Q2 energy matching to " << Q2lowlimit << ";" << Q2uplimit << ";" << Q2binVal << endl;
				#ifdef MANUALSTOP
				getchar();
				#endif
				#endif
			}
		}
		
		f.close();
	}
	else
	{
		cout << "cannot open file with adress: " << adress << endl;
		getchar();
		f.close();
	}
	#ifdef DEBUGREADGRID
	cout << "std::ifstream closed now" << endl;
	#ifdef MANUALSTOP
	getchar();
	#endif
	#endif
	
	//now array size is known. Create arrays, read the file again and now read in the values.
	sizeOfArrays = counter;
	#ifdef DEBUGREADGRID
	cout << sizeOfArrays << "	size of Arrays" << endl;
	#ifdef MANUALSTOP
	getchar();
	#endif
	#endif
	counter = 0;
	
	range = new double[sizeOfArrays];
	rangeTOP = new double[sizeOfArrays];
	unnormalizedvalue = new double[sizeOfArrays];
	intermediatevalue = new double[sizeOfArrays];
	normalizedvalue = new double[sizeOfArrays];
	binAssignedValue = new double[sizeOfArrays];
	f.open(adress);
	while(getline(f,line))
	{
		#ifdef DEBUGREADGRID //HARD COPIED FROM ABOVE
		cout << "whole line	" << line << endl;
		#endif
		for(int i=0; i < dataperline; ++i)
		{
			size_t pos;
			if((i+1)<dataperline)
			{
				pos = line.find(delimiter);
				entry[i]= line.substr(0,pos);
				line = line.substr(pos+1,line.length()-(pos+1));
			}
			else{entry[i] = line;}
			#ifdef DEBUGREADGRID
			cout << "entry[" << i << "]	" << entry[i] << "	current line std::string:	" << line << endl;
			cout << "pos: " << pos << endl;
			#endif
			
		}
		
		#ifdef DEBUGREADGRID 
		cout << stod(entry[0])-Q2lowlimit << "	" << stod(entry[2])-Q2uplimit << "	" << stod(entry[3])-Q2binVal << endl;
		#endif
		
		if((Q2lowlimit==stod(entry[0])) && (Q2uplimit==stod(entry[2])) && (Q2binVal==stod(entry[3]))) // entry[1] should be a "-" (Q2min	-	Q2max	Q2...) 
		{
			#ifdef DEBUGREADGRID
			cout << "MATCHING VALUES, ADDING BIN" << endl;
			
			for(int i=0; i<dataperline; ++i)
			{
				if((i!=1)&&(i!=5))
				cout << stod(entry[i]) << "	as double, and std::string was " << entry[i] << endl;
			}
			#endif
			
			range[counter] = stod(entry[4]);
			rangeTOP[counter] = stod(entry[6]);
			binAssignedValue[counter] = stod(entry[7]);
			unnormalizedvalue[counter] = 0.;
			normalizedvalue[counter] = 0.;
			intermediatevalue[counter] = 0.;
			
			++counter;
		}
	}
	f.close();
	return;
}

