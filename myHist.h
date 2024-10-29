#ifndef MYHIST
#define MYHIST

/*
Changed on 22.03.2021
features: 
	Designed for logarithmic DIS Sigma_reduced-plots/histograms. CAUTION if use for other purpose.
	read Histogram (ranges) in
	create Logarithmic Histogram
	create Linear Histogram (does not work properly, has to be checked)
	Unlike previous version this works with two ranges: range[] and rangeTOP[], 
	Use constructors to initialize values. To avoid user( that is, my own) errors
	Normalization is Full, also in Q2
	
*/
class myHist{
	private:
	
		//Variables
		int sizeOfArrays;
		double *range;
		double *rangeTOP; 		   // for reading bins like xmin - xmax; xmin - max (and it could be that xmax_i isnot xmin_i+1) 
								   // I consider it useful to define second range. Possibly it is not needed, but i will work with this in any case. 
		double *unnormalizedvalue; // counter for "events in bin"
		double *intermediatevalue; // intermediate computation value (such that original counter can be stored)
		double *normalizedvalue;   // final value for sigma_Reduced
		double *binAssignedValue;  // x - value for the xbin
		double Q2bin , Q2binMin, Q2binMax;
		double mySigmaNorm = 0.;
		int failedEntries = 0;
		int acceptedEntries = 0;
		std::string infoString = "_info";
		const char delimiter = '	'; //TAB
		std::string unwantedsymbols = " ,"; // space, komma
		
		bool valuesAreNormalized = false;
		bool valuesAre2ReducedCrossSection = false;
		
		//Create Bins:
		void createBinsLog(double lowBorder, double highBorder,int binsPerDecade);           // creates 10-log bin
		void createBinsLin(double lowBorder, double highBorder,int bins);                    // creates linear bin. Not sure if its correct. Has to be checked.
		void readBins(std::string adressOfBins);                                                  // reads in bins from document
		void readToBinAssiciatedValues(std::string adress);  
		void readBinsPLUSasVal(std::string adress, const double &Q2low, const double &Q2high, const double &Q2bin); // reads in whole histogram structure (depending on Q2 bin too!)                                      // reads bin mean value from document
		
		//Intrinsic computation: normalization,..
		void normalizeValues();                                                              //divide by this bin size. But only the bin size of the bins "here". not by Q2 bin, only x (in "standard" usage)
		void acceptedValue();                                                                // raises counter for acceptedValues
		void crossSection2ReducedCrossSection(const double s, double& constFactor);          //Does only make sense for xB; multiplies values by factor(x,Q2,s,alpha_em,y)
		
		//Error prevention
		void checkInitialisation();     													 // checks if initialization was ok
		void checkSourceFileforsymbol(std::string symbols, std::string adress);                        // checks if there are hidden space, kommas whatevar 
																							 // (specified in "unwantedsymbols-std::string") and one does not see it.
																							 // This function should prevent error coming from this source                                  
		//Info and Output:     
		std::string printHistLine(int &i);                                                        // for output, simplifies life. returns the ith line of the histogram   (full Line)
		std::string printHistLineSelective(int &i, vector<std::string> & wantedVariables);             // allows to select which variables one wants to print
		void createInfo(std::string &adress);													 // create info file (obsolete since this should be done by Infoclass)
		void printMembers();                                                                 // print all members (variables)
		template<class TYPE>
		void printArray(const int &size,TYPE *array);
		void cantAcceptValue(double &measured);                                              // for debugging. Message when value can not be read into Histogram
																							 // (e.g. when its not in Histogram range)
		
		//gain factors: (called inside public function:)
		void readCrossSectionNormalization(double factor);                                   // reads in the const factor needed for reduced cross section
		
	public:
		// constructur
		myHist(const std::string type, const std::string &gridAdress, const double &Q2, const double &Q2min, const double &Q2max); //constructor; creates Grid
		
		// to add a value which raises the counters
		void addValue(double measured);                              // raises counter for respective bin
		void addValueWeighted(double value, double measured);        // adds a value in the to the respective bin (instead of just "+1")
		
		// not needed
		void printHistogramToAdress(std::string adress);                   //does what it says.    
		void printHistogramToAdressNN(std::string adress);                 //NN means non-normalized. Does what it says. 
		
		// call Normalization
		void fullyNormalize(const double sigmaNormR,const double sMandelstam,double& prefactor);          //to reduced cross section and normalize wrt. bins and event counter
		
		// print results                         // for info class 
		void printInOpenOfstreamSelective(ofstream &f, vector<std::string> wantedValues);         // for info class
		void printInOpenOfstreamHist(ofstream &f, vector<std::string> option);                    // for info class
		void columnInfo(ofstream &f);                                                        // prints the meaning of the columns in the histogram
		void printHistogramToTerminal();                                                     // does what it says
		
		// for info checks
		double returnLowestBin();                                                            // returns the lower limit for measured values (the lowest bin limit)
		double returnHighestBin();                                                           // -""- the hightest bin limit
		void printCounterOfBin(int i);                                                       // for info and debug purposes
		void printValueOfBin(int i);                                                       // for info and debug purposes
		
		//Public functions until tests finished, then shift them to private:
		
};

#endif
