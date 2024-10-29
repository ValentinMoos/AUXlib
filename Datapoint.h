#ifndef DATAPOINT
#define DATAPOINT

#include"Aux.h"

class Datapoint{
	private:
		char delimiter = ',';
		unsigned nSlots = 26;
		std::string point_id;
		
		int process_id;
		int weight_process;
		
		double s; //GeV^2
		double targetMass; //GeV
		double productMass; //GeV
		double Q2, Q2min, Q2max; //GeV
		double xB, xBmin, xBmax;
		double z, zmin, zmax;
		double pT, pTmin, pTmax; //GeV
		
		double ycut_min = 0., ycut_max = 1.;
		double W2cut_min = -1., W2cut_max  = -1.; //GeV^2
		
		double xSec; //this typically also has a unit (barn)
		double Uncorr_Err0;
		double Th_Factor;
		
		bool FiducialCuts;
		
		
		std::string present_process_id() const;
		std::string present_weight_process() const;
		//std::string present_Bool(bool);

		
		void compute_Th_Factor();
		
		
		
	public:
		
		//return
		std::string give_point_id();
		
		int give_process_id();
		int give_weight_process();
		
		double give_s();
		double give_targetMass();
		double give_productMass();
		
		double give_Q2();
		double give_Q2min();
		double give_Q2max();
		
		double give_xB();
		double give_xBmin();
		double give_xBmax();
		
		double give_z();
		double give_zmin();
		double give_zmax();
		
		double give_pT();
		double give_pTmin();
		double give_pTmax();
		
		double give_ymin();
		double give_ymax();

		double give_W2min();
		double give_W2max();
		
		double give_xSec();
		double give_Uncorr_Err0();
		double give_Th_Factor();
		
		bool give_FiducialCuts();
		
		
		//set values
		void set_point_id(std::string);
		
		void set_process_id(int);
		void set_weight_process(int);
		
		void set_s(double);
		void set_targetMass(double);
		void set_productMass(double);
		
		void set_Q2(double);
		void set_Q2min(double);
		void set_Q2max(double);
		
		void set_xB(double);
		void set_xBmin(double);
		void set_xBmax(double);
		
		void set_z(double);
		void set_zmin(double);
		void set_zmax(double);
		
		void set_pT(double);
		void set_pTmin(double);
		void set_pTmax(double);
		
		void set_ymin(double);
		void set_ymax(double);

		void set_W2min(double);
		void set_W2max(double);
		
		void set_xSec(double);
		void set_Uncorr_Err0(double);
		void set_Th_Factor(double);
		
		void set_FiducialCuts(bool);
		
		
		std::string printform() const;
		void write(std::ofstream&);
		
		//std::ostream& operator<<(std::ostream& os);
	
};

std::ostream& operator<<(std::ostream&, const Datapoint&);

#endif
