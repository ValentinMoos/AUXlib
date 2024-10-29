#include "Datapoint.h"


std::string Datapoint::present_process_id() const{
	std::string str = "1 - 1 - " + std::to_string(process_id);
	return str;
}

std::string Datapoint::present_weight_process() const{
	std::string str = "1 - 1- " + std::to_string(weight_process);
	return str;
}

/*std::string Datapoint::present_Bool(bool x){
	return TrueFalse(x);
}*/

std::string Datapoint::printform() const{

	std::vector<std::string> v(nSlots,"");
	
	v[0]  = point_id;
	v[1]  = present_process_id();
	v[2]  = present_weight_process();
	v[3]  = std::to_string(s);
	v[4]  = std::to_string(Q2);
	v[5]  = std::to_string(Q2min);
	v[6]  = std::to_string(Q2max);
	v[7]  = std::to_string(xB);
	v[8]  = std::to_string(xBmin);
	v[9]  = std::to_string(xBmax);
	v[10] = std::to_string(z);
	v[11] = std::to_string(zmin);
	v[12] = std::to_string(zmax);
	v[13] = std::to_string(pT);
	v[14] = std::to_string(pTmin);
	v[15] = std::to_string(pTmax);
	v[16] = std::to_string(xSec);
	v[17] = std::to_string(Uncorr_Err0);
	v[18] = std::to_string(Th_Factor);
	v[19] = TrueFalse(FiducialCuts);
	v[20] = std::to_string(ycut_min);
	v[21] = std::to_string(ycut_max);
	v[22] = std::to_string(W2cut_min);
	v[23] = std::to_string(W2cut_max);
	v[24] = std::to_string(targetMass);
	v[25] = std::to_string(productMass);
	
	std::string line = vec_to_line(v, delimiter);
	return line;
	
}


void Datapoint::write(std::ofstream& f){
	//std::string line = printform();
	f << printform() << std::endl;
	
}

void Datapoint::compute_Th_Factor(){
	Th_Factor = 1.;
}

void Datapoint::set_point_id(std::string s){
	point_id = s;
}

void Datapoint::set_process_id(int q){
	process_id = q;
}

void Datapoint::set_weight_process(int q){
	weight_process = q;
}

void Datapoint::set_s(double q){
	s = q;
}
void Datapoint::set_targetMass(double q){
	targetMass = q;
}

void Datapoint::set_productMass(double q){
	productMass = q;
}
void Datapoint::set_Q2(double q){
	Q2 = q;
}
void Datapoint::set_Q2min(double q){
	Q2min = q;
}
void Datapoint::set_Q2max(double q){
	Q2max = q;
}

void Datapoint::set_xB(double q){
	xB = q;
}
void Datapoint::set_xBmin(double q){
	xBmin = q;
}
void Datapoint::set_xBmax(double q){
	xBmax = q;
}

void Datapoint::set_z(double q){
	z = q;
}
void Datapoint::set_zmin(double q){
	zmin = q;
}
void Datapoint::set_zmax(double q){
	zmax = q;
}

void Datapoint::set_pT(double q){
	pT = q;
}
void Datapoint::set_pTmin(double q){
	pTmin = q;
}
void Datapoint::set_pTmax(double q){
	pTmax = q;
}

void Datapoint::set_ymin(double q){
	ycut_min = q;
}
void Datapoint::set_ymax(double q){
	ycut_max = q;
}

void Datapoint::set_W2min(double q){
	W2cut_min = q;
}
void Datapoint::set_W2max(double q){
	W2cut_max = q;
}

void Datapoint::set_xSec(double q){
	xSec = q;
}
void Datapoint::set_Uncorr_Err0(double q){
	Uncorr_Err0 = q;
}
void Datapoint::set_Th_Factor(double q){
	Th_Factor = q;
}

void Datapoint::set_FiducialCuts(bool q){
	FiducialCuts = q;
}
////#####################################################


std::string Datapoint::give_point_id(){
	return point_id;
}

int Datapoint::give_process_id(){
	return process_id;
}

int Datapoint::give_weight_process(){
	return weight_process;
}

double Datapoint::give_s(){
	return s;
}
double Datapoint::give_targetMass(){
	return targetMass;
}

double Datapoint::give_productMass(){
	return productMass;
}

double Datapoint::give_Q2(){
	return Q2;
}
double Datapoint::give_Q2min(){
	return Q2min;
}
double Datapoint::give_Q2max(){
	return Q2max;
}

double Datapoint::give_xB(){
	return xB;
}
double Datapoint::give_xBmin(){
	return xBmin;
}
double Datapoint::give_xBmax(){
	return xBmax;
}

double Datapoint::give_z(){
	return z;
}
double Datapoint::give_zmin(){
	return zmin;
}
double Datapoint::give_zmax(){
	return zmax;
}

double Datapoint::give_pT(){
	return pT;
}
double Datapoint::give_pTmin(){
	return pTmin;
}
double Datapoint::give_pTmax(){
	return pTmax;
}

double Datapoint::give_ymin(){
	return ycut_min;
}
double Datapoint::give_ymax(){
	return ycut_max;
}

double Datapoint::give_W2min(){
	return W2cut_min;
}
double Datapoint::give_W2max(){
	return W2cut_max;
}

double Datapoint::give_xSec(){
	return xSec;
}
double Datapoint::give_Uncorr_Err0(){
	return Uncorr_Err0;
}
double Datapoint::give_Th_Factor(){
	return Th_Factor;
}

bool Datapoint::give_FiducialCuts(){
	return FiducialCuts;
}


//OPERATOROVERLOADING
std::ostream& operator<<(std::ostream& os, const Datapoint& obj){
	os << obj.printform();
	
	return os;
}


