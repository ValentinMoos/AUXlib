#include"statlib.h"

std::vector<double> onesigmastats(const std::vector<double>& distribution){
	const double one_sigma = 0.682689492; // hard copy.
	double average = 0, dev_m, dev_p;
	
	unsigned int l = distribution.size();
	
	std::vector<double> copy1 = order_vector(distribution, '<'); //order vector increasingly
	
	const int lowerlim = floor(l * (1-one_sigma)/2 )+1; // such that index in INSIDE confidence region
    const int upperlim = floor(l * (1+one_sigma)/2 ); // such that index in INSIDE confidence region
    
    dev_m = copy1[lowerlim];
    dev_p = copy1[upperlim];
    
    for(double v : distribution){
    average += v;
    }
	average /= l;
	
	std::vector<double> av_min_pls {average, dev_m, dev_p};
	
	return av_min_pls;
	
}
