#include"num_int.h"


// first example integrate a real function of single real parameter.
double simpson1(double (*f)(double), const double minCoordinate, const double maxCoordinate, const double precision, std::vector<double> &fvals, std::vector<double> &xvals, long unsigned int &counter){


    //const int minEvaluations = 10;
    //const int maxEvaluations = 100;

    double intval = 0;
    
    double I0, I1, I2, I3;
    float precision_reached = false;

    double x0 = minCoordinate;
    double x2 = maxCoordinate;
    double x1;

    double f0, f1, f2;

    if(fvals.size() == 0){ // this means this is the first iteration of the function.
        
        f0 = f(x0);
        f2 = f(x2);

        fvals.push_back(f0);
        xvals.push_back(x0);
        
        fvals.push_back(f2);
        xvals.push_back(x2);
        counter += 2;

    }
    else{ // values are already in list.
        f0 = fvals[getIndex(xvals, x0)];
        f2 = fvals[getIndex(xvals, x2)];
    }
/* 
    if(counter > 30){
        std::cout << xvals << std::endl;
        //exit(1);
    } */




    //std::vector<double> evaluated_coordinates {xl, xu};
    //std::vector<double> evaluated_fvalues {fl, fu};

    //one step where we do need to evaluate the function:
    x1 = (x0 + x2)/2;
    if(onList(x1, xvals)){
        f1 = fvals[getIndex(xvals, x1)];
    }
    else{
        ++counter; // count calls of fucntion
        f1 = f(x1);
    }
    I0 = (f0+f2)*(x2-x0)/2;
    I1 = (f0+f1)*(x1-x0)/2;
    I2 = (f1+f2)*(x2-x1)/2;
    I3 = I1 + I2;

    if( abs(I3 - I0) < precision){
        
        //complete.
        //std::cout << "precision reached in segment " << x0 << " " << x2 << " and returns " << I0 << std::endl << "intersegments " << x0 << " " << x1 << "     returns " << I1 << std::endl << "              " << x1 << " " << x2 << "     returns " << I2 << std::endl;
        //std::cout << "ratio I3/I0 is " << abs(I3/I0) << " and precision demanded was " << precision << std::endl;
        return I0;
    }

    else{
        //this means we need to split the integral in two, and compute each section individually. Before we add the gained info: the middle point:
        xvals.push_back(x1);
        fvals.push_back(f1);
        I1 = simpson1(f, x0, x1, precision, fvals, xvals, counter); // left part
        I2 = simpson1(f, x1, x2, precision, fvals, xvals, counter); // right part

        //if both functions return here it means both sides could be computed with the demanded precision.
        I0 = I1 + I2;
        return I0;
    }

/* 


    double fm;

    for(int i = 0; precision_reached; ++i){

        xm = (xl + xu) / 2;
        //function values
        //fl = f(xl);
        //fu = f(xu);
        fm = f(xm);


        I = (fl + fu)*(xl-xu);
        Il = (fl + fm)*(xm-xl);
        Iu = (fu + fm)*(xu-xm);
        I1 = Il + Iu;
         

        if( abs(Inext / I - 1) < precision){
            precision_reached = true;
        }
        else{ // prepare next iteration.
            //add new insights..
            evaluated_coordinated.append(xm);
            evaluated_fvalues.append(fm);



        }    


    } */



}




std::vector<double> accept_reject(double (*f)(double), const double xmin, const double xmax, const double ymin, const double ymax, const long unsigned int N){
    std::random_device rd_x; // random number
    std::random_device rd_y; // random number
    std::cout << rd_x() << std::endl;
    std::cout << rd_y() << std::endl;

    std::mt19937 gen_x(rd_x());
    std::mt19937 gen_y(rd_y());

    std::uniform_real_distribution<> uniform_x(xmin, xmax); // uniform distribution in this range.
    std::uniform_real_distribution<> uniform_y(ymin, ymax); // uniform distribution in this range.

    std::vector<double> distribution;

    for(long unsigned int i = 0; i < N; ++i){
        double x = uniform_x(gen_x);
        double y = uniform_y(gen_y);
        //std::cout << x << " " << y << " " << f(x) << std::endl;
        if(y < f(x)){
            distribution.push_back(x);
        }
        //else : reject.
        else{
            --i; // such that i obtain N in total.
        }
    
    }

    return distribution;    
}