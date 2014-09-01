
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

#include <iostream>
#include <R_ext/Rdynload.h>
using namespace std;


bool is_overlap(double x1, double y1, double sw1, double sh1, NumericMatrix boxes){
    double x2, y2, sw2, sh2;
    bool overlap = false;
    int nElements = boxes.nrow();
    for(int i=0; i < nElements; i++){
        x2 = boxes(i, 0);
        y2 = boxes(i, 1);
        sw2 = boxes(i, 2);
        sh2 = boxes(i, 3);
        
        if(x2 == -1 && y2 == -1){
            break;
        }
        if(x2 == 2 && y2 == 2){
            return true;
        }
        
        if (x1 < x2){
            overlap = (x1 + sw1) > x2;
        }
        else{
            overlap = (x2 + sw2) > x1;
        }
        
        if (y1 < y2){
            overlap = (overlap && ((y1 + sh1) > y2));
        }
        else{
            overlap = (overlap && ((y2 + sh2) > y1));
        }
        if(overlap){
            return true;
        }
    }
    return false;
}

// definition

extern "C" SEXP findCoordinates( SEXP width, SEXP height ){
    BEGIN_RCPP
    Rcpp::NumericVector width1(width);
    Rcpp::NumericVector height1(height);
    Rcpp::NumericMatrix boxes(width1.size(), 4);
    boxes.fill(-1);

    Environment stats("package:stats");
    Function runif = stats["runif"];
    Environment base("package:base");
    Function cos = base["cos"];
    Function sin = base["sin"];

    double x, y, r, theta;
    double thetaStep = 0.1;
    double rStep = 0.05;


    for(int i=0; i<width1.size(); i++){
        x = 0.5;
        y = 0.5;
        r = 0;
    
        double rand = as<double>(runif(1));
        theta = rand * 2 * 3.141593;
        bool isOverlapped = true;
        while(isOverlapped){
            if(!is_overlap(x - 0.5 * width1(i), y - 0.5 * height1(i), width1(i), height1(i), boxes) && (x - 0.5 * width1(i) > 0) && (y - 0.5 * height1(i) > 0) && (x + 0.5 * width1(i) < 1)  && (y + 0.5 * height1(i) < 1)){
                boxes(i, 0) = x - 0.5 * width1(i);
                boxes(i, 1) = y - 0.5 * height1(i);
                boxes(i, 2) = width1(i);
                boxes(i, 3) = height1(i);
                isOverlapped = false;
            }
            else{
                if(r > 0.7071068){
                    boxes(i, 0) = 2;
                    boxes(i, 1) = 2;
                    boxes(i, 2) = 2;
                    boxes(i, 3) = 2;
                
                    isOverlapped = false;
                }
                else{
                    theta += thetaStep;
                    r += rStep * thetaStep / (2 * 3.141593);
                    x = 0.5 + r * as<double>(cos(theta));
                    y = 0.5 + r * as<double>(sin(theta));
                }
            }
        }    
    }

    boxes(_, 0) = boxes(_, 0) + 0.5 * boxes(_, 2);
    boxes(_, 1) = boxes(_, 1) + 0.5 * boxes(_, 3);

    return Rcpp::wrap(boxes);
    END_RCPP
}




// definition

extern "C" SEXP findCoordinates_left( SEXP width, SEXP height ){
    BEGIN_RCPP
    Rcpp::NumericVector width1(width);
    Rcpp::NumericVector height1(height);
    Rcpp::NumericMatrix boxes(width1.size(), 4);
    boxes.fill(-1);

    Environment stats("package:stats");
    Function runif = stats["runif"];
    Environment base("package:base");

    double x, y, sign, sign2;
    double vertStep = 0.02;
    double horizStep = 0.02;
    
    double rand = as<double>(runif(1));
    sign = rand - 0.5 ;


    for(int i=0; i<width1.size(); i++){
        x = 0.01;
        y = 0.5;
        sign = -sign;
        sign2 = sign;
        
        bool isOverlapped = true;
        while(isOverlapped){
            if(!is_overlap(x, y - 0.5 * height1(i), width1(i), height1(i), boxes) && (x > 0) && (y - 0.5 * height1(i) > 0) && (x + width1(i) < 1)  && (y + 0.5 * height1(i) < 1)){
                boxes(i, 0) = x;
                boxes(i, 1) = y - 0.5 * height1(i);
                boxes(i, 2) = width1(i);
                boxes(i, 3) = height1(i);
                isOverlapped = false;
            }
            else{
                if(x > 1){
                    boxes(i, 0) = 2;
                    boxes(i, 1) = 2;
                    boxes(i, 2) = 2;
                    boxes(i, 3) = 2;

                    isOverlapped = false;
                }
                else{
                    if(sign2 < 0){
                        y = y - vertStep;
                    }
                    else{
                        y = y + vertStep;
                    }

                    if((y > 1) || (y < -1)){
                        y = 0.5;
                        x = x + horizStep;
                        sign2 = -sign2;
                    }
                }
            }
        }    
    }

    boxes(_, 0) = boxes(_, 0);
    boxes(_, 1) = boxes(_, 1) + 0.5 * boxes(_, 3);

    return Rcpp::wrap(boxes);
    END_RCPP
}


// definition

extern "C" SEXP findCoordinates_left_top( SEXP width, SEXP height ){
    BEGIN_RCPP
    Rcpp::NumericVector width1(width);
    Rcpp::NumericVector height1(height);
    Rcpp::NumericMatrix boxes(width1.size(), 4);
    boxes.fill(-1);

    Environment stats("package:stats");
    Function runif = stats["runif"];
    Environment base("package:base");

    double x, y;
    double vertStep = 0.02;
    double horizStep = 0.02;

    for(int i=0; i<width1.size(); i++){
        x = 0.01;
        y = 1;
        
        bool isOverlapped = true;
        while(isOverlapped){
            if(!is_overlap(x, y - 0.5 * height1(i), width1(i), height1(i), boxes) && (x > 0) && (y - 0.5 * height1(i) > 0) && (x + width1(i) < 1)  && (y + 0.5 * height1(i) < 1)){
                boxes(i, 0) = x;
                boxes(i, 1) = y - 0.5 * height1(i);
                boxes(i, 2) = width1(i);
                boxes(i, 3) = height1(i);
                isOverlapped = false;
            }
            else{
                if(x > 1){
                    boxes(i, 0) = 2;
                    boxes(i, 1) = 2;
                    boxes(i, 2) = 2;
                    boxes(i, 3) = 2;

                    isOverlapped = false;
                }
                else{
                    y = y - vertStep;
                    
                    if((y < -1)){
                        y = 1;
                        x = x + horizStep;
                    }
                }
            }
        }    
    }

    boxes(_, 0) = boxes(_, 0);
    boxes(_, 1) = boxes(_, 1) + 0.5 * boxes(_, 3);

    return Rcpp::wrap(boxes);
    END_RCPP
}

// register methods
R_CallMethodDef CallMethods[] = {
    {"findCoordinates", (DL_FUNC) &findCoordinates, 2},
    {"findCoordinates_left", (DL_FUNC) &findCoordinates_left, 2},
    {"findCoordinates_left_top", (DL_FUNC) &findCoordinates_left_top, 2}
};

void R_init_GOsummaries(DllInfo *info){
    R_registerRoutines(info, NULL, CallMethods, NULL, NULL);
}
