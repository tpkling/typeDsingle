//
// Created by Thomas Kling on 11/29/18.
//

#include "ray.h"
#include "Coords.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

const double b21 = 1.0/5.0;
const double b31 = 3.0/40.0;
const double b41 = 3.0/10.0;
const double b51 = -11.0/54.0;
const double b61 = 1631.0/55296.0;

const double b32 = 9.0/40.0;
const double b42 = -9.0/10.0;
const double b52 = 5.0/2.0;
const double b62 = 175.0/512.0;

const double b43 = 6.0/5.0;
const double b53 = -70.0/27.0;
const double b63 = 575.0/13824.0;

const double b54 = 35.0/27.0;
const double b64 = 44275.0/110592.0;

const double b65 = 253.0/4096.0;

const double c1  = 37.0/378.0;
const double c2  = 0.0;
const double c3  = 250.0/621.0;
const double c4  = 125.0/594.0;
const double c5  = 0.0;
const double c6  = 512.0/1771.0;

const double c1s = 2825.0/27648.0;
const double c2s = 0.0;
const double c3s = 18575.0/48384.0;
const double c4s = 13525.0/55296.0;
const double c5s = 277.0/14336.0;
const double c6s = 1.0/4.0;

const double EPS = 1.0e-4; // values were about 1e-3
const double SAFTEY = 0.90;


ray::ray(int i, int N) {

    h_ = 0.0001;

    coords_.Pt_ = 1.0;
    coords_.x_ = 0.0;
    coords_.y_ = 0.0;
    coords_.z_ = 0.0;
    coords_.alpha_ = 0.3;

    double py1 = Py_lim1(coords_.Pt_, coords_.alpha_, coords_.x_);
    double py2 = Py_lim2(coords_.Pt_, coords_.alpha_, coords_.x_);

    coords_.Py_ = py1 + float(i)*(py2-py1)/float(N);
    coords_.Pz_ = 0.0;
    // double pzlim = Pz_lim(coords_.Pt_, coords_.alpha_, coords_.x_, coords_.Py_);
    // coords_.Pz_ = 0.5*pzlim;

    coords_.vx_ = PxSet(coords_.Pt_, coords_.alpha_, coords_.x_, coords_.Py_, coords_.Pz_);
    coords_.dt_ = 0.1; // setting as a dummy positive value at the start
    coords_.lambda_ = 0.0;

}

ray::~ray() {
    //default destructor
}

double ray::Py_lim1(double pt, double alpha, double x) {
    return - pt - alpha * x * pt;
}

double ray::Py_lim2(double pt, double alpha, double x) {
    return + pt - alpha * x * pt;
}

double ray::Pz_lim(double pt, double alpha, double x, double py){
    return sqrt(pt*pt*(1.0 - alpha*alpha * x*x) - py*py -2.0*alpha * x * pt*py);
}

double ray::PxSet(double pt, double alpha, double x, double py, double pz){
    return sqrt(pt*pt*(1.0-alpha*alpha*x*x) - py*py - pz*pz - 2.0*alpha*x*py*pt);
}

void ray::step() {
    Coords newCoords;
    double err[5];

    computeRay(coords_, newCoords, err, h_);
    double changeDir = newCoords.dt_ * coords_.dt_; // if negative, time direction flipped

    double errmax = 0.0;

    for (int j = 0; j<5; j++){
        if (fabs(err[j])>errmax){
            errmax = fabs(err[j]);
        } // this block picks out the biggest error value (in absolute value)
    }

    double errRatio = EPS/errmax;

    if(errRatio > 1) {
        // GOOD STEP
        goodStepCount_++;
        coords_ = newCoords;

        totalError_ = totalError_+  errmax;
        double h1 = SAFTEY * h_ * pow(errRatio,0.25); // from old NRC book
        double h2 = h_ * 1.1; // decreased from 10*h on Sept 9 to limit step size growth (Further on Sept 16th to debug)
        h_ = (h1 > h2) ? h2 : h1;
        if (h_>0.005) {h_ = 0.005;}

        coords_.lambda_ = coords_.lambda_ + h_;


    } else {
        //BAD STEP
        badStepCount_++;
        h_ = SAFTEY * h_ * pow(errRatio, 0.2); // .2 from old NRC book
        if (h_ <= 5.0e-6) {
            h_ = 5.0e-6;
        } // 5e-6

    }
} // closes step

void ray::computeRay(Coords initCoords, Coords &newCoords, double *err, double h) {

    double t_0, x_0, y_0, z_0, vx_0;
    double alpha, pt, py, pz;
    double k1t, k1x, k1y, k1z, k1vx;
    double k2t, k2x, k2y, k2z, k2vx;
    double k3t, k3x, k3y, k3z, k3vx;
    double k4t, k4x, k4y, k4z, k4vx;
    double k5t, k5x, k5y, k5z, k5vx;
    double k6t, k6x, k6y, k6z, k6vx;

    double midx, midvx;

    t_0		  = initCoords.t_;
    x_0       = initCoords.x_;
    y_0       = initCoords.y_;
    z_0		  = initCoords.z_;
    vx_0	  = initCoords.vx_;
    pt        = initCoords.Pt_;
    py        = initCoords.Py_;
    pz        = initCoords.Pz_;
    alpha     = initCoords.alpha_;

    k1t      = h * ft(pt, py, alpha, x_0);
    k1x      = h * vx_0;
    k1y      = h * fy(pt, py, alpha, x_0);
    k1z      = h * pz;
    k1vx     = h * vx_dot(pt, py, alpha, x_0);

    midx  = x_0  + b21*k1x;
    midvx = vx_0 + b21*k1vx;

    k2t      = h * ft(pt, py, alpha, midx);
    k2x      = h * midvx;
    k2y      = h * fy(pt, py, alpha, midx);
    k2z      = h * pz;
    k2vx     = h * vx_dot(pt, py, alpha, midx);

    midx  = x_0  + b31*k1x  + b32*k2x;
    midvx = vx_0 + b31*k1vx + b32*k2vx;

    k3t      = h * ft(pt, py, alpha, midx);
    k3x      = h * midvx;
    k3y      = h * fy(pt, py, alpha, midx);
    k3z      = h * pz;
    k3vx     = h * vx_dot(pt, py, alpha, midx);

    midx  = x_0  + b41*k1x  + b42*k2x  + b43*k3x;
    midvx = vx_0 + b41*k1vx + b42*k2vx + b43*k3vx;

    k4t      = h * ft(pt, py, alpha, midx);
    k4x      = h * midvx;
    k4y      = h * fy(pt, py, alpha, midx);
    k4z      = h * pz;
    k4vx     = h * vx_dot(pt, py, alpha, midx);

    midx  = x_0  + b51*k1x  + b52*k2x  + b53*k3x  + b54*k4x;
    midvx = vx_0 + b51*k1vx + b52*k2vx + b53*k3vx + b54*k4vx;

    k5t      = h * ft(pt, py, alpha, midx);
    k5x      = h * midvx;
    k5y      = h * fy(pt, py, alpha, midx);
    k5z      = h * pz;
    k5vx     = h * vx_dot(pt, py, alpha, midx);

    midx  = x_0  + b61*k1x  + b62*k2x  + b63*k3x  + b64*k4x  + b65*k5x;
    midvx = vx_0 + b61*k1vx + b62*k2vx + b63*k3vx + b64*k4vx + b65*k5vx;

    k6t      = h * ft(pt, py, alpha, midx);
    k6x      = h * midvx;
    k6y      = h * fy(pt, py, alpha, midx);
    k6z      = h * pz;
    k6vx     = h * vx_dot(pt, py, alpha, midx);

    newCoords.dt_= c1*k1t   + c2*k2t      + c3*k3t      + c4*k4t      + c5*k5t      + c6*k6t;
    newCoords.t_ = t_0      + newCoords.dt_;
    newCoords.x_ = x_0      + c1*k1x      + c2*k2x      + c3*k3x      + c4*k4x      + c5*k5x      + c6*k6x;
    newCoords.y_ = y_0      + c1*k1y	  + c2*k2y	    + c3*k3y	  + c4*k4y      + c5*k5y	  + c6*k6y;
    newCoords.z_ = z_0      + c1*k1z      + c2*k2z      + c3*k3z      + c4*k4z      + c5*k5z      + c6*k6z;
    newCoords.vx_ = vx_0	+ c1*k1vx	  + c2*k2vx		+ c3*k3vx	  + c4*k4vx		+ c5*k5vx	  + c6*k6vx;
    newCoords.Py_ = initCoords.Py_;
    newCoords.Pt_ = initCoords.Pt_;
    newCoords.Pz_ = initCoords.Pz_;
    newCoords.lambda_ = initCoords.lambda_;
    newCoords.alpha_ = initCoords.alpha_;

//    std::cout << "init t:" << initCoords.t_ << "new t:" << newCoords.t_ << std::endl;

    err[0] = (c1 - c1s)*k1t   + (c2-c2s)*k2t      + (c3-c3s)*k3t      + (c4-c4s)*k4t      + (c5-c5s)*k5t      + (c6-c6s)*k6t;
    err[1] = (c1 - c1s)*k1x   + (c2-c2s)*k2x      + (c3-c3s)*k3x      + (c4-c4s)*k4x      + (c5-c5s)*k5x      + (c6-c6s)*k6x;
    err[2] = (c1 - c1s)*k1y   + (c2-c2s)*k2y      + (c3-c3s)*k3y      + (c4-c4s)*k4y      + (c5-c5s)*k5y      + (c6-c6s)*k6y;
    err[3] = (c1 - c1s)*k1z   + (c2-c2s)*k2z      + (c3-c3s)*k3z      + (c4-c4s)*k4z      + (c5-c5s)*k5z      + (c6-c6s)*k6z;
    err[4] = (c1 - c1s)*k1vx  + (c2-c2s)*k2vx     + (c3-c3s)*k3vx     + (c4-c4s)*k4vx     + (c5-c5s)*k5vx     + (c6-c6s)*k6vx;
} // closes computeray

void ray::runray() {

    while(coords_.t_ < 10.0){
        step();
        stepcount_++;
        if(stepcount_ % 10 ==0) {
            cout << stepcount_ << "  " << coords_.t_ << "  " << coords_.x_ << "  " << coords_.y_ << endl;
        }
        if(stepcount_ > 10000){
            cout<<"took too many steps" <<endl;
            break;
        }
    }
}  // closes runray

double ray::fy(double pt, double py, double alpha,  double x) {
    return py + alpha * pt * x;
}

double ray::ft(double pt, double py, double alpha,  double x) {
    return pt*(1.0 - alpha*alpha * x*x ) - alpha * x * py;
}

double ray::vx_dot(double pt, double py, double alpha,  double x) {
    return -alpha * pt * (py + alpha * pt * x);
}