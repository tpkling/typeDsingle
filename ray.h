//
// Created by Thomas Kling on 11/29/18.
//

#ifndef TYPEDSINGLE_RAY_H
#define TYPEDSINGLE_RAY_H

#include "Coords.h"


class ray {

public:

    ray(int i, int N);
    ~ray();
    void step();
    void runray();

    void computeRay(Coords initCoords, Coords &newCoords, double err[], double h);
    double fy(double pt, double py, double alpha,  double x);
    double ft(double pt, double py, double alpha,  double x);
    double vx_dot(double pt, double py, double alpha,  double x);

    double Py_lim1(double pt, double alpha, double x);
    double Py_lim2(double pt, double alpha, double x);
    double Pz_lim(double pt, double alpha, double x, double py);
    double PxSet(double pt, double alpha, double x, double py, double pz);

    Coords coords_;
    double h_;
    double totalError_ = 0.0;
    int stepcount_ = 0;
    int goodStepCount_ = 0;
    int badStepCount_ = 0;

};


#endif //TYPEDSINGLE_RAY_H
