#pragma once

#include <vector> 

class Point {
    private:
    std::vector<double> xyz_;
    double x_;
    double y_;
    double z_; 

    public:
    Point(double n1, double n2, double n3);
    Point(std::vector<double>); 
    double get_x() { return x_; }
    double get_y() { return y_; }
    double get_z() { return z_; }

}; 