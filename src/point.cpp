#include "point.h"

Point::Point(double x, double y, double z):
    xyz_({x,y,z}), x_(x), y_(y), z_(z)  {}
Point::Point(std::vector<double> vec):
    xyz_(vec), x_(vec[0]), y_(vec[1]), z_(vec[2]) {}