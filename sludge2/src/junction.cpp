#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction() throw() :
next(0),
prev(0),
pos(),
klo(0),
khi(0),
bubble(0),
alpha(0),
curvature(0),
n()
{}

