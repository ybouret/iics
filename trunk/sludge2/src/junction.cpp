#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction() throw() :
next(0),
prev(0),
vertex(),
klo(0),
khi(0),
bubble(0),
alpha(0),
curvature(0),
pressure(0),
t(),
n(),
gradP_t(0)
{}

