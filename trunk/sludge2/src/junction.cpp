#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction() throw() : next(0), prev(0), pos(), bubble(0), curvature(0) {}

