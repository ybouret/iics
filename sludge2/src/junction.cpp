#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction( const Vertex &at ) throw() : next(0), prev(0), pos(at) {}

