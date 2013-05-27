#include "tracer.hpp"

Tracer:: Tracer() throw() :
prev(0),next(0),
pos()
{}

Tracer:: ~Tracer() throw() {}

Tracer:: Ring:: Ring() throw() {}

Tracer:: Ring:: ~Ring() throw() { auto_delete(); }
