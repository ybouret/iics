#include "marker.hpp"

Marker::~Marker() throw() {}

Marker:: Marker() throw() : coord(), next(0), prev(0) {}

void Marker:: reset() throw()
{
    coord.ldz();
}