#include "../segmenter.hpp"



void Segmenter:: find_bracketing_junctions(const Spot *spot) const
{
    assert(spot);
    
    //--------------------------------------------------------------------------
    // Use spot info on the grid
    //--------------------------------------------------------------------------
    const Coord   klo    = spot->klo;
    const Coord   kup    = spot->kup;
    const Tracer *tracer = spot->handle; assert(tracer);
    const Vertex  v      = tracer->vertex;
   
#if 0
    const Junctions &HSegLo = Horz(klo.y);
    const Junctions &HSegUp = Horz(kup.y);
    const Junctions &VSegLo = Vert(klo.x);
    const Junctions &VSegUp = Vert(kup.x);
#endif
    
    
}