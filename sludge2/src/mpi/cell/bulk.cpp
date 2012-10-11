#include "../cell.hpp"

void Cell:: build_bulk()
{
    for(unit_t j=Y.lower;j<Y.upper;++j)
    {
        const Array1D &B_j    = B[j];
        Array1D       &Bulk_j = Bulk[j];
        const Array1D &B_jp   = B[j+1];
        
        for(unit_t i=X.lower;i<X.upper;++i)
        {
            Real &bulk = Bulk_j[i];
            bulk = 0;
            const unit_t ip = i+1;
            if( B_j[i]   <= 0) ++bulk;
            if( B_j[ip]  <= 0) ++bulk;
            if( B_jp[i]  <= 0) ++bulk;
            if( B_jp[ip] <= 0) ++bulk;
        }
    }
}