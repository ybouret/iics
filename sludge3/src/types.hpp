#ifndef SLUDGE_TYPE_INCLUDED
#define SLUDGE_TYPE_INCLUDED 1

#include "yocto/math/v2d.hpp"
#include "yocto/hashing/sha1.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace math;

typedef double        Real;
typedef v2d<Real>     Vertex;
typedef hashing::sha1 HasherType; //!< base class
#define               REAL_TYPE MPI_DOUBLE

//! mostly used to debug MPI
class Hasher : public HasherType
{
public:
    typedef unsigned KeyType;
    
    explicit Hasher() throw();
    virtual ~Hasher() throw();
    void operator()( const void *addr, size_t size) throw();
    
    KeyType getKey() throw(); //!< get key and reset
    template <typename T>
    inline void operator()( const T &obj ) throw()
    {
        run(&obj,sizeof(obj));
    }
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Hasher);
};


#endif

