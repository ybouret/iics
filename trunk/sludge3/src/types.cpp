#include "types.hpp"


Hasher::  Hasher() throw() {}
Hasher:: ~Hasher() throw() {}

Hasher::KeyType Hasher:: getKey() throw()
{
    const Hasher::KeyType k = key<KeyType>();
    set();
    return k;
}

void Hasher:: operator()( const void *addr, size_t size) throw()
{
    assert( !(0==addr&&size>0) );
    run(addr,size);
}