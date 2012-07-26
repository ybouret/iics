#include "segmenter.hpp"

static inline void __display_segments( const Segment::List &a, Real ya, const Segment::List &b , Real yb)
{
    const Segment *sa = a.head;
    const Segment *sb = b.head;
    
    fprintf(stderr,"A@%8.4f:",ya);
    while(sa)
    {
        Junction *J = sa->handle;
        fprintf(stderr," (%8.4f,%8.4f)", J->vertex.x, J->vertex.y );
        sa=sa->next;
    }
    fprintf(stderr,"\n");
    
    fprintf(stderr,"B@%8.4f:",yb);
    while(sb)
    {
        Junction *J = sb->handle;
        fprintf(stderr," (%8.4f,%8.4f)", J->vertex.x, J->vertex.y );
        sb=sb->next;
    }
    fprintf(stderr,"\n");
    
}

void Segmenter:: merge_pbc( Segment::List &a, Real ya, Segment::List &b , Real yb)
{
    static const Real xtol = 10 *  numeric<Real>::ftol;
    
    Segment::List new_a(s_cache);
    Segment::List new_b(s_cache);
    
    if(false)
    {
        fprintf(stderr,"Old Segments\n");
        __display_segments(a, ya, b, yb);
    }
    //--------------------------------------------------------------------------
    //
    // common part
    //
    //--------------------------------------------------------------------------
    while( a.size>0 && b.size>0 )
    {
        Segment   *sa = a.pop_front();
        Segment   *sb = b.pop_front();
        Junction  *Ja = sa->handle;
        Junction  *Jb = sb->handle;
        const Real dx = Ja->vertex.x - Jb->vertex.x;
        fprintf( stderr, "dx=%g\n", dx);
        if( Fabs(dx)<=xtol )
        {
            //------------------------------------------------------------------
            // these are the same: move directly into new
            //------------------------------------------------------------------
            new_a.push_back(sa);
            new_b.push_back(sb);
        }
        else
        {
            //------------------------------------------------------------------
            // they are different
            //------------------------------------------------------------------
            if( dx < 0 )
            {
                assert(Ja->vertex.x < Jb->vertex.x );
                // move sa, restore sb, and move a copy of sa into sb
                new_a.push_back(sa);
                b.push_front(sb);
                Jb = junctions.append();
                Jb->copy(Ja);
                Jb->vertex.y = yb;
                new_b.attach(Jb);
            }
            else
            {
                //move sb, restore sa, and move a copy of sb into sa
                new_b.push_back(sb);
                a.push_front(sa);
                Ja = junctions.append();
                Ja->copy(Jb);
                Ja->vertex.y = ya;
                new_a.attach(Ja);
            }
        }
    }
    
    //--------------------------------------------------------------------------
    //
    // remaining a
    //
    //--------------------------------------------------------------------------
    while( a.size > 0 )
    {
        Segment   *sa = a.pop_front();
        new_a.push_back(sa);
        Junction *Jb = junctions.append();
        Jb->copy( sa->handle );
        Jb->vertex.y = yb;
        new_b.attach(Jb);
    }
    
    //--------------------------------------------------------------------------
    //
    // remaining b
    //
    //--------------------------------------------------------------------------
    while( b.size > 0 )
    {
        Segment   *sb = b.pop_front();
        new_b.push_back(sb);
        Junction *Ja = junctions.append();
        Ja->copy( sb->handle );
        Ja->vertex.y = ya;
        new_a.attach(Ja);
    }
    
    if(false)
    {
        fprintf(stderr,"New Segments\n");
        __display_segments(new_a, ya, new_b, yb);
    }
    
    new_a.swap_with(a);
    new_b.swap_with(b);
    
}
