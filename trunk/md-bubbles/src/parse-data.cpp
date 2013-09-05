#include "yocto/exception.hpp"
#include "yocto/fs/local-fs.hpp"
#include <iostream>
#include "yocto/string/tokenizer.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"

using namespace yocto;

typedef vector<string> Words;

static inline bool ReadLine( string &line, ios::istream &fp )
{
    line.clear();
    return fp.read_line(line) >= 0;
}

static void Output( const Words &words )
{
    for(size_t i=1; i <= words.size(); ++i)
    {
        std::cout << words[i];
        if(i<words.size()) std::cout << " ";
    }
    std::cout << std::endl;
}

static inline bool is_sep( char C ) { return C == ' ' || C == '\t'; }

int main( int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        //======================================================================
        //
        // Open File
        //
        //======================================================================
        if(argc<=1)
            throw exception("usage: %s data.dat", prog);
        
        Words          words;
        string         line;
        ios::icstream  fp( argv[1] );
        
        while( ReadLine(line,fp) )
        {
            //std::cerr << line << std::endl;
            words.free();
            tokenizer::split(words, line, is_sep);
            if(words.size()>0 && words[1] == "Step" )
            {
                std::cout << "#"; Output(words);
                while( ReadLine(line,fp) )
                {
                    words.free();
                    tokenizer::split(words, line, is_sep);
                    if( words.size() > 0 && words[1] == "Loop" )
                        break;
                    Output(words);
                }
            }
            
        }
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception in " << prog << std::endl;
    }
    return 1;
}
