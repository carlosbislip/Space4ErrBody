#ifndef GET_STUFF_H
#define GET_STUFF_H

#include <iostream>
#include <ctime>
#include <sstream>

std::vector< std::string > getFiles ( const std::string &filename );

std::vector< double > getData ( const std::string &filename );

std::string getCurrentDateTime( const bool useLocalTime = false );

#endif // GET_STUFF_H
