#ifndef GET_STUFF_H
#define GET_STUFF_H

#include <iostream>
#include <ctime>
#include <sstream>

std::vector< std::string > getFiles ( const std::string &filename );

std::vector< double > getData ( const std::string &filename );

std::string getCurrentDateTime( const bool useLocalTime = false );

double getAngularDistance ( const double &lat_i,
    const double &lon_i,
    const double &lat_f,
    const double &lon_f);

double getHeadingToTarget ( const double &lat_i,
    const double &lon_i,
    const double &lat_f,
    const double &lon_f);


#endif // GET_STUFF_H
