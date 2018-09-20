#include <iostream>
#include <ctime>
#include <sstream>


std::vector< std::string > getFiles ( const std::string &filename );

std::vector< double > getData ( const std::string &filename );

std::string getCurrentDateTime( const bool useLocalTime = false );

std::vector< double > getFitness (const std::vector< double > &input_data_, const std::vector< double > &output_settings_, const std::string &outputSubFolder_, const std::vector< double > &x);


