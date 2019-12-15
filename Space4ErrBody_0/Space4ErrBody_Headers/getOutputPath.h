#ifndef BISLIP_OUTPUTPATH_H
#define BISLIP_OUTPUTPATH_H

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Basics/utilities.h"

namespace bislip { namespace Variables {

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "getOutputPath.h" ).length( ) );
    std::string outputPath = reducedPath + "executables/SimulationOutput/";

    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

}

}
#endif // BISLIP_OUTPUTPATH_H
