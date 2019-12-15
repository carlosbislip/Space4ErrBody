#include "Space4ErrBody_Headers/Space4ErrBodyProblem.h"
#include "Space4ErrBody_Headers/runInputArgumentEvaluator.h"
#include "Space4ErrBody_Headers/getPagmoAlgorithm.h"
#include <Tudat/Bislip/bislipUtilities.h>
#include <Tudat/Bislip/bislipVariables.h>
#include <Tudat/Bislip/bislipConstraints.h>
#include <Tudat/Bislip/bislipVehicleParameters.h>

#include <Tudat/Bislip/bislipDecisionVectorEvaluation.h>
#include <Tudat/Bislip/createTUDATSettings.h>
#include <Tudat/Bislip/updateGuidance.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    std::string primaryArgument = argv[ 1 ];

    if( primaryArgument == "batch"  || primaryArgument == "Batch" || primaryArgument == "BATCH" )
    {
        const std::vector< std::vector< std::string > > runList = bislip::utilities::readFileToVectorOfStringVectors( argv[ 2 ] , "" );

        std::cout << runList.size() << " Runs to Evaluate" << std::endl;
        int p = 0;
        for( unsigned int i = 0; i < runList.size(); i++ )
        {
              try
              {
                std::cout << "  Starting run " << i + 1 << " / " << runList.size() << std::endl;
                std::cout << " " << std::endl;

                bislip::runInputArgumentEvaluator( runList.at( i ) );

                std::cout << " " << std::endl;
                std::cout << "  Completed run " << i + 1 << " / " << runList.size() << std::endl;
                std::cout << " " << std::endl;
                p+= 1;
              }
              catch (...) // any exception happens
              {
                std::cout << "Something bad happened" << std::endl;
                throw;
              }
        }
        std::cout<< p << " Runs Evaluated" << std::endl;
    }
    else
    {
        //try
       // {
            std::cout<< "Single Run to Evaluate" << std::endl;
            std::cout << " " << std::endl;
            std::vector< std::string > runInputArguments( argv, argv + argc );
            runInputArguments.erase( runInputArguments.begin() );

            bislip::runInputArgumentEvaluator( runInputArguments );

            std::cout << " " << std::endl;
            std::cout<< "Single Run Evaluated" << std::endl;
            std::cout << " " << std::endl;
       // }
       // catch (...) // any exception happens
       // {
        //   throw( "Some error occured" );
       // }
    }
    std::cout<< "Closing Executable" << std::endl;
    return 0;
}
