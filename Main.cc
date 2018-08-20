/*
 * Main.cc
 *
 *  Created on: 2016-02-05
 *      Author: kassy
 */

#include <fcntl.h>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include "comms.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
	srand(genRandomSeed());

	int numCols = 0; // the number of symbol nodes to be created
	int numRows = 0; // the number of check nodes to be created
	int symbolDegreeDistribution = 0;

	bool nonGreedy = false;
	int maxDepth = 2;

    // Get the user input for the simulation
    getUserInput(numCols, numRows, symbolDegreeDistribution);

    ProgressiveEdgeGrowth PEG (numRows, numCols);

    int  checkNode = 0;
    int  depth = 0;
    int  girth = 0;
    int  oldCardinality = 0;
    int  newCardinality = 0;

    bool expandGraph = true;
    bool systematicParity = false;

    printf("Connecting symbol node        ");

    // Connect each symbol, one at a time
    while(!systematicParity)
    {
    	for(int symbolNode = 0; symbolNode < numCols; symbolNode++)
    	{
    		printf("\b\b\b\b\b\b\b%7d", symbolNode);
    		fflush(0);

    		// Connect each symbol edge, one at a time
    		for(int edge = 0; edge < symbolDegreeDistribution; edge++)
    		{
    			if(edge == 0) // If this is the first edge of the symbol
    			{
    				// Find the check node with the lowest degree.
    				PEG.SetMinimumIndices();
    				checkNode = PEG.GetIndexAtRandom(); // fetch a random check node from possibilities
    				PEG.ConnectEdge(checkNode, symbolNode); // connect an edge
    			}
    			else
    			{
    				PEG.ResetIndicators(); // reset the symbol and check node indicators

    				// Expand a subgraph from the current symbol node up to depth l
    				// Set up depth 0
    				depth = 0;
    				PEG.SetSymbolNodeIndicator(symbolNode);
    				PEG.SetCheckNodeIndicator();
    				oldCardinality = PEG.GetCardinality();
    				expandGraph = true;

    				while(expandGraph)
    				{
    					depth++;
    					PEG.SetSymbolNodeIndicator();
    					PEG.SetCheckNodeIndicator();
    					newCardinality = PEG.GetCardinality();

    					// if cardinality stops increasing, but is less than the number of check nodes (i.e. no cycle created)
    					if( (oldCardinality == newCardinality) && (newCardinality < numRows) )
    					{
    						depth--;
    						expandGraph = false;
    					}
    					else if( (oldCardinality < numRows) && (newCardinality == numRows) )
    					{
    						depth--;
    						PEG.SetCheckNodeIndicator(true);
    						expandGraph = false;
    					}else if(nonGreedy && (depth == maxDepth))
    						expandGraph = false;

    					oldCardinality = newCardinality;
    				} // while(expandGraph)

    				PEG.SetPossibleIndices();
    				checkNode = PEG.GetIndexAtRandom(); // fetch a random check node from possibilities
    				PEG.ConnectEdge(checkNode, symbolNode); // connect an edge
    			} // else
    		} // for(edge < symbolDegreeDistribution)
    	} // for(symbolNode < numCols)

    	// Check if this parity check matrix can be put into systematic form
    	systematicParity = PEG.CheckSystematic();

    }// while(!systematicParity)

    // Write the parity check matrix to a file that is readable by MATLAB
    PEG.PrintToFile();
    PEG.DisplayCheckDegreeDistribution();

    // Find the girth of the graph
    for(int symbolNode = 0; symbolNode < numCols; symbolNode++)
    {
    	depth = 0;
    	PEG.ResetIndicators();
		PEG.SetSymbolNodeIndicator(symbolNode);
		PEG.SetCheckNodeIndicator();
		oldCardinality = PEG.GetSymbolCardinality();
		expandGraph = true;

    	while(expandGraph)
    	{
    		depth++;
			PEG.SetSymbolNodeIndicator();
			PEG.SetCheckNodeIndicator();
			newCardinality = PEG.GetSymbolCardinality();

			if(oldCardinality == newCardinality)
			{
				girth = 0;
				PEG.SetSymbolGirth(symbolNode, girth);
				expandGraph = false;
			}
			if( (oldCardinality < numCols) && (newCardinality == numCols) )
			{
				depth--;
				girth = 2*(depth + 2);
				PEG.SetSymbolGirth(symbolNode, girth);
				expandGraph = false;
			}
    	}// while
    }// for

    PEG.DisplayGirthInformation();

	return 0;
}



