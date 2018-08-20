#ifndef COMMS_H
#define COMMS_H

using namespace arma;
using namespace std;


long genRandomSeed();

void getUserInput(int & numSymbolNodes, int & numCheckNodes, int & symbolDegreeDistribution);//, bool & lookaheadEnhanced, bool & nonGreedy, bool & random);

bool checkForYesOrNo();

class ProgressiveEdgeGrowth
{
	public:

		ProgressiveEdgeGrowth(int _rows, int _cols);
		~ProgressiveEdgeGrowth();

		void 	ConnectEdge(int _row, int _col);
		void 	ResetIndicators();
		bool 	CheckSystematic();

		void 	SetSymbolNodeIndicator(int symbolNode);
		void 	SetSymbolNodeIndicator();
		void 	SetCheckNodeIndicator();
		void 	SetCheckNodeIndicator(bool previous);
		void 	SetMinimumIndices();
		void 	SetPossibleIndices();
		void 	SetSymbolGirth(int symbol, int girth);

		int		GetIndexAtRandom();
		int		GetCardinality();
		int		GetSymbolCardinality();
		void 	PrintToFile();
		void 	DisplayCheckDegreeDistribution();
		void 	DisplayGirthInformation();

	private:
		int rows;
		int cols;

		int rrefRows;

		int** rowSupport;
		int** colSupport;
		int*  rowOnes; // equivalent to the check degree distribution
		int*  colOnes;

		int* checkNodeIndicator;
		int* previousCheckNodeIndicator;
		int* symbolNodeIndicator;

		int* possibleIndices;

		int* symbolGirth;

		int numPossibleIndices;
		int cardinality;

		int** H;
		int** H_systematic;

		int		getMinimum(bool rowSubSet);
		void 	createHMatrix();
		void 	rrefMod2();
		void	swapRows(int row1, int row2);
		void	mod2AddRows(int Row1, int Row2);
		bool 	checkForSystematic();
};

#endif
