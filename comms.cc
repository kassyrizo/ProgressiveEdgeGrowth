#include <armadillo>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <set>

#include "comms.h"

long genRandomSeed()
{
  long randnum;
  int fd = open ("/dev/urandom", O_RDONLY);
  if (fd != -1) {
    read( fd, &randnum, 4 );
    close(fd);
  }
  else{
    printf("ERROR: Can't read /dev/urandom.\n");
    exit(1);
  }

  if(randnum>0) randnum*=-1;

  if(!randnum){
    printf("ERROR: Zero seed.\n");
    exit(1);
  }

  return randnum;
}

void getUserInput(int & numSymbolNodes, int & numCheckNodes, int & symbolDegreeDistribution) //, bool & lookaheadEnhanced, bool & nonGreedy, bool & random)
{
    cout << "Please enter the desired number of symbol nodes >> ";
    cin >> numSymbolNodes;
    cout << "Please enter the desired number of check nodes >> ";
    cin >> numCheckNodes;
    cout << "Please enter the desired symbol degree distribution >> ";
    cin >> symbolDegreeDistribution;

    return;
}

bool checkForYesOrNo()
{
    while(1)
    {
        string yesOrNo;
        cin >> yesOrNo;

        if((yesOrNo.compare("Y") == 0) || (yesOrNo.compare("y") == 0) )
            return true;
        else if( (yesOrNo.compare("N") == 0) || (yesOrNo.compare("n") == 0))
            return false;

        cout << "You have entered an invalid option. Please try again >> " << flush;
    }

    return true;
}

ProgressiveEdgeGrowth::ProgressiveEdgeGrowth(int _rows, int _cols)
{
	// Parity matrix dimensions
	rows = _rows;
	cols = _cols;

	// Create 2D dynamic array for rowSupport and colSupport
	// Dynamic memory allocation for row support and column support of H matrix
	rowSupport = new int*[rows](); // each pointer in rowSupport points to an array that holds the column index of the ones in the H matrix for that particular row
	colSupport = new int*[cols](); // each pointer in colSupport points to an array that holds the row index of the ones in the H matrix for that particular column

	H = new int*[rows]();
	H_systematic = new int*[rows]();

	rowOnes = new int[rows](); // hold the number of ones in each row
	colOnes = new int[cols](); // hold the number of ones in each column

	for(int i = 0; i < rows; i++)
	{
		rowSupport[i] = new int[rowOnes[i]]();
		H[i] = new int[cols]();
		H_systematic[i] = new int[cols]();
	}

	for(int i = 0; i < cols; i++)
		colSupport[i] = new int[colOnes[i]]();

	checkNodeIndicator = new int[rows]();
	previousCheckNodeIndicator = new int[rows]();
	symbolNodeIndicator = new int[cols]();
	symbolGirth = new int[cols]();

	possibleIndices = new int[rows];
}

ProgressiveEdgeGrowth::~ProgressiveEdgeGrowth()
{
	// Free up the memory from dynamic allocations
	for(int i = 0; i < rows; i++)
	{
		delete [] rowSupport[i];
		delete [] H[i];
		delete [] H_systematic[i];
	}
	delete [] rowSupport;
	delete [] H;
	delete [] H_systematic;

	for(int i = 0; i < cols; i++)
		delete [] colSupport[i];
	delete [] colSupport;

	delete [] rowOnes;
	delete [] colOnes;

	delete [] checkNodeIndicator;
	delete [] previousCheckNodeIndicator;
	delete [] symbolNodeIndicator;
	delete [] symbolGirth;

	delete [] possibleIndices;
}

void ProgressiveEdgeGrowth::ConnectEdge(int _row, int _col)
{
	rowOnes[_row]++;
	colOnes[_col]++;

	// Add the element to rowSupport
	int* tempRow = new int[ rowOnes[_row] ]();

	for(int i = 0; i < rowOnes[_row] - 1; i++)
		tempRow[i] = rowSupport[_row][i];

	tempRow[ rowOnes[_row] - 1 ] = _col;

	delete[] rowSupport[_row];

	rowSupport[_row] = tempRow;

	// Add the element to colSupport
	int* tempCol; tempCol = new int[ colOnes[_col] ]();

	for(int i = 0; i < colOnes[_col] - 1; i++)
		tempCol[i] = colSupport[_col][i];

	tempCol[ colOnes[_col] - 1 ] = _row;

	delete[] colSupport[_col];

	colSupport[_col] = tempCol;

	/*cout << "rowSupport: ";
	for(int i = 0; i < rows; i++)
	{
		cout << "row " << i;
		for(int j = 0; j < rowOnes[i]; j++)
			cout << " " << rowSupport[i][j];
		cout << "\n            ";
	}
	cout << endl;

	cout << "colSupport: ";
	for(int i = 0; i < cols; i++)
	{
		cout << "col " << i;
		for(int j = 0; j < colOnes[i]; j++)
			cout << " " << colSupport[i][j];
		cout << "\n            ";
	}
	cout << endl;*/


	return;
}

void ProgressiveEdgeGrowth::ResetIndicators()
{
	// Reset the indicators
	for(int i = 0; i < rows; i++)
	{
		checkNodeIndicator[i] = 0;
		previousCheckNodeIndicator[i] = 0;
	}
	for(int i = 0; i < cols; i++)
		symbolNodeIndicator[i] = 0;

	return;
}

bool ProgressiveEdgeGrowth::CheckSystematic()
{
	// Create H matrix
	createHMatrix();

	// Try to put matrix into [Ik P'] form
	rrefMod2();

	// Check if the first part of the matrix is in the proper form
	if( checkForSystematic() )
		return true;
	else
	{
		// Delete the old matrix, and make new for a new potentially systematic matrix
		for(int i = 0; i < rows; i++)
		{
				delete [] rowSupport[i];
				delete [] H[i];
				rowOnes[i] = 0;
		}

		for(int i = 0; i < rrefRows; i++)
			delete [] H_systematic[i];

		delete [] rowSupport;
		delete [] H;
		delete [] H_systematic;

		for(int i = 0; i < cols; i++)
		{
			delete [] colSupport[i];
			colOnes[i] = 0;
		}
		delete [] colSupport;

		// Create new, reset versions
		rrefRows = rows;
		rowSupport = new int*[rows]();
		H = new int*[rows]();
		H_systematic = new int*[rows]();
		colSupport = new int*[cols]();

		for(int i = 0; i < rows; i++)
		{
			rowSupport[i] = new int[rowOnes[i]]();
			H[i] = new int[cols]();
			H_systematic[i] = new int[cols]();
		}

		for(int i = 0; i < cols; i++)
			colSupport[i] = new int[rows]();

		ResetIndicators();


		return false;
	}
}

void ProgressiveEdgeGrowth::SetSymbolNodeIndicator(int symbolNode)
{
	symbolNodeIndicator[symbolNode] = 1;

	/*cout << "symbolNodeIndicator: ";
	for(int i = 0; i < cols; i++)
		cout << symbolNodeIndicator[i] << " ";
	cout << endl;*/

	return;
}

void ProgressiveEdgeGrowth::SetSymbolNodeIndicator()
{
	// Iterate through the checkNodeIndicator set to find all
	// check nodes under the current graph condition

	for(int i = 0; i < rows; i++)
	{
		if(checkNodeIndicator[i] == 1)
		{
			// Find all attached symbol nodes, and indicate that they have been found at this depth
			for(int j = 0; j < rowOnes[i]; j++)
				symbolNodeIndicator[ rowSupport[i][j] ] = 1;
		}
	}

	/*cout << "symbolNodeIndicator: ";
	for(int i = 0; i < cols; i++)
		cout << symbolNodeIndicator[i] << " ";
	cout << endl;

	return;*/
}

void ProgressiveEdgeGrowth::SetCheckNodeIndicator()
{
	// First, set the previous indicator so it can be rolled back if need be
	for(int i = 0; i < rows; i++)
		previousCheckNodeIndicator[i] = checkNodeIndicator[i];

	// Iterate through the symbolNodeIndicator set to find all
	// symbol nodes under the current graph condition
	for(int i = 0; i < cols; i++)
	{
		if(symbolNodeIndicator[i] == 1)
		{
			// Find all attached check nodes, and indicate that they have been found at this depth
			for(int j = 0; j < colOnes[i]; j++)
				checkNodeIndicator[ colSupport[i][j] ] = 1;
		}
	}

	/*cout << "checkNodeIndicator: ";
	for(int i = 0; i < rows; i++)
		cout << checkNodeIndicator[i] << " ";
	cout << endl;*/

	return;
}

void ProgressiveEdgeGrowth::SetCheckNodeIndicator(bool previous)
{
	if(previous)
	{
		for(int i = 0; i < rows; i++)
			checkNodeIndicator[i] = previousCheckNodeIndicator[i];
	}

	/*(cout << "checkNodeIndicator: ";
	for(int i = 0; i < rows; i++)
		cout << checkNodeIndicator[i] << " ";
	cout << endl;*/

	return;
}

void ProgressiveEdgeGrowth::SetPossibleIndices()
{
	// Find the possible indices based on the check node indicator first
	numPossibleIndices = 0;

	for(int i = 0; i < rows; i++)
	{
		if(checkNodeIndicator[i] == 0)
		{
			numPossibleIndices++;
			int* temp = new int[numPossibleIndices]();

			for(int j = 0; j < numPossibleIndices - 1; j++)
				temp[j] = possibleIndices[j];
			temp[numPossibleIndices - 1] = i;

			delete[] possibleIndices;
			possibleIndices = temp;
		}
	}

	// Find the minimum edges in the degree distribution
	int minimum = getMinimum(true);

	// Erase indices from possibleIndices that do not have the minimum
	for(int i = 0; i < numPossibleIndices; i++)
	{
		if( rowOnes[ possibleIndices[i] ] != minimum)
		{
			numPossibleIndices--;
			int* temp = new int[numPossibleIndices];

			for(int j = 0; j < numPossibleIndices; j++)
				temp[j] = j < i ? possibleIndices[j] : possibleIndices[j+1];

			delete[] possibleIndices;
			possibleIndices = temp;
			i--;
		}
	}

	/*"\npossibleIndices: ";
	for(int i = 0; i < numPossibleIndices; i++)
		cout << possibleIndices[i] << " ";
	cout << endl;*/

	return;
}

void ProgressiveEdgeGrowth::SetMinimumIndices()
{
	int minimum = getMinimum(false);
	numPossibleIndices = 0;

	// Count the number of check nodes that have the minimum value
	for(int i = 0; i < rows; i++)
	{
		if(rowOnes[i] == minimum)
		{
			numPossibleIndices++;

			int* temp = new int[numPossibleIndices]();

			for(int j = 0; j < numPossibleIndices - 1; j++)
				temp[j] = possibleIndices[j];

			temp[numPossibleIndices - 1] = i;

			delete[] possibleIndices;

			possibleIndices = temp;
		}
	}

	/*"\npossibleIndices: ";
	for(int i = 0; i < numPossibleIndices; i++)
		cout << possibleIndices[i] << " ";
	cout << endl;*/

	return;
}

void ProgressiveEdgeGrowth::PrintToFile()
{
	// Create the file
	ofstream H_file;
	stringstream stream;

	stream << "H_" << cols << "_" << rows << ".m";
	const char* fileName = stream.str().c_str();
	H_file.open(fileName);
	H_file << "H = [ ";

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
			H_file << H[i][j] << " ";

		if(i == (rows - 1))
			H_file << "];";
		else
			H_file << ";\n      ";
	}
	H_file.close();

	stream.str("");
	stream.clear();

	stream << "PEG" << rows << "x" << cols << ".txt";
	fileName = stream.str().c_str();
	H_file.open(fileName);

	for(int i = 0; i < rows; i++)
		H_file << rowOnes[i] << " ";

	H_file << "\n";

	for(int i = 0; i < cols; i++)
	{
		for(int j = 0; j < colOnes[i]; j++)
		{
			H_file << left << setw(4) << colSupport[i][j] << " ";
		}
		H_file << "\n";
	}

	H_file.close();

	return;
}

void ProgressiveEdgeGrowth::DisplayCheckDegreeDistribution()
{
	std::set<int> checkDegree(rowOnes, rowOnes + rows);

	int number = 0;

	cout << "\n";
	for( set<int>::iterator i = checkDegree.begin(); i != checkDegree.end(); ++i)
	{
		number = 0;
		for(int j = 0; j < rows; j++)
	    	if(rowOnes[j] == *i)
	    		number++;

		cout << "There are " << number << " check nodes with " << *i << " edges." << endl;
	}
	return;
}

void ProgressiveEdgeGrowth::DisplayGirthInformation()
{
	std::set<int> symGirth(symbolGirth, symbolGirth + cols);

	int girth = 0;

	cout << "\n\n";
	for( set<int>::iterator i = symGirth.begin(); i != symGirth.end(); i++)
	{
		girth = 0;
		for(int j = 0; j < cols; j++)
			if(symbolGirth[j] == *i)
				girth++;

		cout << "There are " << girth << " symbol nodes with local girth " << *i << "." << endl;
	}

	return;
}

void ProgressiveEdgeGrowth::SetSymbolGirth(int symbol, int girth)
{
	symbolGirth[symbol] = girth;
}

int ProgressiveEdgeGrowth::GetIndexAtRandom()
{
	int index = 0;

	index = possibleIndices[ std::rand() % numPossibleIndices ];

	return index;
}

int ProgressiveEdgeGrowth::GetCardinality()
{
	cardinality = 0;

	for(int i = 0; i < rows; i++)
		if(checkNodeIndicator[i] == 1)
			cardinality++;

	return cardinality;
}

int ProgressiveEdgeGrowth::GetSymbolCardinality()
{
	cardinality = 0;

	for(int i = 0; i < cols; i++)
		if(symbolNodeIndicator[i] == 1)
			cardinality++;

	return cardinality;
}

int ProgressiveEdgeGrowth::getMinimum(bool rowSubSet)
{
	int minimum = 0;

	if(rowSubSet) // If the minimum of the possible indices is required
	{
		for(int i = 0; i < numPossibleIndices; i++)
		{
			if(rowOnes[ possibleIndices[i] ] < minimum || i == 0)
				minimum = rowOnes[ possibleIndices[i] ];
		}

	}
	else
	{
		for(int i = 0; i < rows; i++) // If the minimum of the entire check degree distribution is required
		{
			if(rowOnes[i] < minimum || i == 0)
				minimum = rowOnes[i];
		}
	}

	return minimum;
}

void ProgressiveEdgeGrowth::createHMatrix()
{
	// Clear the H matrices
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
			{
				H[i][j] = 0;
				H_systematic[i][j] = 0;
			}
	}

	rrefRows = rows;

	// Create the parity check matrix
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < rowOnes[i]; j++)
		{
			H[i][ rowSupport[i][j] ] = 1;
			H_systematic[i][ rowSupport[i][j] ] = 1;
		}
	}

	return;
}

void ProgressiveEdgeGrowth::rrefMod2()
{
	// This doesn't actually put matrix into reduced row echelon form if [Ik P] form can't be found
	// But that's fine for what this program needs.

	/*for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
			cout << H_systematic[i][j] << " ";
		cout << endl;
	}
	cout << endl;*/

	// Swap rows so that a diagonal of ones is formed
	for(int i = 0; i < rows; i++)
	{
		// Check if the (i,i) element is a 1, if it's not, switch with the next row that has a leading one in the (i,i)th position
		if(H_systematic[i][i] == 0)
		{
			// Check for next row with a one in the proper position
			for(int j = i+1; j < rows; j++)
			{
				if( H_systematic[j][i] == 1)
				{
					swapRows(i,j);
					break;
				}
			}

		}

		/*for(int k = 0; k < rows; k++)
		{
			for(int j = 0; j < cols; j++)
				cout << H_systematic[k][j] << " ";
			cout << endl;
		}
		cout << endl;*/

		// Clear any ones before or after leading one
		for(int j = 0; j < rows; j++)
		{
			if( H_systematic[j][i] == 1 && (j != i))
				mod2AddRows(j,i);
		}

		/*for(int k = 0; k < rows; k++)
		{
			for(int j = 0; j < cols; j++)
				cout << H_systematic[k][j] << " ";
			cout << endl;
		}
		cout << endl;*/
	}

	// Remove any rows of all zero
	for(int i = 0; i < rrefRows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			if(H_systematic[i][j] != 0) // if there are any ones, don't erase the row
				break;
			else if( j == (cols - 1) ) // if there are no ones, erase the row
			{
				// erase the row
				delete [] H_systematic[i];
				rrefRows--;

				// Create temporary matrix to copy new matrix to
				int ** tempH = new int*[rrefRows]();
				for(int k = 0; k < rrefRows; k++)
				{
					tempH[k] = new int[cols]();
					for(int p = 0; p < cols; p++)
					{
						tempH[k][p] = H_systematic[k][p];
					}
					delete [] H_systematic[k]; // delete the copied row
				}

				// delete the handle to the old H_systematic H matrix
				delete [] H_systematic;

				H_systematic = tempH;
				i--; // to make sure we start at the next row and don't skip one.
			}
		}
	}

	/*for(int k = 0; k < rrefRows; k++)
	{
		for(int j = 0; j < cols; j++)
			cout << H_systematic[k][j] << " ";
		cout << endl;
	}
	cout << endl;*/


	return;
}

void ProgressiveEdgeGrowth::swapRows(int row1, int row2)
{
	int temp;

	for(int i = 0; i < cols; i++)
	{
		temp = H_systematic[row1][i];
		H_systematic[row1][i] = H_systematic[row2][i];
		H_systematic[row2][i] = temp;
	}

	return;
}

void ProgressiveEdgeGrowth::mod2AddRows(int row1, int row2)
{
	// This will do Row1 = Row1 + Row2

	for(int i = 0; i < cols; i++)
		H_systematic[row1][i] = H_systematic[row1][i] ^ H_systematic[row2][i];

	return;
}

bool ProgressiveEdgeGrowth::checkForSystematic()
{
	bool systematic = true;

	// Check if each [i][i] element is 1, it should be guaranteed to be 0's in every other column
	for(int i = 0; i < rrefRows; i++)
	{
		if(H_systematic[i][i] != 1)
		{
			systematic = false;
			break;
		}
	}

	return systematic;
}
