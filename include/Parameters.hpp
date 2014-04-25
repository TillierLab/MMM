#ifndef _PARAMETERS_H
#define _PARAMETERS_H


#include <cstdlib>
#include <limits>
#include <iostream>


using namespace std;

struct Parameters
{


	string distfileName1, distfileName2, batchFileName, outfileName, distfileNamesPrefix, stateFileName;
	string reqTaxName;
	int strictSize;
	int minSize;
	int maxTrees;
	int rangeStart;
	int rangeEnd;
	bool tabulateTop;
	bool useTaxInfo;
	bool reqTax;
	bool tabDelimitedTable;
	bool brief;
	bool verbose;
	bool doShowStats;
	bool outputZeros;
	bool stateRestore;
	bool exactTime;
	bool useHardware;
	double allow;
	
	// MMML extensions: New 2010.02.28: RLC:
	bool MMMLmode;
	string treeFileName, ignoreTaxaFileName;
	int minBootstrapSupport;
	bool listTaxa, topCoverage, LverboseMode; 
	// New 2010.02.28: RLC.

	
	double cutoffLower, cutoffUpper;
	bool logDistances;
	
	bool noHeader,onlyHeader, onlyMP, quickMP;
	bool ignoreErrors;
	
	bool aln2pmb, printPmb, printSlice;
	double gapThreshold; //by default gap threshold for alignments is disabled
	unsigned int alnCacheSizeLimit;
	unsigned int distanceMarixCacheSizeLimit;
	unsigned int memoryUsageLimitMegabytes;
	
	
	
	
	bool gzipIn;
	bool gzipOut;


	// Constructor
	Parameters();
	~Parameters(){};
	//bool isParameter( char *arg );
	//bool checkNonBoolParameter( int argNumber, int argc, char **argv, string myArg );
	//bool checkBoolParameter( int argNumber, int argc, char **argv, string myArg );
	//void readCommandLineOptions( int argc, char **argv );
	//void printUsage();
	
	
	
	
	static void printUsage1(const string& pname);
	static void printUsage2(const string& pname);
	static void printUsage(const string& pname);
	static void printVerboseUsage(const string& pname);
	void process(int argc, const char* argv[]);
	
	
};

#endif