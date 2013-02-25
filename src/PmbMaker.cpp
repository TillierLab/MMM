#include <cstdlib>

#include "classCmdLineArgParser.h"
#include "Alignment.hpp"
#include "DistanceMatrix.hpp"

using namespace std;


static void print_usage(string const& pname)
{	cout	<< "Create a PMB distance matrix file from a FASTA alignment (or from an existing matrix - for debugging)" << endl;
	cout	<< "Usage:\n" << pname << "( -a inputAlignmentFilename | -d inputDistancesFilename ) -o outputFilename  [-p inputFilenamePrefix] [-q] [-g coefficient of variataion]" << endl;
	cout	<< "\t-q : use the quick PMB method" << endl;
	cout	<< "\t-g : enable correction for unequal rates of evolution. Not compatible with quick PMB mode. The coefficient of variance parameter must satisfy: 0<CV<1" << endl;
	exit(1);
	
}
//


int main(int argc, const char* argv[])
{	string inFilename, inFilenamePrefix = "", outFilename="outfile";
	double gapThresh = -1; //disable thresholds
	
	CmdLineArgParser options(argc, argv);
	//parse required parameters
	bool haveAlignmentFile = options.parse("-a", &inFilename);
	bool haveDistancesFile = options.parse("-d", &inFilename);
	
	if ((!haveAlignmentFile && !haveDistancesFile) || (haveAlignmentFile && haveDistancesFile))
		print_usage(options.programName());
	
	options.parse("-o", &outFilename);
	options.parse("-p", &inFilenamePrefix);
		
	
	
	//parse optional parameters
	bool quickPmb,gamma;
	double cv;
	if (haveAlignmentFile)
	{	quickPmb = options.parse("-q");
		gamma = options.parse("-g", &cv);
		
		if (gamma && quickPmb)
			print_usage(options.programName());
		
		if (gamma && (cv <= 0.0 || cv >= 1.0))
			print_usage(options.programName());
	}
	
	//verify that all options have been parsed
	if (!options.empty())
	{	cout << "Unrecognized option: " << endl;
		options.print();
		cout << endl;
		print_usage(options.programName());
	}
	
	
	DistanceMatrix distanceMatrix(inFilename);
	if(haveAlignmentFile)
	{	//obtain the alignment
		Alignment alignment(inFilename);
		alignment.read_alignment(inFilenamePrefix+inFilename, gapThresh);
		
		int sequenceCount = alignment.get_nOfSequences();
		int seqLength = alignment.get_nOfColumns();
		
		cout << "Input file: " << inFilenamePrefix <<inFilename << ", sequences: " << sequenceCount << ", length: " << seqLength << endl;
		//cout << alignment.get_header(0) << endl;
		//cout << alignment.get_sequence(0) << endl;
		
		if (gamma)
			distanceMatrix.enable_gamma(cv);
		
		distanceMatrix.create_from_alignment(alignment, quickPmb);
	} else if (haveDistancesFile) 
	{	distanceMatrix.read_from_file(inFilenamePrefix+inFilename);
	}

	
	distanceMatrix.print_to_file(outFilename);
	
	cout << "Finishing" << endl;
	return 0;
}