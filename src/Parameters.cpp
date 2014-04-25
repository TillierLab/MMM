#include "Parameters.hpp"
#include "classCmdLineArgParser.h"

Parameters::Parameters()
{
	

	strictSize = -1;
	minSize = 0;
	 maxTrees = numeric_limits<int>::max();
	rangeStart = -1, rangeEnd = numeric_limits<int>::max();
	
	tabulateTop = false;
	useTaxInfo = false;
	reqTax = false;
	tabDelimitedTable = false;
	brief = false;
	verbose = false;
	doShowStats = false;
	outputZeros = false;
	stateRestore = false;
	exactTime = false;
	useHardware = false;
	
	// MMML extensions: New 2010.02.28: RLC:
	minBootstrapSupport = 1;
	MMMLmode = false, listTaxa = false, topCoverage = false, LverboseMode = false; 
	// New 2010.02.28: RLC.

	
	cutoffLower = 0.000011; //protdist's threshold of 0.000010 fits here
	cutoffUpper = 1.0e+10;
	logDistances = false;
	
	noHeader = false;
	onlyHeader = false;
	onlyMP = false;
	quickMP = false;
	ignoreErrors = false;
	
	aln2pmb = false;
	printPmb = false;
	printSlice = false;
	gapThreshold = -1; //by default gap threshold for alignments is disabled
	alnCacheSizeLimit = 10;
	distanceMarixCacheSizeLimit = 1000;
	memoryUsageLimitMegabytes = 1024;
	
	
	gzipIn = false;
	gzipOut = false;
}	


	
void Parameters::printUsage1(const string& pname)
{
	cout << "Usage 1: " << pname <<
		" -d1 distfile1 -d2 distfile2 -o outfile -a allowance [-brief] " <<
		"[-table [-top]] [-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] " <<
		"[-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb [-gapthresh alignmentGapThreshold] [-printPmb] [-printSlice]] " <<
		"[-v] [-restore stateFileName]" << endl;
}
//

void Parameters::printUsage2(const string& pname)
{
	cout << "Usage 2: " << pname <<
		" -b batchfile -o outfile -a allowance [-brief] " <<
		"[-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] " <<
		"[-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb [-gapthresh alignmentGapThreshold] [-printPmb] [-printSlice]] " <<
		"[-v]" << endl;
}
//

void Parameters::printUsage(const string& pname)
{
	printUsage1(pname);
	cout << "...or..." << endl;
	printUsage2(pname);
	cout << "...or..." << endl;
	cout << "Usage: " << pname << " -help" << endl;
	cout << "MMML extensions, incompatible with -brief or -z: -Lt rootedTreeFileName [-LtopCov] [-Li ignoreTaxaFileName] " <<
		"[-Lb minBootstrapSupport] [-List] [-Lv]" << endl;
	exit(1);
}
//

void Parameters::printVerboseUsage(const string& pname)
{
	cout << endl;
	printUsage1(pname);
	cout << "...or..." << endl;
	printUsage2(pname);
	cout << "\nRequired parameters:\n" << endl;
	cout << "-d1 distfile1 (if -b batchfile is not specified)" << endl;
	cout << "\tSpecifies the file name of the first distance matrix." << endl;
	cout << "-d2 distfile2 (if -b batchfile is not specified)" << endl;
	cout << "\tSpecifies the file name of the second distance matrix." << endl;
	cout << "-b batchfile (if -d1 distfile1 -d2 distfile2 are not specified)" << endl;
	cout << "\tSpecifies the name of a file containing tab-delimited pairs of distance matrix file names, one pair per line." << endl;
	cout << "-o outfile" << endl;
	cout << "\tThe desired name of the output file." << endl;
	cout << "-a allowance" << endl;
	cout << "\tSpecifies a value between 0.0 and 1.0, with smaller values requiring more precise matching, and "
		  << "larger values tolerating looser matching." << endl;
	cout << "\nOptional parameters:\n" << endl;
	cout << "-brief" << endl;
	cout << "\tOutputs a single tab-delimited line of output for the pair of distance matrices." << endl;
	cout << "-table (valid only when -b batchfile is not specified)" << endl;
	cout << "\tProduces a tab-delimited matrix that lists the "
		  << "fit-weighted proportion of matches for each entry from the two distance matrices, "
		  << "over all common submatrices. For example, if seq1 matched seq2 in two of ten discovered submatrices, with "
		  << "a score of 0.8 and 0.9, respectively, then the entry for seq1_seq2 would be (0.8 + 0.9) / 10 = 0.17" << endl;
	cout << "-top (valid only when -table is specified)" << endl;
	cout << "\tRestricts the -table output to just the best-scoring submatrix." << endl;
	cout << "-u" << endl;
	cout << "\tOnly attempts to match entries from the pair of distance matrices that have the same tag in their name. A tag "
		  << "ends with a pipe (\"|\") character, e.g. Ecoli in Ecoli|DnaA, or Hsapiens in Hsapiens|myc" << endl;
	cout << "-req reqTaxName" << endl;
	cout << "\tRequires that a submatrix include at least one entry with the specified reqTaxName, e.g. Hsapiens." << endl;
	cout << "-sz strictSize" << endl;
	cout << "\tReturns common submatrices of size strictSize, if any. Does not attempt to find larger submatrices." << endl;

	cout << "-p distfileNamesPrefix" << endl;
	cout << "\tThe prefix is appended to the path of all distance matrices." << endl;
	cout << "-z" << endl;
	cout << "\tForces the output of zeros for combinations that have not reached the minSize (to be used in 'brief' mode)." << endl;
	cout << "-minsize minSize" << endl;
	cout << "\tReturns common submatrices of at least minSize, if any. Will not return smaller submatricies. Use this option to save memory with large matrices." << endl;
	cout << "-maxtrees maxTrees" << endl;
	cout << "\tLimits the number of returned trees. Saves memory." << endl;
	cout << "-rangeStart start -rangeEnd end" << endl;
	cout << "\tFor distributed computing, launch this program multiple times with non-overlapping ranges, [start, end), "
		  << "start >= 0 and end <= distfile1's matrix size (auto-adjusted when -b is specified). Incompatible with -rndOrder." << endl;
	cout << "-quickmp" << endl;
	cout << "\tIf specified, the quickly but naively estimated match potential is included in batch output." << endl;
	cout << "-v" << endl;
	cout << "\tOutputs the size of the largest common submatrix found so far, to stdout." << endl;
	cout << "-onlymp" << endl;
	cout << "\tDoes not perform MMM analysis, only prints MatchPotential (brief and MMML modes only)." << endl;
	cout << "-restore stateFileName" << endl;
	cout << "\tResumes work from a previously interrupted MMM session" << endl;
	cout << "-exacttime" << endl;
	cout << "\tDo not round running time to nearest second, in brief mode" << endl;
	cout << "-upper cutoffUpper" << endl;
	cout << "\tAll distances above cutoffUpper are set to zero." << endl;
	cout << "-lower cutoffLower" << endl;
	cout << "\tAll distances below cutoffLower are set to zero." << endl;
	cout << "-log" << endl;
	cout << "\tReplaces the distances with their logs: dst=ln(dst)+1 if dst > 1." << endl;
	cout << "-noheader" << endl;
	cout << "\tOmits the header from the output file (brief and MMML modes only)." << endl;
	cout << "-onlyheader" << endl;
	cout << "\tPrints the header to the output file and exits (brief and MMML modes only)." << endl;
	cout << "-ignoreerrors" << endl;
	cout << "\tDoes not abort the run even in the case of program errors (intended for debugging only)." << endl;
	
	cout << "MMML extensions:" << endl;
	cout << "\t-Lt rootedTreeFileName (required)" << endl;
	cout << "\t\tFile name for a Newick style tree file encompassing all of the taxa to be encountered in the MMM run." << endl;
	cout << "\t-LtopCov (optional)" << endl;
	cout << "\t\ttop-coverage style output." << endl;
	cout << "\t-Li ignoreTaxaFileName (optional)" << endl;
	cout << "\t\tfile listing taxon names to ignore in the MMML analysis." << endl;
	cout << "\t-Lb minBootstrapSupport (optional, default=1)" << endl;
	cout << "\t\tan integer between 0 and 100 representing the threshold below which to collapse branches in parsing the treefile." << endl;
	cout << "\t-List (optional)" << endl;
	cout << "\t\tif specified, lists the taxa in each clade in the output." << endl;
	cout << "\t-Lv (optional)" << endl;
	cout << "\t\tif specified, outputs the clade sets." << endl;
	
	cout << "Alignment processing mode:" << endl;
	cout << "\tPermits the recognition of .aln files as FASTA alignments and processing them into PMB matrices." << endl;
	cout << "\tWhen enabled, simply specify a .aln filename in place of distance matrix file. This works in both command-line and batch modes." << endl;
	cout << "\tA named alignment slice can be defined as comma-separated ranges of columns added after the filename." << endl;
	cout << "\tWhen using multiple ranges to define a slice, the ranges must go in increasing order, and must not overlap." << endl;
	cout << "\tAn \"inverse\"slice (remove the slice itself and keep all other columns) can be specified by adding 'r' after the range. Do not use with multiple ranges." << endl;
	cout << "\tFormat: <filename>|[range name]:<start>-<end>[r][,<start>-<end> ...]" << endl;
	cout << "\tExamples:" << endl;
	cout << "\t\talignment1.aln" << endl;
	cout << "\t\talignment2.aln|pfam12820:373-546" << endl;
	cout << "\t\talignment3.aln|disorder:16-90,150-150,430-519" << endl;
	cout << "\t\talignment4.aln|:300-400r" << endl;
	cout << "\t-aln2pmb"  << endl;
	cout << "\t\tEnables alignment processing mode. Required for the recognition of following options:" << endl;
	cout << "\t-gapthresh alignmentGapThreshold" << endl;
	cout << "\t\tDefines the fraction of gaps (symbols - or X) beyond which sequences will be excluded from alignments (and slices)" << endl;
	cout << "\t-printpmb" << endl;
	cout << "\t\toutputs generated distance matrix files as name.pmbq" << endl;
	cout << "\t-printslice" << endl;
	cout << "\t\toutputs generated alignement slice files as name.aln" << endl;
	
	exit(1);
}
//
	
void Parameters::process(int argc, const char* argv[])	
{	CmdLineArgParser options(argc, argv);
	doShowStats = options.parse("-stats");
	if (options.parse("-help")) printVerboseUsage(options.programName());
	
	bool haveDist1 = options.parse("-d1", &distfileName1);
	bool haveDist2 = options.parse("-d2", &distfileName2);
	bool haveBatch = options.parse("-b", &batchFileName);

	stateRestore = options.parse("-restore", &stateFileName);

	if (!doShowStats || !haveDist1 || !haveDist2)
	{
		if ((!haveBatch && (!haveDist1 || !haveDist2)) ||
			!options.parse("-o", &outfileName) || !options.parse("-a", &allow))
		{
			printUsage(options.programName());
		}
	}

	if (haveBatch && stateRestore)
	{
		printUsage(options.programName());
	}

	brief = options.parse("-brief"); // optional; edit New 2010.03.01: RLC, fixed by abezginov
	tabDelimitedTable = (batchFileName.empty() && options.parse("-table")); // optional, precluded if -b is selected
	if (!batchFileName.empty() && !distfileName1.empty()) printUsage(options.programName());
	tabulateTop = options.parse("-top"); // optional;
	if (tabulateTop && !tabDelimitedTable) printUsage(options.programName());
	useTaxInfo = options.parse("-u"); // optional
	options.parse("-sz", &strictSize); // optional
	if (strictSize != -1 && strictSize < 3) {
		cout << "strictSize must be 3 or more." << endl;
		exit(1);
	}
	if (options.parse("-minsize", &minSize) && minSize < 3) // optional
	{
		cout << "minSize must be 3 or more." << endl;
		exit(1);
	}
	options.parse("-maxtrees", &maxTrees); // optional
	verbose = options.parse("-v"); // optional
	options.parse("-req", &reqTaxName); // optional
	doShowStats = options.parse("-stats");
	options.parse("-p", &distfileNamesPrefix);
	outputZeros = options.parse("-z");
	exactTime = options.parse("-exacttime");
	useHardware = options.parse("-hard");
	options.parse("-lower", &cutoffLower);
	options.parse("-upper", &cutoffUpper);
	logDistances = options.parse("-log"); 
	ignoreErrors = options.parse("-ignoreerrors");
	
	
	reqTax = !reqTaxName.empty();
	if (reqTax && !useTaxInfo)
	{	cout << "\nError: Must use taxon info when requiring a taxon." << endl;
		exit(1);
	}
	
	if (reqTax)
		cout << "Warning: -req option is currently broken (may return lower scores than it should). Use with caution." << endl;
		
		
	// MMML extensions: New 2010.02.28: RLC:
	options.parse("-Lt", &treeFileName); // optional
	if (!treeFileName.empty()) {
		if (brief) printUsage(options.programName()); // incompatible options
		MMMLmode = true;
		// Note: outfileName now refers to MMML output
		options.parse("-Li", &ignoreTaxaFileName); // optional
		options.parse("-Lb", &minBootstrapSupport); // optional
		listTaxa = options.parse("-List"); // optional
		topCoverage = options.parse("-LtopCov"); // optional
		LverboseMode = options.parse("-Lv"); // optional
	}
	// New 2010.02.28: RLC.
		
	// abezgino 2011.01.21
	if (brief)
	{	quickMP = options.parse("-quickmp");
		onlyMP = options.parse("-onlymp");
		if (onlyMP) outputZeros = true;
	}
	
	
	if (brief || MMMLmode)
	{	noHeader = options.parse("-noheader"); 
		onlyHeader = options.parse("-onlyheader");
		if (onlyHeader && noHeader)
		{	cout << "Incompatible options: -noheader -onlyheader" << endl;
			printUsage(options.programName());
		}
	}
	
	
	const int foundRangeStart = options.parse("-rangeStart", &rangeStart); // optional
	const int foundRangeEnd = options.parse("-rangeEnd", &rangeEnd); // optional
	if (foundRangeStart ^ foundRangeEnd) printUsage(options.programName());
	
	aln2pmb = options.parse("-aln2pmb"); 
	if (aln2pmb)
	{	options.parse("-gapthresh", &gapThreshold);
		printPmb = options.parse("-printpmb");
		printSlice = options.parse("-printslice");
		
		options.parse("-alnCacheSize", alnCacheSizeLimit);
		if (alnCacheSizeLimit < 2)
		{	cout << "Alignment cache size must be at least 2" << endl;
			printUsage(options.programName());
		}
		
		options.parse("-matrixCacheSize", distanceMarixCacheSizeLimit);
		if (distanceMarixCacheSizeLimit < 2)
		{	cout << "Matrix cache size must be at least 2" << endl;
			printUsage(options.programName());
		}
		
	}
	
	gzipIn = options.parse("-gzipin"); 
	gzipOut = options.parse("-gzipout"); 
	
	options.parse("-memoryUsageLimit", memoryUsageLimitMegabytes);
	
	
	if (!options.empty())
	{	cout << "Unrecognized option:" << endl;
		options.print();
		cout << endl;
		printUsage(options.programName());
	}
}
//
	
	
