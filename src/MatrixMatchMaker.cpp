// MatrixMatchMaker
// Version 2010.10.31
// (c) 2008-2010, Robert L. Charlebois and Elisabeth R. M. Tillier

// Do not distribute

#include "classCmdLineArgParser.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "mmm_algorithm.h"
#include "timer.h"
#include "sighand.h"
#include "MMML.h"

#include "Alignment.hpp"
#include "DistanceMatrix.hpp"

#include "Parameters.hpp"

#include "getRSS.h"




using namespace std;
using namespace boost::iostreams;

struct CmpTreesByWeight
{
	bool operator()(const WeightedProteinPairVec& a, const WeightedProteinPairVec& b)
	{
		return a.second > b.second;
	}
};



static size_t findBiggestPossibleMatch(const vector<string>& taxLabels1, const vector<string>& taxLabels2,
									   const bool useTaxInfo);
static void processTaxonLabels(vector<string> taxLabels, const string reqTaxName, vector<int>& taxIndices,
										 vector<string>& taxLabelsGlobal, int& reqTaxIndexGlobal);
static void showStats(vector<int>& taxIndices1, vector<int>& taxIndices2, vector<string>& taxLabels);
static int calcTaxonCombinations(vector<int>& taxIndices1, vector<int>& taxIndices2);
static string makeStateSaveFileName(const string& first, const string& second);
static void processDistances(vector<double>& distances, double cutoffLower, double cutoffUpper, bool logDistances);
static void handleErrors(bool ignoreErrors, filtering_stream<output>& outfile);
static void rearrangeMatrix(vector<string>& names, vector<string>& taxLabels, vector<double>& distances, const string reqTaxName);



static unsigned int getDistMatrix(const string& which, const string& distfileName,
									  list <DistanceMatrix>& distMatricesCache, list <Alignment>& alignmnentsCache,
									  vector <double>& distances, vector <string>& names, vector <string>& taxLabels,
									  bool useTaxInfo, Parameters params);
static Alignment const& getAlignment(const string& inFileName, list <Alignment>& alignmnentsCache, Parameters params);
static bool find_in_alignments_cache(list <Alignment>& cache, string const& name);
static bool find_in_dist_matrices_cache(list <DistanceMatrix>& cache, string const& name);
static Alignment getSlice(stringstream& distFilenameSS, Alignment const& alignment, double gapThreshold);


int main(int argc, const char* argv[])
{


	cout << "Started" << endl;
	// Read arguments from the command line:
	Parameters params;
	params.process(argc, argv);
	
	cout << "params done" << endl;
	
	//populate the list of pairs of distance matrices to run
	vector<pair<string, string> > theDistFileNames;
	if (params.batchFileName.empty())
	{	theDistFileNames.push_back(pair<string, string>(params.distfileName1, params.distfileName2));
	} 
	else
	{	filtering_stream<input> batchFile;
		if (params.gzipIn)
		{	batchFile.push(gzip_decompressor());
		}
		try
		{	batchFile.push(file_descriptor_source(params.batchFileName));
		} catch (exception& e)
		{	cout << "\nError: cannot open file " << params.batchFileName << "\nReason:" << e.what() << endl;
			exit(1);
		}
		
		string thePair;
		string::size_type p;
		getline(batchFile, thePair, '\n');
		while (batchFile)
		{	p = thePair.find('\t');
			if (p == string::npos) break;
			theDistFileNames.push_back(pair<string, string>(thePair.substr(0, p), thePair.substr(p + 1)));
			getline(batchFile, thePair, '\n');
		}

		batchFile.strict_sync();
		batchFile.reset();
	}

	// MMML extensions: New 2010.02.28: RLC:
	const set<string> ignoreTaxa(identifyTaxaToIgnore(params.ignoreTaxaFileName));
	
	const pair<TreeType, map<string, set<string> > > theClades(
	retrieveRefCladesFromFile(params.treeFileName, params.minBootstrapSupport));
	const TreeType& refSets = theClades.first; // for convenience
	const map<string, set<string> >& refClades = theClades.second;
	// These are empty if we're not in MMML mode...
	// MMML extensions: New 2010.02.28: RLC.

	ofstream tablefile; // New 2010.01.23: RLC edit
	filtering_stream<output> outfile;


	if (params.gzipOut)
	{	outfile.push(gzip_compressor());
	}
	try
	{	outfile.push(file_descriptor_sink(params.outfileName.c_str()));
	} catch (exception& e)
	{	cout << "\nError: cannot open file " << params.outfileName << "\nReason:" << e.what() << endl;
		exit(1);
	}



	
	// MMML extensions: New 2010.02.28: RLC:
	if (params.MMMLmode)
	{	if (params.MMMLmode && params.LverboseMode) writeTree(outfile, refSets, params.treeFileName);
	}
	// MMML extensions: New 2010.02.28: RLC.
	

	
	
	//abezgino 2011.01.21 : print header when in batch mode
	if (!params.noHeader)
	{	if (params.MMMLmode) writeHeader(outfile, params.topCoverage, params.listTaxa);
		if (params.brief)
		{	outfile << "#Matrix1\tMatrix2\ttaxonCombinations";
			if (params.quickMP) outfile <<"\tquickMP";
			outfile << "\tMP";	
			if (!params.onlyMP) outfile << "\ttime\tscore\tRMSD\ttrees";
			outfile << endl;
		}
			
		if (params.onlyHeader)
		{	outfile.strict_sync();
			outfile.reset();
			exit (0);
		}
	}

	if (params.useHardware)
	{
		if (!CHardAlgorithm::initHardware())
		{
			cout << "Could not initialize max clique hardware" << endl;
			return -1;
		}
	}

	list <DistanceMatrix> distMatricesCache;
	list <Alignment> alignmentsCache;
	
// Batch loop:
	vector<pair<string, string> >::const_iterator batchIt, batchItEnd = theDistFileNames.end();
	for (batchIt = theDistFileNames.begin(); batchIt != batchItEnd; ++batchIt)
	{
		// Print matrix file names to stdout
		cout << batchIt->first << '\t' << batchIt->second; 
	
		// Brief format output: start each row with the two matrix file names
		if(params.brief && params.outputZeros) 
		{	
			outfile << batchIt->first << '\t' << batchIt->second;
		}

		// Read and process the distance matrices:
		vector<double> distances1, distances2;
		vector<string> names1, names2, taxLabels1, taxLabels2, taxLabels;

		
		int num1, num2;
		try
		{	num1 = getDistMatrix("1:", batchIt->first, distMatricesCache, alignmentsCache, distances1, names1, taxLabels1, true, params);
			num2 = getDistMatrix("2:", batchIt->second, distMatricesCache, alignmentsCache, distances2, names2, taxLabels2, true, params);
		} catch (exception &e)
		{	cout << "\tError: Could not obtain at least one distance matrix: " << e.what() << endl;
			cerr << "\tError: Could not obtain at least one distance matrix: " << e.what() << endl;
			if(params.brief && params.outputZeros) 
			{	outfile << endl;
			}
			continue;
		}
		
		if (params.reqTax)
		{	rearrangeMatrix(names1, taxLabels1, distances1, params.reqTaxName);
			//rearrangeMatrix(names2, taxLabels2, distances2, reqTaxName);

			//reqTax=false; //debug - check whether simply rearranging the matrices affects the results
		}
		#ifdef VERBOSE
		cout << "reqTax done" << endl;
		#endif
		// Assign unique integers to taxon labels
		vector<int> taxIndices1, taxIndices2;
		int reqTaxIndex;
		processTaxonLabels(taxLabels1, params.reqTaxName, taxIndices1, taxLabels, reqTaxIndex);
		processTaxonLabels(taxLabels2, params.reqTaxName, taxIndices2, taxLabels, reqTaxIndex);
		
		#ifdef VERBOSE
		cout << "processTaxonLabels done" << endl;
		#endif
		
		// showstats: just show information about matrices and do nothing
		if (params.doShowStats)
		{	showStats(taxIndices1, taxIndices2, taxLabels);
			continue;
		}
	
		processDistances(distances1, params.cutoffLower, params.cutoffUpper, params.logDistances);
		processDistances(distances2, params.cutoffLower, params.cutoffUpper, params.logDistances);

		#ifdef VERBOSE
		cout << "processDistances done" << endl;
		#endif
		
		// Initialize the algorithm
		CMMMAlgorithm* algo = NULL;
		try 
		{	algo = new CHardAlgorithm(params.allow, distances1, distances2,
				taxIndices1, taxIndices2, params.useTaxInfo, reqTaxIndex, 
				params.reqTax, params.strictSize, params.maxTrees, params.minSize, params.useHardware);
		}catch (exception &e)
		{	cout << "\tError: Could not initialize the algorithm: " << e.what() << endl;
			cerr << "\tError: Could not initialize the algorithm: " << e.what() << endl;
			if(params.brief && params.outputZeros) 
			{	outfile << endl;
			}
			sighand_restore();
			delete algo;
			continue;
		}
		
		#ifdef VERBOSE
		cout << "algo done" << endl;
		#endif
		
		// Calculate match potential and taxon combinations
		int taxonCombinations = calcTaxonCombinations(taxIndices1, taxIndices2);
		
		#ifdef VERBOSE
		cout << "calcTaxonCombinations done" << endl;
		#endif
		
		//abezgino 2011.01.21: calculate both match potentials (should not add too much overhead)
		int quickMatchPotential = findBiggestPossibleMatch(taxLabels1, taxLabels2, params.useTaxInfo);
		
		#ifdef VERBOSE
		cout << "findBiggestPossibleMatch done" << endl;
		#endif
		
		int matchPotential = algo->calcMatchPotential();
		#ifdef VERBOSE
		cout << "calcMatchPotential done" << endl;
		#endif
		
		//abezgino 2011.01.21
		//debug code to ensure that new matchPotential is never larger than quickMatchPotential
		if (matchPotential > quickMatchPotential)
		{	matchPotential = quickMatchPotential;
			#ifdef VERBOSE
			cout << batchIt->first << '\t' << batchIt->second << "\tError: matchPotential=" << matchPotential << " > quickMatchPotential=" << quickMatchPotential;
			cerr << batchIt->first << '\t' << batchIt->second << "\tError: matchPotential=" << matchPotential << " > quickMatchPotential=" << quickMatchPotential;
			handleErrors(params.ignoreErrors, outfile);
			#endif
	
		}
		// Brief format output: print taxon combinations and match potential before analysis begins
		if(params.brief && params.outputZeros)
		{
			//abezgino: print known information about current combination before the actual analysis. Make taxonCombinations optional? 
			//abezgino 2011.01.21
			outfile << '\t' << taxonCombinations;
			if (params.quickMP) outfile << '\t' << quickMatchPotential;
			outfile << '\t' << matchPotential;
			if (params.onlyMP)
			{	cout << endl;
				outfile << endl;
				continue;
			} else
			{	cout << flush;
				outfile << flush;
			}
		}
		
		
		
		// Check that the required taxon name, if specified, is in both matrices:
		if (	params.reqTax &&
				(	find(taxLabels1.begin(), taxLabels1.end(), params.reqTaxName) == taxLabels1.end() ||
					find(taxLabels2.begin(), taxLabels2.end(), params.reqTaxName) == taxLabels2.end()
				)
			)
		{	if (params.brief && params.outputZeros) outfile << "\t0\t0\t0\t0" << endl; //format: time	score	rmsd	trees
			delete algo;
			continue; // New 2010.01.23: RLC bug fix
		}

		// Conduct the analysis:
		Timeval startTime, endTime;
		WeightedProteinPairVecVec resultVec;

		int wkRangeStart = max(0, min(params.rangeStart, num1));
		int wkRangeEnd = min(params.rangeEnd, num1);
		int maxScore = 0;

		timer_gettime(&startTime);
	
		if (params.stateRestore)
		{
			if (!algo->stateRestore(params.stateFileName))
			{	cout << "\tError loading state file " << params.stateFileName << endl;
				exit(1);
			}
		}
		
		string saveFile = makeStateSaveFileName(batchIt->first, batchIt->second);
		sighand_init(algo, saveFile);
		
		try
		{
			maxScore = algo->launchAnalysis(wkRangeStart, wkRangeEnd);	
			
			if (maxScore < 0)
			{
				sighand_restore();
				cout << "\tLast reached index: " << algo->lastReachedIndex() << endl;
				delete algo;
				continue;
			}
			algo->calcScores(resultVec);
			algo->printBonusStats(outfile, params.brief);
		}
		catch (bad_alloc&)
		{
			cout << "\tOut of memory (bad_alloc)" << endl;
			sighand_restore();
			delete algo;
			continue;
		}
		
		timer_gettime(&endTime);
		double runtimeExact = timer_diff(&startTime, &endTime);
		double runtime = (params.exactTime) ? runtimeExact : floor(runtimeExact + 0.5);
		cout << "\tTime: " << runtime;

		if (params.maxTrees >= 0 && algo->maxTreesReached())
		{
			cout << "\tWarning: result was truncated by -maxTrees to " << params.maxTrees;
		}
		cout << endl;
		
		sighand_restore();
		delete algo;

		if (params.outputZeros && maxScore == 0) 
		{	if (params.brief)
			{	outfile << '\t' << runtime << "\t0\t0\t0" << endl; //format:time	score	rmsd	trees
			} else if (params.MMMLmode)
			{	MMMdata* theMMMdata = new MMMdata(batchIt->first, batchIt->second, params.distfileNamesPrefix, 0, ignoreTaxa, refClades, params.listTaxa);
				if (theMMMdata->bad())
					cout << "Error constructing MMMdata object for " << batchIt->first << " - " << batchIt->second << endl;
				else
					theMMMdata->produceSummaryOutput(outfile);
				delete theMMMdata;
			}
			continue;
		}
		
		
		
		// Generate the output files:
		if ((params.strictSize < 3 && maxScore > 0) || (maxScore == params.strictSize))
		{
			if (maxScore < 3)
			{
				//cout << "Warning: maxScore=" << maxScore << "< 3" << flush;
				cerr << batchIt->first << '\t' << batchIt->second << "\tError: maxScore=" << maxScore << " < 3";
				handleErrors(params.ignoreErrors, outfile);
			}
			
			const size_t numTrees = resultVec.size();
			
			sort(resultVec.begin(), resultVec.end(), CmpTreesByWeight());
			
			// Debug code: returned trees should not be larger than biggest possible match
			if (maxScore > matchPotential)
			{
				cerr << batchIt->first << '\t' << batchIt->second << "\tError: maxScore=" << maxScore << " > matchPotential=" << matchPotential;
				handleErrors(params.ignoreErrors, outfile);
			}
			//assert(maxScore <= matchPotential);
			
			// Output summary
			if (params.brief)
			{
				if (!params.outputZeros)
				{	
					outfile << batchIt->first << '\t' << batchIt->second << '\t' << taxonCombinations; 
					if (params.quickMP) outfile << '\t' << quickMatchPotential;
					outfile << '\t' << matchPotential;
				}
				
				double rmsd = 0.0;
				if (!resultVec.empty())
				{
					rmsd = resultVec.front().second;
				}
				outfile << '\t' << runtime << '\t' << maxScore << '\t' << rmsd << '\t' << numTrees << endl;
				if (!params.tabDelimitedTable) continue; //nothing more to do in brief mode except when in table mode. Do we actually need table mode at all now?
			}
			else if (!params.MMMLmode) // New 2010.02.28: RLC
			{
				outfile << "infile 1: " << batchIt->first << '\n'; // New 2010.01.23: RLC, to support MMML
				outfile << "infile 2: " << batchIt->second << '\n'; // New 2010.01.23: RLC, to support MMML
				outfile << "maxScore = " << maxScore << " (max possible score = " << matchPotential << ")\n";
				outfile << "number of matching trees = " << numTrees << '\n';
			}

			// MMML extensions: New 2010.02.28: RLC:
			MMMdata* theMMMdata = 0;
			if (params.MMMLmode) {
				theMMMdata = new MMMdata(batchIt->first, batchIt->second, params.distfileNamesPrefix, numTrees, ignoreTaxa, refClades, params.listTaxa);
				if (theMMMdata->bad()) {
					cout << "Error constructing MMMdata object for " << batchIt->first << " - " << batchIt->second << endl;
					delete theMMMdata;
					continue;
				}
			}
			// MMML extensions: New 2010.02.28: RLC.

			// Output the trees
			WeightedProteinPairVecVec::const_iterator vmIt, vmItEnd = resultVec.end();
			ProteinPairVec::const_iterator mIt, mItEnd;
			int t = 1;
			vector<vector<float> > counts(num1, vector<float>(num2));
			for (vmIt = resultVec.begin(); vmIt != vmItEnd; ++vmIt, ++t)
			{
				mItEnd = vmIt->first.end();
				// New 2010.02.28: RLC, modified to feed MMML where desired
				if (params.MMMLmode)
				{	set<string> theTaxa1, theTaxa2;
					for (mIt = vmIt->first.begin(); mIt != mItEnd; ++mIt)
					{	theTaxa1.insert(names1[mIt->first].substr(2));
						theTaxa2.insert(names2[mIt->second].substr(2));
					}
					theMMMdata->processSubmatricesFromMMM(theTaxa1, theTaxa2, t, (params.topCoverage ? 0 : &outfile));
				} else
				{	// New 2010.02.28: RLC.
					double theWtScore = vmIt->second;
					if (!params.brief) outfile << "Submatrix " << t << ": weighted score = " << theWtScore << '\n'; // New 2010.01.23: RLC edit
					for (mIt = vmIt->first.begin(); mIt != mItEnd; ++mIt)
					{	if (!params.brief) outfile << names2[mIt->second] << "->" << names1[mIt->first] << '\n';
						if (t == 1 || !params.tabulateTop) counts[mIt->first][mIt->second] += 1.0F;
					}
				} // New 2010.02.28: RLC
			}
			// MMML extensions: New 2010.02.28: RLC:
			if (params.MMMLmode)
			{	theMMMdata->produceSummaryOutput(outfile);
				delete theMMMdata;
			}
			// MMML extensions: New 2010.02.28: RLC.

			if (params.tabDelimitedTable)
			{
				if (!tablefile.is_open()) // New 2010.01.23: RLC bug fix: as it was, the file kept getting overwritten
				{
					tablefile.open((params.outfileName + ".tab").c_str());
				}
				tablefile << numTrees << '\t' << maxScore << "\n----";
				for (t = 0; t < num1; ++t) 
				{
					tablefile << '\t' << names1[t];
				}
				
				const size_t numTreeNorm = params.tabulateTop ? 1 : numTrees;
				for (int u = 0; u < num2; ++u)
				{
					tablefile << '\n' << names2[u];
					for (t = 0; t < num1; ++t) tablefile << '\t' << counts[t][u] / numTreeNorm;
				}
				tablefile << endl;
				// New 2010.01.23: RLC bug fix: moved closing of file out of loop lest it be repeatedly overwritten
			}
		}
	}

	if (params.useHardware)
	{
		CHardAlgorithm::closeHardware();
	}

	tablefile.close(); // New 2010.01.23: RLC bug fix
	outfile.strict_sync();
	outfile.reset();
	
	return 0;
}
//





//abezgino
//rearrange the matrices to ensure that the sequences for the reqTax are placed at the top
static void rearrangeMatrix(vector<string>& names, vector<string>& taxLabels, vector<double>& distances, const string reqTaxName)
{	
	unsigned int size = taxLabels.size();
	vector<unsigned int> map; //mapping of original indices to the rearranged ones
	map.reserve(size);
			
	unsigned int firstNonReqTaxIdx = 0;
	
	//rearrange taxLabels and names
	for (unsigned int i = 0; i < size; i++)
	{	map.push_back(i);
		if (taxLabels[i] == reqTaxName)
		{	swap(map[i], map[firstNonReqTaxIdx]);
			swap(taxLabels[i], taxLabels[firstNonReqTaxIdx]);
			swap(names[i], names[firstNonReqTaxIdx]);
			firstNonReqTaxIdx++;
		}
	}
	
	//now rearrange distances
	//use the new in-place approach
	for (unsigned int i = 0; i < size; i++)
	{	for (unsigned int j = 0; j < size; j++)
		{	//make sure that each element gets swapped only once
			if(map[i] > i || (map[i] == i && map[j] > j))
			{ swap(distances[i*size+j], distances[map[i]*size+map[j]]);
			}
		}
	}
	
	return;
}
//

//abezgino
//convert taxon labels to integers for easier manipulation
static void processTaxonLabels(vector<string> taxLabels, const string reqTaxName, vector<int>& taxIndices,
										 vector<string>& taxLabelsGlobal, int& reqTaxIndexGlobal)
{	for (unsigned int i = 0; i < taxLabels.size(); i++)
	{
		bool found = false;
		
		unsigned int j;
		for (j = 0; j < taxLabelsGlobal.size(); j++)
		{	if (taxLabels[i] == taxLabelsGlobal[j])
			{	found = true;
				taxIndices.push_back(j);
				break;
			}
		}
		
		if (!found)
		{	//since the above loop was allowed to run to the end then j == taxLabelsGlobal.size();
			taxLabelsGlobal.push_back(taxLabels[i]);
			taxIndices.push_back(j);
			if (reqTaxName == taxLabels[i])
				reqTaxIndexGlobal = j;
		}
	}
}
//

size_t findBiggestPossibleMatch(const vector<string>& taxLabels1, const vector<string>& taxLabels2,
								const bool useTaxInfo)
{
	const size_t num1 = taxLabels1.size();
	const size_t num2 = taxLabels2.size();

	if (!useTaxInfo) return min(num1, num2);

	map<string, pair<size_t, size_t> > taxCounts;
	size_t i;

	for (i = 0; i < num1; ++i)
	{
		++taxCounts[taxLabels1[i]].first;
	}
	for (i = 0; i < num2; ++i)
	{
		++taxCounts[taxLabels2[i]].second;
	}

	size_t answer = 0;
	map<string, pair<size_t, size_t> >::const_iterator mIt, mItEnd = taxCounts.end();
	for (mIt = taxCounts.begin(); mIt != mItEnd; ++mIt)
	{
		answer += min(mIt->second.first, mIt->second.second);
	}
	return answer;
}
//

static void showStats(vector<int>& taxIndices1, vector<int>& taxIndices2, vector<string>& taxLabels)
{
	cout << "Matrix 1 size: " << taxIndices1.size() << endl;
	cout << "Matrix 2 size: " << taxIndices2.size() << endl;
	cout << "Number of unique taxons: " << taxLabels.size() << endl;

	unsigned int compat = calcTaxonCombinations(taxIndices1, taxIndices2);
	unsigned int total = taxIndices1.size() * taxIndices2.size();
	float ratio = (float)compat/(float)total;

	cout << "Compatible/total pairings: " << compat << '/' << total << " (" << ratio << ")" << endl;
}
//

static int calcTaxonCombinations(vector<int>& taxIndices1, vector<int>& taxIndices2)
{
	//abezgino: calculate the number of combinations for all taxon names.
	map<int, int> taxCounts1, taxCounts2;
	int taxonCombinations = 0;

	for (vector<int>::iterator it = taxIndices1.begin();  it < taxIndices1.end(); it++ )
		taxCounts1[*it]++;

	for (vector<int>::iterator it = taxIndices2.begin();  it < taxIndices2.end(); it++ )
		taxCounts2[*it]++;

	for (map<int, int>::iterator it = taxCounts1.begin();  it != taxCounts1.end(); it++ )
	{	
		int currentTaxCount1 = it->second;
		int currentTaxCount2 = taxCounts2[it->first];

		taxonCombinations += currentTaxCount1 * currentTaxCount2;
	}

	return taxonCombinations;
}
//

static string makeStateSaveFileName(const string& first, const string& second)
{
	string::size_type slashpos1, slashpos2;
	string result;

	slashpos1 = first.find_last_of("/\\") + 1;
	slashpos2 = second.find_last_of("/\\") + 1;
	result = first.substr(slashpos1) + '_' + second.substr(slashpos2) + ".sav";

	return result;
}
//

//process a raw distance matrixs 
//for now only remove distances above and below specified thresholds
static void processDistances(vector<double>& distances, double cutoffLower, double cutoffUpper, bool logDistances)
{	bool errorReported = false;
	for (unsigned int i = 0; i < distances.size(); i++ )
	{
		//cout << distances[i] <<" \t";
		if (distances[i] < cutoffLower) distances[i] = 1.0e-20;
		if (distances[i] > cutoffUpper) distances[i] = 1.0e-20;
		if (abs(distances[i] - 0.000010) < 1.0e-20 && !errorReported)
		{	cout << "\tWarning: distance of 0.000010 found. Forgot to zero protdist distances?";
			errorReported = true;
		}
		//cout << distances[i] <<"\n";
		if (logDistances && distances[i] > 1)
		{	distances[i] = log(distances[i])+1;
		}
	}
	return;
}
//


static void handleErrors(bool ignoreErrors, filtering_stream<output>& outfile)
{	if (ignoreErrors)
	{	cerr << "\tIgnored ..." << endl;
	}
	else
	{	cerr << endl;
		cout << endl;
		if (outfile.good())
		{	outfile << endl;
			outfile.strict_sync();
			outfile.reset();
		}
		exit(1);
	}
	return;
}
//

static unsigned int getDistMatrix(const string& which, const string& distFilename,
									  list <DistanceMatrix>& distMatricesCache, list <Alignment>& alignmnentsCache,
									  vector <double>& distances, vector <string>& names, vector <string>& taxLabels,
									  
									  bool useTaxInfo, Parameters params
									  //bool aln2pmb, bool printPmb, bool printSlice, double gapThreshold
									  
									  )
{
	//parameters
	bool quickPmb=true;// printPmb=false, printSlice=false; 
	
	//try to find the matrix in the cache
	if(!find_in_dist_matrices_cache(distMatricesCache, distFilename))
	{
		#ifdef VERBOSE
		cout << "\nWarning: matrix [" << distFilename << "] is not cached. Trying to read from file. Cache size=" << distMatricesCache.size() << '/' << params.distanceMarixCacheSizeLimit << " (usage=" << ceil(getCurrentRSS() / 1024 / 1024) << " MB)" << endl;
		#endif
		distMatricesCache.push_front(DistanceMatrix(distFilename));
		DistanceMatrix& distanceMatrixTmp = distMatricesCache.front();
		//distanceMatrixTmp.set_name(distFilename);
		
		
		stringstream distFilenameSS(distFilename);
		string filename;	
	
		getline(distFilenameSS, filename, '|');
		bool subAln = !distFilenameSS.eof();
		
		string extension = filename.substr(filename.find_last_of('.')+1);
		
		if (params.aln2pmb && extension == "aln")
		{
			#ifdef VERBOSE
			cout << "\tWarning: extension [" << extension <<"] found. Calculating PMB distnce matrix from alignment ..." << endl;
			#endif
			
			Alignment const& alignment = getAlignment(filename, alignmnentsCache, params);
			#ifdef VERBOSE
			cout << "Alignment length: " << alignment.get_nOfColumns() << endl;
			//cout << alignment.get_sequence(1) << endl;
			#endif
			if(subAln)
			{	Alignment alignmentSlice = getSlice(distFilenameSS, alignment, params.gapThreshold);
				#ifdef VERBOSE
				cout << "Slice length: " << alignmentSlice.get_nOfColumns() << endl;
				//cout << alignmentSlice.get_sequence(1) << endl;
				#endif
				if(params.printSlice)
				{	string filepath = params.distfileNamesPrefix+distFilename+".slice";
					ifstream file(filepath.c_str()); 
					if (!file)
						try
						{	alignmentSlice.print_to_file(filepath);
						}
						catch (exception &e)
						{	cout << "\tWarning: Skipped slice output: " << e.what() << endl;
							cerr << "\tWarning: Skipped slice output: " << e.what() << endl;
						}
					else
						file.close();
				}
				distanceMatrixTmp.create_from_alignment(alignmentSlice, quickPmb);
			} else
			{	distanceMatrixTmp.create_from_alignment(alignment, quickPmb);
			}
			
			if (params.printPmb)
			{	string filepath = params.distfileNamesPrefix+distFilename+".pmbq";
				ifstream file(filepath.c_str()); 
				if (!file)
					try
					{	distanceMatrixTmp.print_to_file(filepath);
					}
					catch (exception &e)
					{	cout << "\tWarning: Skipped pmb output: " << e.what() << endl;
						cerr << "\tWarning: Skipped pmb output: " << e.what() << endl;
					}
				else
					file.close();
			}
		} else
		{
			#ifdef VERBOSE
			if(extension == "aln")
				cout << "\tWarning: extension [" << extension <<"] found. Forgot to add -aln2pmb option?" << endl;
			#endif
			distanceMatrixTmp.read_from_file(params.distfileNamesPrefix+distFilename);
		}
		
		while(distMatricesCache.size() > params.distanceMarixCacheSizeLimit || ceil(getCurrentRSS() / 1024 / 1024) > params.memoryUsageLimitMegabytes)
		{	
			#ifdef VERBOSE
			cout << "\tWarning: cache has too many elements. Removing the last one [" << distMatricesCache.back().get_name() << "]" << endl;
			#endif
			distMatricesCache.pop_back();
		}
	}
	
	//at this point the first element of cache list must be the matrix.
	
	DistanceMatrix& distanceMatrix = distMatricesCache.front();
	
	//fill the output variables
	unsigned int nOfSequences = distanceMatrix.get_nOfSequences();
	names.resize(nOfSequences);
	taxLabels.resize(nOfSequences);
	for (unsigned int i = 0; i < nOfSequences; i++)
	{
		string header = distanceMatrix.get_header(i);
		if (useTaxInfo)
		{	taxLabels[i] = header.substr(0, header.find('|')); // If not found, selects the entire header
		}
		names.at(i) = which + header;
	}
	distanceMatrix.get_distances_as_1d(distances);
	
	return nOfSequences; //can be zero?
}
//

static Alignment const& getAlignment(const string& inFilename, list <Alignment>& alignmnentsCache, Parameters params)
{	
	//if(! Alignment::find_in_cache(alignmnentsCache, inFilename))
	if(!find_in_alignments_cache(alignmnentsCache, inFilename))
	{
		#ifdef VERBOSE
		cout << "\tWarning: alignment [" << inFilename << "] not cached. Trying to read from file. Cache size=" << alignmnentsCache.size() << '/' << params.alnCacheSizeLimit << endl;
		#endif
		alignmnentsCache.push_front(Alignment(inFilename));
		//cout << "name=[" << alignmnentsCache.front().get_name() << "]" << endl;
		//if (alignmnentsCache.front().get_name().empty())
		//	cout << "empty name\n"; 
		//alignmnentsCache.front().set_name(inFilename);
		alignmnentsCache.front().read_alignment(params.distfileNamesPrefix+inFilename, params.gapThreshold); //even if this fails, the alignment name is already cached, so the file won't be attempted to get re-read several times over.

		if(alignmnentsCache.size() > 2 && alignmnentsCache.size() > params.alnCacheSizeLimit)
		{
			#ifdef VERBOSE
			cout << "\tWarning: cache has too many elements. Removing the last one [" << alignmnentsCache.back().get_name() << "]" << endl;
			#endif
			alignmnentsCache.pop_back();
		}
	}
	//at this point the front element of cache list must be the alignment of interest.
	
	return alignmnentsCache.front(); 
}
//


static bool find_in_alignments_cache(list <Alignment>& cache, string const& name)
{	bool found = false;
	unsigned int i = 0;
	for(list <Alignment>::iterator it = cache.begin(); it != cache.end() && !found; it++)
	{	i++;
		if (it->get_name() == name)
		{	found = true;
			#ifdef VERBOSE
			cout << "\talignment [" << name << "] is cached at position " << i << " of " << cache.size() << endl;
			#endif
			//move this element to the front (as the most recently used)
			cache.splice(cache.begin(), cache, it);
		}
	}
	return found;
}
//

static bool find_in_dist_matrices_cache(list <DistanceMatrix>& cache, string const& name)
{	bool found = false;
	unsigned int i = 0;
	for(list <DistanceMatrix>::iterator it = cache.begin(); it != cache.end() && !found; it++)
	{	i++;
		if (it->get_name() == name)
		{	found = true;
			#ifdef VERBOSE
			cout << "\tmatrix [" << name << "] is cached at position " << i << " of " << cache.size() << endl;
			#endif
			//move this element to the front (as the most recently used)
			cache.splice(cache.begin(), cache, it);
		}
	}
	return found;
}
//


static Alignment getSlice(stringstream& distFilenameSS, Alignment const& alignment, double gapThreshold = -1)
{	string name;
	getline(distFilenameSS, name, ':');
	
	#ifdef VERBOSE
	cout << "\tWarning: a sub-alignment found: name=[" << name << "]:";
	#endif
	
	vector<pair<int, int > > ranges;
	
	while (!distFilenameSS.eof())
	{	int start=0, end=0;
		distFilenameSS >> start;
		if(distFilenameSS.peek() != '-')
			throw runtime_error("could not parse range for sub-alignment: " + name);
		distFilenameSS.ignore(1);
		distFilenameSS >> end;

		#ifdef VERBOSE
		cout << "[" << start <<"]-[" << end << "]; ";
		#endif
		
		if(!distFilenameSS.eof() && distFilenameSS.peek() == 'r') //"reverse-slice" mode, where only the slice itself is removed
		{	distFilenameSS.ignore(1);
			#ifdef VERBOSE
			cout << 'r';
			#endif
			if (start > 1)
				ranges.push_back(pair<int, int>(1,start-1));
			if (end < alignment.get_nOfColumns())
				ranges.push_back(pair<int, int>(end+1,alignment.get_nOfColumns()));
		} else
		{	ranges.push_back(pair<int, int>(start,end));
		}
		
		if(!distFilenameSS.eof() && distFilenameSS.peek() == ',')
		{	distFilenameSS.ignore(1);
		} else
		{	break; //finished reading the ranges, ignore the rest (if any)
		}
	}
	#ifdef VERBOSE
	cout << endl;
	#endif
	if (distFilenameSS.fail() || ranges.size() < 1)
		throw runtime_error("could not parse sub-alignment: " + name);
	return alignment.get_slices(ranges, gapThreshold);
}
//
