#ifndef CLASS_ALIGNMENT
#define CLASS_ALIGNMENT

#include <cstring>
//#include <stdio.h>
#include <iostream>
#include <fstream>
//#include <stdlib.h>
//#include <string>
//#include <vector>
#include <map>
//#include <sstream>

#include "ProteinGroup.hpp"

#define BUFFER_SIZE 65536

using namespace std;

class Alignment : public ProteinGroup {

private:
	
	int nOfColumns_;
	int joinedAt_; //if the alignment is made by joining 2 alignments then it refers to the length of the first alignment
	
	vector <string> sequences_;
	
	void alignment_init(void){nOfColumns_=-1; joinedAt_=0;};
	void init(void){ProteinGroup::init(); alignment_init();}
	static void prepare_sequence(char *buffer);
	static int sequence_count_gaps(string const& sequence);

	//Alignment(int const joinedAt){init(); joinedAt_ = joinedAt;}

public:
	//Alignment(string const& inFilename, double gapThreshold = -1) {init();read_alignment(inFilename, gapThreshold);}
	Alignment(string const& name, int joinedAt = 0) : ProteinGroup(name) {alignment_init(); joinedAt_ = joinedAt;};
	~Alignment(){}
	
	void read_alignment(string const& inFilename, double gapThreshold = -1);
	//void read_alignment(string const& inFilename) {read_alignment(inFilename, -1);}
	
	void add_sequence(string header, string sequence, double gapThreshold = -1);
	Alignment append_alignment(Alignment const& other);
		
	Alignment get_slices(vector<pair<int, int > > const& ranges, double gapThreshold = -1) const;
	Alignment get_slice(int start, int end, double gapThreshold = -1) const {pair <int,int> range (start,end); vector<pair<int, int > > ranges; ranges.push_back(range); return get_slices(ranges, gapThreshold);}
	void print_to_file(string const& outFilepath) const;

	int get_nOfColumns() const {validate(); return nOfColumns_;}
	int get_joinedAt() const {validate(); return joinedAt_;}
	string get_sequence(int index) const {validate(); return sequences_.at(index);}
};
//
#endif
