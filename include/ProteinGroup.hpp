#ifndef CLASS_PROTEIN_GROUP
#define CLASS_PROTEIN_GROUP

#include <string>
#include <vector>
#include <sstream>
#include <list>
#include <stdexcept>
#include <typeinfo>


using namespace std;

class ProteinGroup {

protected:
	bool valid_;
	string name_;
	int nOfSequences_;
	vector <string> headers_;
	void init(void){valid_=false; nOfSequences_=0; headers_.clear();}

public:
	ProteinGroup(string const& name){name_ = name; init();}
	~ProteinGroup(){}
	
	//void set_valid(bool valid) {valid_ = valid;}
	//void set_name(string name) {name_ = name;}
		
	void validate(void) const;
	bool is_valid(void) const {return valid_;}
	string get_name(void) const {return name_;}
	int get_nOfSequences() const {validate(); return nOfSequences_;}
		
	string get_header(int index) const {validate(); return headers_.at(index);}
	string get_taxon(int index) const;
	void throw_error(string const& message) const {throw runtime_error(get_name() + ' ' + message);}
	
	
	//doesn't translate to derrived classes: http://yosefk.com/c++fqa/inheritance-proper.html#fqa-21.4
	//static bool find_in_cache(list <ProteinGroup>& cache, string const& name);
};
//
#endif
