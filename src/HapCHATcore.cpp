#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>

#include "readset.h"
#include "columniterator.h"
#include "basic_types.h"
#include "entry.h"
using namespace std;

class HapCHATcore {

private: 
ReadSet* readset;
ColumnIterator iterator;
public:
	HapCHATcore(ReadSet* read_set) : readset(read_set),iterator(*readset,nullptr){};
	HapCHATcore():iterator(this->prova(),nullptr){};
//create a readset to try the class
const ReadSet prova(){
	ReadSet* r=new ReadSet();
	Read* read1 = new Read ("0",254,0,0);
	read1->addVariant(9010253,0,4);
	read1->addVariant(9010806,0,10);
	read1->addVariant(9010983,0,5);
	read1->addVariant(9011045,0,14);
	read1->addVariant(9011536,0,14);
	read1->addVariant(9011585,0,13);
	read1->addVariant(9011721,0,13);
	r->add(read1);

	Read* read2 = new Read ("1",254,0,0);
	read2->addVariant(9010253,1,12);
	read2->addVariant(9010569,1,4);
	read2->addVariant(9010806,1,9);
	r->add(read2);

	Read* read3 = new Read ("2",254,0,0);
	read3->addVariant(9010253,1,13);
	read3->addVariant(9010569,1,11);
	read3->addVariant(9010983,0,1);
	read3->addVariant(9011045,0,6);
	read3->addVariant(9011536,1,4);
	read3->addVariant(9011585,1,10);
	read3->addVariant(9011721,1,13);
	read3->addVariant(9012775,1,14);
	read3->addVariant(9013095,0,10);
	read3->addVariant(9013969,1,13);
	r->add(read3);
	
	Read* read4 = new Read ("3",254,0,0);
	read4->addVariant(9010253,0,13);
	read4->addVariant(9010806,0,9);
	read4->addVariant(9010983,0,14);
	read4->addVariant(9011045,0,11);
	read4->addVariant(9011536,0,12);
	read4->addVariant(9011585,0,11);
	r->add(read4);
	this->readset=r;
	const ReadSet e = *r;
	return e;
};
//take the current column 
//the line commented is where i got the error
Column getColumn(){
unique_ptr<vector<const Entry *> > next_input_column;
//next_input_column = iterator.get_next();

Column column;
return column;
};
};
