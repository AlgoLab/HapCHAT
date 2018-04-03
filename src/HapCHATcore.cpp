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
typedef struct Iterator{
ReadSet p=ReadSet();
ColumnIterator it;
Iterator(ReadSet* readset):it(*readset,nullptr){};
Iterator():it(p,nullptr){};

}Iterator;
class HapCHATcore {

private: 
ReadSet* readset;
Iterator* iterator;
public:
	//HapCHATcore(ReadSet* read_set) : readset(read_set),iterator(readset){};
	HapCHATcore(){
//create a readset to try the class
	readset=new ReadSet();
	Read* read1 = new Read ("0",254,0,0);
	read1->addVariant(9010253,0,4);
	read1->addVariant(9010806,0,10);
	read1->addVariant(9010983,0,5);
	read1->addVariant(9011045,0,14);
	read1->addVariant(9011536,0,14);
	read1->addVariant(9011585,0,13);
	read1->addVariant(9011721,0,13);
	readset->add(read1);

	Read* read2 = new Read ("1",254,0,0);
	read2->addVariant(9010253,1,12);
	read2->addVariant(9010569,1,4);
	read2->addVariant(9010806,1,9);
	readset->add(read2);

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
	readset->add(read3);
	
	Read* read4 = new Read ("3",254,0,0);
	read4->addVariant(9010253,0,13);
	read4->addVariant(9010806,0,9);
	read4->addVariant(9010983,0,14);
	read4->addVariant(9011045,0,11);
	read4->addVariant(9011536,0,12);
	read4->addVariant(9011585,0,11);
	readset->add(read4);
	readset->reassignReadIds();
	Iterator it(readset);
	this->iterator=new Iterator(readset);
};
//take the current column 
//the line commented is where i got the error
Column getColumn(){
unique_ptr<vector<const Entry *> > next;
//if(iterator->it.has_next()) { cout <<"has next" << endl; } // sanity check
next= iterator->it.get_next();
vector<const Entry*>* p=next.release();
Column column;
for(unsigned int i=0;i<p->size();i++){
column.push_back(Entry(
			p->at(i)->get_read_id(),
			p->at(i)->get_allele_type(),
			p->at(i)->get_phred_score()
			));
}

return column;
};
bool hasNext(){
return iterator->it.has_next(); 
}
void reset(){
iterator->it.jump_to_column(0);
};
void print(Column column){

	cout << "column: ";
	for(unsigned int i=0; i<column.size();i++){
	cout << column[i].get_read_id()<<","<<column[i].get_allele_type()<<","<<column[i].get_phred_score()<<";";
	};
	cout << endl;







};
};
