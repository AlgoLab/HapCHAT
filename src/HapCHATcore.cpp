//include c++ libraries
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <ios>

//include Readset/HapChat libraries
#include "readset.h"
#include "columniterator.h"
#include "basic_types.h"
#include "entry.h"




using namespace std;

//Struct iterator for readset
typedef struct Iterator{
	ColumnIterator it;
	Iterator(ReadSet* readset):it(*readset,nullptr){};
	Iterator():it(*new ReadSet(),nullptr){};

}Iterator;

class HapCHATcore {

	private: 
		Iterator* iterator;
		bool end;
	public:
		//standard constructor
		HapCHATcore(ReadSet* read_set){
		read_set->reassignReadIds();
		this->iterator=new Iterator(read_set);
		end=false;
		};
		HapCHATcore(string filename){
		 ifstream input;
		 int position,allele;
		 unsigned int phred,nread=0;
		 ReadSet* readset =new ReadSet();
		 Read* read;
		 bool flag;
		 stringstream sline;
  		string entry;
  		string token;
		 try {
     			 input.open(filename, ios::in);
    		 } catch(exception & e) { 
      			cerr << "ERROR: failing opening the input file: " << filename << "\": " << e.what() << endl;
      			exit(EXIT_FAILURE);
    		 }

    if(!input.is_open()) {
    		  cerr << "ERROR: failing opening the input file: " << filename << endl;
    	  exit(EXIT_FAILURE);
    	}
    	string line;
    	while(!input.eof()){
    	getline(input, line, '\n');
    	read=new Read(to_string(nread),0,0,0);
    	nread++;
    	if(!line.empty()){
    			sline=stringstream(line);
    			flag=true;
    			while(flag){
    			getline(sline,entry,':');
    			stringstream sentry(entry);
    			sentry >> token;
    			if(!token.compare("#") == 0){
    			position=atoi(token.c_str());
    			sentry >> token;
    			sentry >> token;
					allele = atoi(token.c_str());    			
    			sentry >> token;
    			phred = atoi(token.c_str());
    			read->addVariant(position,allele,phred);
    			}else{flag=false;}
    			}
    		readset->add(read);    		
    	}   	
    	}
    	end=false;
			this->iterator=new Iterator(readset);
			readset->reassignReadIds();
					};
					
		//take the current column 
		Column getColumn(){
			if(hasNext()){
				unique_ptr<vector<const Entry *> > next;
				//if(iterator->it.has_next()) { cout <<"has next" << endl; } // sanity check
			
				next= iterator->it.get_next();
				vector<const Entry*>* p=next.release();
				Column column;
		
				for(unsigned int i=0;i<p->size();i++){
					column.push_back(*p->at(i));
				}
				return column;
			}
				end=true;
    		return Column(0, Entry(-1, Entry::BLANK, 0));
		};
	
		//return true if there is other column
		bool hasNext(){
			return iterator->it.has_next();
		}
		//set the pointer to the first column
		void reset(){
			iterator->it.jump_to_column(0);
			end=false;
		};
		
		//print the readid,allele,quality of each entry of the column(for test)
		void print(Column column){
			cout << "column: ";
			for(unsigned int i=0; i<column.size();i++){
			cout << column[i].get_read_id()<<","<<column[i].get_allele_type()<<","<<column[i].get_phred_score()<<";";
			};
			cout << endl;
		};
		
		unsigned int columnCount(){
			return iterator->it.get_column_count();

		};
		bool isEnded(){
			return end;
		}
};
