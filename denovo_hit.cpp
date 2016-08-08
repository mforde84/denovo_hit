#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <memory>

using namespace std;

void create_temp(char *curr_dir){
	
	//create working directories
	
	string sys_string = "mkdir -p " + string(curr_dir) + "/denovo-tmp/secondary;";
	system(sys_string.c_str());

} 

void delete_temp(char *curr_dir){
	
	//delete working directories
	
	string sys_string = "rm -rf " + string(curr_dir) + "/denovo-tmp;";
	system(sys_string.c_str());

} 

vector <string > get_header(char *filename){
	
	//gets first line of vcf file
	
	vector <string > header_vector;
	ifstream vcf_file(filename);
	string line;
	if (vcf_file.is_open()){
		getline(vcf_file,line);
		istringstream iss(line);
		string token;
		while(getline(iss,token,'\t')){
			header_vector.push_back(token);
		}
	}
	vcf_file.close();
	return header_vector;

}

vector <ofstream *> create_tmp_files(char *curr_dir, vector <string > &vcf_header){
	
	//creates vector of unique streams for output to temp files
	
	vector <ofstream *> streams;
	streams.resize (vcf_header.size());

	for (int i=0; i < streams.size(); i++){
		streams[i] = new ofstream();
	}
	
	string string_curr = string(curr_dir);
	string open_file = string_curr + "/denovo-tmp/calls.txt";
	streams[0]->open(open_file.c_str());
	
	int file_int=1;
	for (int x = 9; x < vcf_header.size(); x++){
		open_file = string_curr + "/denovo-tmp/" + vcf_header[x] + ".txt";
		streams[file_int]->open(open_file.c_str());
		file_int++;
	}
	
	return streams;
	
}

void close_streams(vector <ofstream *> streams){

	//closes output streams
	
	for (int i = 0; i < streams.size(); i++){
		(*streams[i]).close();
	}

}

void dump_streams(vector <ofstream *> streams, char *filename){

	//splits vcf files to individual temp files lists cell string and line number in vcf
	
	ifstream vcf_file(filename);
	string line;

	while(getline(vcf_file,line)){
		int call_desc = 0;
		istringstream iss(line);
		string token;
		while(getline(iss,token,'\t')){
			int index_pos = call_desc - 8;
			if (call_desc <= 8){
				if (call_desc == 8){
					*streams[0] << token << "\n";
				}else{
					*streams[0] << token << "\t";
				}
			}else{
				*streams[index_pos] << token << "\n";
			}
			call_desc++;
		}	
	}
	vcf_file.close();

}

void spot_denovo(char* curr_dir, char *filename){

	//iterates through pedigree file and writes out lines where the child has a non-null string when both parents are null
	
	ifstream child;
	ifstream p1;
	ifstream p2;
	ifstream ped_file(filename);
	string indexpath = string(curr_dir) + "/denovo-tmp/index.txt";
	string child_line, p1_line, p2_line, ped_line, token;
	
	while(getline(ped_file,ped_line)){
		
		istringstream iss(ped_line);
		int stream_i=0;
		int line_count=1;
		string outpath = string(curr_dir) + "/denovo-tmp/secondary/";
		string outchild = "";
		
		while(getline(iss,token,'\t')){
			string sys_string = string(curr_dir) + "/denovo-tmp/" + token + ".txt";
			if (stream_i==0){
				child.open(sys_string.c_str());
				outchild += token + "_";
			}else if (stream_i==1){
				p1.open(sys_string.c_str());
				outchild += token + "_";
			}else{
				p2.open(sys_string.c_str());
				outchild += token + ".txt";
			}
			stream_i++;
		}
		outpath += outchild;
	
		ofstream hit_stream(outpath.c_str());
		while(getline(child,child_line)){
			getline(p1,p1_line);
			getline(p2,p2_line);
			if ((child_line != ".") and (p1_line == ".") and (p2_line == ".")){
				hit_stream << line_count << "\t" << child_line << "\n";
			}
			line_count++;
		}
		
		hit_stream.close();
		child.close();
		p1.close();
		p2.close();
		
	}

}

void print_denovo(char* curr_dir){

	//combines all hits, and outputs the corresponding line in the call file along with the pedigree information and call stats
	
	string sys_string = string(curr_dir) + "/denovo-tmp/";
	string secondary_string = sys_string + "secondary/";
	string calls_string = sys_string + "calls.txt";
	string raw_calls = sys_string + "raw_denovo.txt";
	string combo_string = sys_string + "combination.awked";
	
	string basher = "for f in $(find " + secondary_string + " -name \"*.txt\"); do sed \"s|^|$f\\t|g\" $f | sed \"s/.*secondary//g\" | sed \"s/.txt//g\" >> ./denovo-tmp/combination.awked; done;";
	system(basher.c_str()); 
	
	string call_line;
	vector <string > calls_holder;
	ifstream calls(calls_string.c_str());
	while(getline(calls,call_line)){
		calls_holder.push_back(call_line);
	}
	
	string hit_line;
	vector <string > hit_holder;
	ofstream raw_hits(raw_calls.c_str());
	ifstream hits(combo_string.c_str());	
	while(getline(hits,hit_line)){
		string token;
		istringstream iss(hit_line);
		int count=0;
		string string_hold[3];
		while(getline(iss,token,'\t')){
			string_hold[count] = token;
			count++;
		}
		int lookup = stoi(string_hold[1]) - 1;
		raw_hits << calls_holder[lookup] << "\t" << string_hold[0] << "\t" << string_hold[2] << "\n";
	}
	raw_hits.close();
	hits.close();

}

void process_denovo (char* curr_dir){

	//sort and aggregate calls

	int hits = 1;
	int count = 0;
	string sort_line;
	string compare_hold[11];
	string compare_token;
	string sys_string = string(curr_dir) + "/denovo-tmp/";
	string sort_calls = sys_string + "raw_denovo.sort.txt";
	string agg_calls = string(curr_dir) + "/potential_denovo.txt";
	string t1, t2, ped_agg, ped_stats;
	
	string basher = "sort denovo-tmp/raw_denovo.txt > denovo-tmp/raw_denovo.sort.txt;";
	system(basher.c_str()); 

	ofstream agg_hits(agg_calls.c_str());
	ifstream sort_hits(sort_calls.c_str());
	getline(sort_hits,sort_line);
	istringstream compare_iss(sort_line);
	while(getline(compare_iss,compare_token,'\t')){
		compare_hold[count] = compare_token;
		count++;
	}
	count = 0;
	ped_agg = compare_hold[9];
	ped_stats = compare_hold[10];
	t1 = compare_hold[0] + compare_hold[1] + compare_hold[2] + compare_hold[3] + compare_hold[4];
	
	agg_hits << "CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT" << "\t" << "PEDIGREE" << "\t" << "CALLSTATS" << "\t" << "COUNTS" << "\n";
	
	while(getline(sort_hits,sort_line)){
		string token;
		istringstream iss(sort_line);
		string string_hold[11];
		while(getline(iss,token,'\t')){
			string_hold[count] = token;
			count++;
		}
		count = 0;
		t2 = string_hold[0] + string_hold[1] + string_hold[2] + string_hold[3] + string_hold[4];
		if (t1 == t2){
			hits++;
			ped_agg = ped_agg + "~" + string_hold[9]; 
			ped_stats = ped_stats + "~" + string_hold[10];
		} else {
			agg_hits << compare_hold[0] << "\t" << compare_hold[1] << "\t" << compare_hold[2] << "\t" << compare_hold[3] << "\t" << compare_hold[4] << "\t" << compare_hold[5] << "\t" << compare_hold[6] << "\t" << compare_hold[7] << "\t" << compare_hold[8] << "\t" << ped_agg << "\t" << ped_stats << "\t" << hits << "\n"; 
			t1 = string_hold[0] + string_hold[1] + string_hold[2] + string_hold[3] + string_hold[4];
			hits = 1;
			ped_agg = string_hold[9];
			ped_stats = string_hold[10];
		}
		for (int i = 0; i < 11; i++){
			compare_hold[i] = string_hold[i];
		}
	}
	agg_hits.close();
	sort_hits.close();
	
}

int main (int argc, char *argv[]){
	
	//get current directory
	char *curr_dir = get_current_dir_name();
	
	//clean up any previous temp files
	delete_temp(curr_dir);
	
	//create temp directory
	create_temp(curr_dir);
	
	//get headers
	vector <string > vcf_header = get_header(argv[1]);
	
	//create temp files
	vector <ofstream *> tmp_files = create_tmp_files(curr_dir,vcf_header);
	
	//split vcf to temp files and close streams
	dump_streams(tmp_files, argv[1]);
	close_streams(tmp_files);
	
	//output all denovo hits to temp file
	spot_denovo(curr_dir, argv[2]);
	
	//merge and output call stats
	print_denovo(curr_dir);
	
	//sort and output raw and unique calls
	process_denovo(curr_dir);
	
	//clean up temp files
	delete_temp(curr_dir);
	
}		
