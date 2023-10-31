#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <sstream>
#include <cstring>
#include <chrono>
#include <experimental/filesystem>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

// compile: g++ checkSNPs.cc -std=c++11 -o checkSNPs

/* Input: two SAM files and two .bed files
.bed file:
chromosome mismatch_location1 mismatch_location1+1 base1 base2

First parse the .bed file on chromosome x to find a list:
a) locations on chromosome x within alignment that should have a SNP
Note: MD and alignment are all based on REFERENCE



two SAM files aligned with minimap2
1st SAM: UC96 reads aligned to Salinas
2nd SAM: UC96 reads aligned to UC96

READNAME CHR_ALIGNED START END MMLOC, MMLOC, ...
*/

map<string, int> chromosome_names = {
{"SALH_C_Chr_9", 9},
{"SALH_C_Chr_8", 8},
{"SALH_C_Chr_7", 7},
{"SALH_C_Chr_6", 6},
{"SALH_C_Chr_5", 5},
{"SALH_C_Chr_4", 4},
{"SALH_C_Chr_3", 3},
{"SALH_C_Chr_2", 2},
{"SALH_C_Chr_1", 1},
{"SERH_U_Chr_9", 9},
{"SERH_U_Chr_8", 8},
{"SERH_U_Chr_7", 7},
{"SERH_U_Chr_6", 6},
{"SERH_U_Chr_5", 5},
{"SERH_U_Chr_4", 4},
{"SERH_U_Chr_3", 3},
{"SERH_U_Chr_2", 2},
{"SERH_U_Chr_1", 1}
};

//static map<string, int>* chromosome_names_directory[] = {&chromosome_names_1, &chromosome_names_2};

static int chromosome_sizes[] = { // (measured in millions of bases)
258,
242,
332,
419,
380,
211,
215,
351,
235
};


string chr;
string prev_chr;
int loc;
int integer;
int current_chr_num = 0;
char c;

int numPrints = 0;

string readname;
int start;
int end;

int total_Salinas = 0;
int total_US96UC23 = 0;
int total_Recombinant = 0;
int total_Bad_Recombinant = 0;
int total_Discarded = 0;
int total_Contradict = 0;
int total_Bad_Mismatch_Number = 0;

int alignment_outcomes[] = {-1, -1};

void countSalinas() {
	cout << "Read is SALINAS" << endl;
	total_Salinas++;
	return;
}
void countUS96UC23() {
	cout << "Read is US96UC23" << endl;
	total_US96UC23++;
	return;
}

int main(int argc, char* argv[]) {

	auto begin = chrono::high_resolution_clock::now();
	// Note: change to 2 input requirements

	if(argc < 2) {
		//cout << "Usage: ./checkSNPS snp_file.bed samfile1.SAM samfile2.SAM" << endl;
		cout << "Usage: ./checkSNPS snp_file1.bed snp_file2.bed samfile1.SAM samfile2.SAM" << endl;
		cout << "    snp_file1: align US to Salinas" << endl;
		cout << "    samfile1: align US to Salinas" << endl;
		exit(1);
	}

	cout << "1. Allocating arrays..." << endl;
	static bool** mismatch_locations_1 = (bool**) malloc(9 * sizeof(bool*));
	static bool** mismatch_locations_2 = (bool**) malloc(9 * sizeof(bool*));
	for(int i = 0; i < 9; i++) {
		mismatch_locations_1[i] = (bool*) malloc(chromosome_sizes[i] * 1000000 * sizeof(bool));
		mismatch_locations_2[i] = (bool*) malloc(chromosome_sizes[i] * 1000000 * sizeof(bool));
		cout << "Successfully allocated arrays for chr: " << i << endl;
		memset(mismatch_locations_1[i], 0, chromosome_sizes[i] * 1000000 * sizeof(bool));
		memset(mismatch_locations_2[i], 0, chromosome_sizes[i] * 1000000 * sizeof(bool));
		// for(int j = 0; j < chromosome_sizes[i] * 1000000; j++) {
		// 	mismatch_locations_1[i][j] = 0;
		// 	mismatch_locations_2[i][j] = 0;
		// }
		//cout << " > finished initializing" << endl;
	}

	// at each location in this array, it states whether or not there is an expected SNP
	static bool** mismatch_locations[] = {mismatch_locations_1, mismatch_locations_2};
	//ifstream snp_file(argv[1]);
	
	cout << "2. Reading in SNP_files..." << endl;
	// Read in BOTH SNP (bed) files
	for(int i = 0; i < 2; i++) {
		current_chr_num = 0;
		numPrints = 0;
		ifstream snp_file(argv[1 + i]);
		if(!snp_file) {
			cout << "SNP_file " << i << " not good" << endl;
			exit(1);
		}
		cout << " Reading SNP file number " << i + 1 << endl;

		while(snp_file) {
			// read in the SNP_file
			snp_file >> chr;
			if(chr != prev_chr) {
				cout << "Found SNP at position (last SNP):" << loc << " in chromosome " << current_chr_num << endl;
				cout << "Reading SNPs for chromosome " << chr << endl;
				current_chr_num ++;
				if(current_chr_num > 9) {
					cout << "Finished reading all 9 chromosomes" << endl;
					break;
				}
				numPrints = 0;
			}

			// here, it matters which order the files are given
			// snp_file >> loc;
			// mismatch_locations_2[current_chr_num - 1][loc + 1] = 1;
			// if(numPrints % 100000 == 0) {
			// 	cout << "Found SNPs at positions: " << loc + 1;
			// }

			// snp_file >> loc;
			// mismatch_locations_1[current_chr_num - 1][loc + 1] = 1;
			// if(numPrints % 100000 == 0) {
			// 	cout << " and " << loc + 1 << " in chromosome " << current_chr_num << endl;
				
			// }

			snp_file >> loc;
			// Should this be after reading the second location??
			mismatch_locations[i][current_chr_num - 1][loc + 1] = 1;
			if(numPrints % 100000 == 0) {
			 	cout << "Found SNPs at positions: " << loc + 1 << " in chromosome " << current_chr_num << endl;
			}
			snp_file >> loc;

			numPrints ++;

			snp_file >> c;
			snp_file >> c;
			prev_chr = chr;
		}
	}

	// PARSE, CHECK ALIGNED READS
	ifstream sam_file_1(argv[3]);
	ifstream sam_file_2(argv[4]);
	if(!sam_file_1 || !sam_file_2) {
		cout << "ERROR: one or both SAM file not good" << endl;
		exit(1);
	}

	string dir_name = argv[3];
	dir_name = dir_name + ".RRSUM";

	mkdir((dir_name).c_str(), 0777);

	ifstream* sam_files[] = {&sam_file_1, &sam_file_2};

	stringstream MDS; // MD string stream

	stringstream OSS; // output string stream
	
	set<int> mm_locations;
	string line;
	int num;

	string name;
	string seq;
	string MD;
	string chromosome;

	string prev_name;

	char base;
	int num_MM;
	int start_pos;
	int end_pos;
	int pos;

	// Taken from mdparser2, parses the MD tag of the SAM file
	for(int i = 0; i < 11; i++) {
		getline(sam_file_1, line);
		getline(sam_file_2, line);
	}


	cout << "3. Processing SAM Reads for both files" << endl;
	cout << dir_name << endl;
	ofstream Recombinant_reads(dir_name + "/Recombinant_Reads.txt");
	ofstream US96UC23_reads(dir_name + "/US96UC23_Reads.txt");
	ofstream Salinas_reads(dir_name + "/Salinas_Reads.txt");
	ofstream AllDiscarded_reads(dir_name + "/AllDiscarded_Reads.txt");
	ofstream AllBadMismatchNumber_reads(dir_name + "/AllBadMismatchNumber_reads.txt");

	ofstream Expected4MM_numMM(dir_name + "/Expected4MM_numMM.txt");
	ofstream Expected5MM_numMM(dir_name + "/Expected5MM_numMM.txt");
	ofstream Expected6MM_numMM(dir_name + "/Expected6MM_numMM.txt");
	ofstream Expected7MM_numMM(dir_name + "/Expected7MM_numMM.txt");
	ofstream Expected8MM_numMM(dir_name + "/Expected8MM_numMM.txt");

	ofstream* ExpectedxMM_numMM[] = {&Expected4MM_numMM,&Expected5MM_numMM,&Expected6MM_numMM,&Expected7MM_numMM,&Expected8MM_numMM};

	if(!Recombinant_reads) {
		cout << "cannot open ofstream" << endl;
		exit(1);
	}

	cout << "Output Directory: " << dir_name << endl;
	//cout << "File: " << dir_name + "/Recombinant_Reads.txt" << endl;

	//Recombinant_reads << "hi" << endl;
	//exit(1);


	// NOTE: both SAM files must have the reads in the exact same order; when being parsed here, it treats the next alignments in both
	// 	sam files as being for the same read
	int j = 0;
	while(sam_file_1 && sam_file_2) {

		for(int s = 0; s < 2; s++) { // changed to 2! Test both SAM files in a row
			// 1. PARSE MD
			mm_locations.clear();
			*sam_files[s] >> name;
			// Ignore secondary alignments (only take first)
			while(name == prev_name && sam_files[s]) {
				getline((*sam_files[s]), line);
				*sam_files[s] >> name;
			}
			if(!sam_files[s]) {
				cout << "[DONE] Finished processing reads from at least one sam file" << endl;
				break;
			}
			
			OSS << endl << "####### Processing read " << name << " for parent s = " << s << ", parsing SAM..." << endl;
			*sam_files[s]>> num; //  bitwise flag
			*sam_files[s]>> chromosome; // which chromosome aligned to
			*sam_files[s]>> start_pos;  // position on chromosome
			*sam_files[s]>> num;  // alignment score

			*sam_files[s]>> line; // consume rest of CIGAR string

			*sam_files[s]>> line; // MRNM (for paired reads)
			*sam_files[s]>> num; //  MPOS (for paired reads)
			*sam_files[s]>> num; //  Inferred insert size

			*sam_files[s]>> seq; //  sequence of read
			*sam_files[s]>> line; // quality scores

			//extra field with flags
			*sam_files[s]>> line; // NM: number of mismatches (of any type, SNP, DEL, INS)
			*sam_files[s]>> line; // ms:
			*sam_files[s]>> line; // AS:
			*sam_files[s]>> line; // nn:
			*sam_files[s]>> line; // tp:
			*sam_files[s]>> line; // cm:
			*sam_files[s]>> line; // s1:
			*sam_files[s]>> line; // s2:
			*sam_files[s]>> MD; // de:
			while(MD.substr(0, 2) != "MD")
				*sam_files[s]>> MD;

			getline(*sam_files[s], line); //consume the rest of the optional fields

			num_MM = 0;
			pos = 0;

			MDS << MD << endl;
			for(int j = 0; j < 5; j++) // remove the "MD:Z:" characters
				MDS >> base;

			while(true) {
				if(MDS >> num) {
					pos += num;
					//cout << "read in number: " << num << endl;
				} else {
					//MDS.clear();
					MDS.clear();
					if(MDS >> base) {
						//cout << "read in base: " << base << endl;
					}
					else {
						// cout << "failed to read base, MD string complete" << endl;
						break;
					}
					if(base == '^') {
						//cout << " >> read in ^" << endl;
						while(true) {
							MDS.get(base);
							if(isdigit(base)) {
								MDS.unget();
								//cout << "     >> read in number, breaking" << endl;
								break;
							}
							pos++;
							
							//cout << "      >> read in character" << endl;
						}
					} else {
						// found a mismatch
						num_MM ++;
						pos ++;

						// mdparser2: add the start position to the position
						//   although this will increase processing in the MD parse, this number will need to be looked up at some point
						mm_locations.insert(pos + start_pos);
					}
				}
			}

			MDS.clear();
			end_pos = pos + start_pos;
			//output_file << name << '\t' << "LENGTH:" << seq.length() << '\t' << "TOTALMM:" << num_MM << '\t' << "LOCATIONS:" << OSS.str() << endl;
			//output_file << name << '\t' << chromosome << '\t' << start_pos << '\t' << end_pos << OSS.str() << endl;

			// CHECK MD

			int mismatch = 0;
			int match = 0;
			double mismatch_rate;

			//int chromosome_aligned = (*chromosome_names_directory[s])[chromosome];
			int chromosome_aligned = chromosome_names[chromosome];

			int prev_result = -1;
			int num_switches = 0;

			int expected_mismatches = 0;

			OSS << "  >>Aligned to chromosome " << chromosome_aligned  << " (called " << chromosome << ")" << ", finding mismatches..." << endl;
			if(chromosome_aligned > 0) {
				for(int b = start_pos + 1; b < end_pos; b++) {
					// iterate through the called SNPs to see if they actually occur
					if(mismatch_locations[s][chromosome_aligned - 1][b] == 1) {
						// should be a mismatch here
						expected_mismatches ++;
						if(mm_locations.find(b) != mm_locations.end()) {
							//found the mismatch in the set
							mismatch++;
							if(s == 0) OSS << "[U]\t";
							if(s == 1) OSS << "[S]\t";
							if(prev_result == 1) num_switches ++;
							prev_result = 0;
							
						} else {
							// no mismatch
							match++;
							if(s == 0) OSS << "[S]\t";
							if(s == 1) OSS << "[U]\t";
							if(prev_result == 0) num_switches ++; //switches from being a mismatch to being a match
							prev_result = 1;
						}
					}
				}

			} else {
				OSS << "Chromosome did not map to chr(1-9), skipping analysis" << endl;
			}

			for(int i = 4; i < 9; i++) {
				if(expected_mismatches == i) {
					//*ExpectedxMM_numMM[i - 4] << name;
					*ExpectedxMM_numMM[i - 4] << mm_locations.size() << endl;
				}
			}
			
			mismatch_rate = ((double)(mismatch) / (double)(match + mismatch));
			OSS << endl << "Finished processing " << end_pos - start_pos << " length read from position " << start_pos << " to " << end_pos 
				<< ": " << mismatch_rate * 100 << "% of expected mismatch locations are mismatches: ";
			bool print_info = false;
			// POSSIBLE FILTER: remove reads with too many mismatches	
			// if(expected_mismatches < mm_locations.size()) {
			// 	// filter out this read.
			// 	alignment_outcomes[s] = 5;
			// 	OSS << "Read is filtered out: more mismatches than expected:" << endl;
			// 	OSS << "   Num Expected:   " << expected_mismatches << endl;
			// 	OSS << "   Num Mismatches: " << mm_locations.size() << endl;
				
			// } else 
			if(mismatch + match < 2){
				print_info = false;
				// too little spanned sites
				alignment_outcomes[s] = 3;
				// total_Discarded++;
				//total_Recombinant++;
			} else if(mismatch_rate <= 0.0001) { // account for floating point precision
				if(s == 0){ alignment_outcomes[0] = 0; OSS << "Read is Salinas" << endl;}
				if(s == 1){ alignment_outcomes[1] = 1; OSS << "Read is US96UC23" << endl;}//countUS96UC23();
				print_info = true;
			} else if(mismatch_rate >= 0.9999) {
				if(s == 0){ alignment_outcomes[0] = 1; OSS << "Read is US96UC23" << endl;}//countUS96UC23();
				if(s == 1){ alignment_outcomes[1] = 0; OSS << "Read is Salinas" << endl;}//countSalinas();
				print_info = true;
			} else{
				// read is potential recombinant
				print_info = true;
				if(num_switches < 2) {
					OSS << "Read is RECOMBINANT" << endl;
					alignment_outcomes[s] = 2;
				} else {
					OSS << "Read is BAD recombinant" << endl;
					alignment_outcomes[s] = 4;
				}
			}
			if(print_info) {
				// print mismatch positions
				OSS << " Actual   Mismatch Locations:" << endl;
				for(auto it = mm_locations.begin(); it != mm_locations.end(); ++it) {
					OSS << '\t' << (*it - start_pos);
				}
				OSS << endl;
				OSS << " Expected Mismatch Locations:" << endl;
				for(int b = start_pos; b < end_pos; b++) {
					if(mismatch_locations[s][chromosome_aligned - 1][b] == 1) {
						OSS << '\t' << b - start_pos;
					}
				}
				OSS << endl;
			}

			if(s == 1) { // changed temporarily to s == 0
				prev_name = name;
				if(alignment_outcomes[0] == 5 || alignment_outcomes[1] == 5) {
					total_Bad_Mismatch_Number++;
					OSS << " ******* Read is discarded, too many Mismatches in one or more read" << endl << endl;
					AllBadMismatchNumber_reads << OSS.str();
				} else if(alignment_outcomes[0] == 4 || alignment_outcomes[1] == 4) {
					total_Bad_Recombinant++;
					OSS << " ******* Read is discarded, one or more alignments have a bad Recombination" << endl << endl;
					AllDiscarded_reads << OSS.str();
				} else if(alignment_outcomes[0] == 3 || alignment_outcomes[1] == 3) {
					total_Discarded++;
					OSS << " ******* Read is discarded, one or more alignments failed to span more than one mismatch sites" << endl << endl;
					AllDiscarded_reads << OSS.str();
				} else if(alignment_outcomes[0] == 0 && alignment_outcomes[1] == 0) { // changed temporarily to or
					total_Salinas++;
					OSS << " ******* Read is SALINAS" << endl << endl;
					Salinas_reads << OSS.str();
				}
				else if(alignment_outcomes[0] == 1 && alignment_outcomes[1] == 1) { // changed temporarily to or
					total_US96UC23++;
					OSS << " ******* Read is US96UC23" << endl << endl;
					US96UC23_reads << OSS.str();
				}
				else if(alignment_outcomes[0] == 2 && alignment_outcomes[1] == 2) { // changed temporarily to or
					total_Recombinant++;
					OSS << " ******* Read is  >>>RECOMBINANT " << endl << endl;
					Recombinant_reads << OSS.str();
				}
				else {
					total_Contradict++;
					OSS << " ******* Read is discarded, the alignments contradict each other: " << endl;
					OSS << "    Alignment s = 0 has value " << alignment_outcomes[0] << endl;
					OSS << "    Alignment s = 1 has value " << alignment_outcomes[1] << endl;
					OSS << endl;
					AllDiscarded_reads << OSS.str();
				}
			}
			
		}
		cout << OSS.str();
		OSS.clear();
		OSS.str("");
		
	}
	int total_reads = total_Salinas + total_US96UC23 + total_Recombinant + total_Bad_Recombinant + total_Discarded + total_Contradict + total_Bad_Mismatch_Number;
	int total_good_reads = total_Salinas + total_US96UC23 + total_Recombinant;
	double attrition_rate = double(total_reads - total_good_reads) / double(total_reads);
	double Recombination_Rate = (double)total_Recombinant / (double)total_good_reads;

	cout << "Finished processing all non-duplicate reads..." << endl;
	cout << "Total Salinas:     " << total_Salinas << endl;
	cout << "Total US96UC23:    " << total_US96UC23 << endl;
	cout << "Total Recombinant: " << total_Recombinant << endl;
	cout << "Observed Recombination Rate: " << Recombination_Rate * 100 << "%" << endl << endl;
	cout << "Total Bad Mismatches (too many mismatches)                                               : " << total_Bad_Mismatch_Number << endl << endl;
	cout << "Total Bad Recombinant (at least one alignment is 'Recombinant' with more than 1 switch)  : " << total_Bad_Recombinant << endl;
	cout << "Total Bad Alignment (at least one alignment <2 mismatch sites spanned or aligned to Chr0): " << total_Discarded << endl;
	cout << "Total Contradictions (alignment conclusions are not in agreement)                        : " << total_Contradict << endl << endl;
	
	cout << "Processed [" << total_reads << "] reads, " <<  attrition_rate * 100 << " % were thrown out" << endl;
	auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(end - begin).count();
    cout << endl << "Elapsed time: " << (double)elapsed / 1000000 << " seconds" << endl;

}

