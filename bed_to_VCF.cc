// Convert Alex's .bed file to a variant call file (VCF)
// specifications: http://samtools.github.io/hts-specs/VCFv4.4.pdf


// chromosome mismatch_location1 mismatch_location2 base1 base2

#include<iostream>
#include<fstream>
#include<string>
#include<map>
using namespace std;
string chrom;
int loc1;
int loc2;
char base1;
char base2;

map<string, string> chromosome_names_1 = {
{"Lsat_1_Salinas_v11_chr9", "Lsat_1_v11_chr9"},
{"Lsat_1_Salinas_v11_chr8", "Lsat_1_v11_chr8"},
{"Lsat_1_Salinas_v11_chr7", "Lsat_1_v11_chr7"},
{"Lsat_1_Salinas_v11_chr6", "Lsat_1_v11_chr6"},
{"Lsat_1_Salinas_v11_chr5", "Lsat_1_v11_chr5"},
{"Lsat_1_Salinas_v11_chr4", "Lsat_1_v11_chr4"},
{"Lsat_1_Salinas_v11_chr3", "Lsat_1_v11_chr3"},
{"Lsat_1_Salinas_v11_chr2", "Lsat_1_v11_chr2"},
{"Lsat_1_Salinas_v11_chr1", "Lsat_1_v11_chr1"}
};

int main(int argc, char* argv[]) {
	if(argc < 2) {
		cout << "Usage: ./bed_to_VCF snp_file.bed" << endl;
		exit(1);
	}
	string filename(argv[1]);
	ifstream bed(filename);
	ofstream vcf(filename + ".vcf");
	//ifstream bed("./SNP/SNP-table_Lsat_1_Salinas_Genome_v11.01.GB.SNPs_US96UC23.CLC.bed");
	//ofstream vcf("./SNP/SNP-table_Lsat_1_Salinas_Genome_v11.01.GB.SNPs_US96UC23.vcf");
	vcf << "##fileformat=VCFv4.2" << endl;
	vcf << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
	int i = 0;
	while(bed) {
		bed >> chrom; if(!bed) break; 
		bed >> loc1;
		bed >> loc2;
		bed >> base1;
		bed >> base2;

		if(chromosome_names_1[chrom] == "") break;
		if(true) {
			vcf << chromosome_names_1[chrom] << '\t';
				  vcf << loc2 << "\t.\t";
			
			vcf << base1 << '\t';
			vcf << base2 << "\t.\tPASS\t." << endl;
			//if(chromosome_names_1[chrom] == "Lsat_1_v11_chr2")
			//	break;
		}
		i++;
		if(i % 10000 == 0)
			cout << i << endl;
	}

	return 0;
}