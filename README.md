# recombination-rate
A tool for classifying aligned F2 PacBio or Nanopore sequencing reads as parental A, parental B, or recombinant.

## checkSNPS.cc
### Usage
**checkSNPs.cc** is compiled on the command line using:
> g++ checkSNPs.cc -o checkSNPs -std=c++11

The executable **checkSNPs** takes the following command line arguments:
> ./checkSNPs _snp_file.bed_ _samfile1.SAM_ _samfile2.SAM_

### _SNP File_
**snp_file.bed** represents the variants between parent A and parent B and consists of tab-seperated lines as follows:
> _chromosome_ _SNPstartlocation_ _SNPendlocation_ _referencebase_ _alternatebase_

which represent that any sequence from parent B aligning to parent A should have a mismatch on chromosome _chromosome_ at position _SNPstartlocation_ changing base _referencebase_ to _alternatebase_.

### _SAM File_
**samfile1.SAM** is the output file of the _minimap2_ aligner aligned to parental assembly A using the options:
> ./minimap2 --MD -ax map-hifi --secondary=no _parentAreference.fasta_  _reads.fasta_ > _samfile1.SAM_

according to the specifications [minimap2](https://lh3.github.io/minimap2/minimap2.html).

## Algorithm
For each read in **samfile1.SAM**, check which mismatches should occur within the length of the aligned read according to **snp_file.bed** iff the read was of parent B origin. Call these expected mismatches.

Parse the MD:Z tag of the SAM file to determine the position of mismatches in the alignment. 

For each expected mismatch location, check if the read actually contains the mismatch. If it does, mark that SNP as of parent B origin. Otherwise, mark it as parent A origin.  

Reads and subsequently categorized as follows:  
1. Discard any reads with more mismatches than expected mismatches. _These are presumed to be sequencing errors and this read can be ignored._  
2. Discard any remaining reads with less than two expected mismatches. _For recombination analysis, these will provide no data. Note that this includes reads that fail to align to any chromosome in the assembly (i.e. chr0)_  
3. If no mismatch sites have mismatches, classify the read as **parent A**.  
4. If all mismatch sites have mismatches, classify the read as **parent B**.  
5. If the read switches between match and mismatch more than once, discard it. _Any read with this behavior is incorrectly marked, as there cannot be more than one recombination_  
6. If the read switches between match and mismatch exactly once, mark it as **recombinant**.  

## Output
In addition to logging a summary to **cout**, checkSNPs generates the following output files in a new folder in the directory of the SAM file. Each contains detailed information about the processing of a group of reads.

- **US96UC23_Reads.txt** Reads marked according to 4 above.
- **Salinas_Reads.txt**  Reads marked according to 3 above.
- **Recombinant_Reads.txt** Reads marked according to 6 above.
- **AllDiscarded_Reads.txt** Reads marked according to 1 and 2 above.
- **AllBadMismatchNumber_reads.txt** Reads marked according to 5 above.

## Benchmarking
In order to ensure that **checkSNPs** categorizes reads, feed it reads that originate from only one parent.
Early tests prove very poor performance:
> When processing 99989 parent B reads, 90123 are identified as parent B, but 33 are incorrectly marked as recombinant. Recombinant reads have mistake mismatches that make some mismatch sites appear to be from parent A.
In order to properly assess recombination rate, the error rate must be vastly reduced. It is estimated that only 3 among 50000 reads should be recombinant.
