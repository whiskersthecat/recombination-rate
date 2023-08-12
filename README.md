# recombination-rate
A tool for marking aligned F2 PacBio or Nanopore sequencing reads as parental A, parental B, or recombinant

**checkSNPs.cc** is compiled on the command line using:
> g++ checkSNPs.cc -o checkSNPs -std=c++11

The executable **checkSNPs** takes the following command line arguments:
> ./checkSNPs _snp_file.bed_ _samfile1.SAM_ _samfile2.SAM_

**snp_file.bed** represents the variants between parent A and parent B and consists of tab-seperated lines as follows:
> _chromosome_ _SNPstartlocation_ _SNPendlocation_ _referencebase_ _alternatebase_
which represent that any sequence from parent B aligning to parent A should have a mismatch on chromosome _chromosome_ at position _SNPstartlocation_ changing base _referencebase_ to _alternatebase_.

**samfile1.SAM** is the output file of the _minimap2_ aligner aligned to parental assembly A using the options:
> ./minimap2 --MD -ax map-hifi --secondary=no _parentAreference_ _reads_ > _samfile1.SAM_
according to the specifications [minimap2](https://lh3.github.io/minimap2/minimap2.html).
