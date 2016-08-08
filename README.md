# denovo_hit

Usage: ./denovo_hit <VCF file> <PED file>

Output: ./potential_denovo.txt

Compile: g++ -std=c++11 denovo_hit.cpp -o denovo_hit

Determines if variants calls are novel for a child given pedigree information. Program can handle both annotated and unannotated VCF files, however output file is not compatible with GATK's snpEff, so you should annotate prior to using denovo_hit.

Please see files for example VCF file, and pedigree file.

Format VCF file (tab delimited):

CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE_ID, ...

Format pedigree file (tab delimited):

SAMPLE_ID_CHILD, SAMPLE_ID_PARENT1, SAMPLE_ID_PARENT2

Output file format (tab delimited):

CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, PEDIGREE, CALL_STATS, COUNTS ...
