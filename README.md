What is a genome?
Genome: set of DNA
Difference between Genomics and genetics
Genetics: deals with heriditory 
Genomics deals with the organisms gimones in entirety  
### DNA: Deoxyribonucleic acid  
DNA contains:  
- A : Adenine  
- T : Thymine  
- C : Cytosine  
- G : guanine  

pairs:  
- A <-> T  
- C <-> G  

### RNA: Ribonucleic acid
- A : Adenylate  
- U : Uridylate  
- C : Cytidylate  
- G : Guanylate  

pairs:  
- A <-> U  
- C <-> G 

### DNA pairing
Paring between the top and bottom layeris inverted

5' ATCGGCTA 3'
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| | | | | | | | | |
3' TAGCCGAT 5'

3' and 5' represent the numner pf sugar  backbones.
5' -> Phosphate group
3' -> hydroxly group

#### GC content
the persentage that amots the amount of C and G in a given DNA
formula: 
(No of C + No of G)/(len of DNA)

### Codons  
They are the codes for amino acids. Pairs of 3 taken in the genome  
https://en.wikipedia.org/wiki/DNA_codon_table
![alt text](
Codons table
```
'M' - START, '_' - STOP
"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
"TGT": "C", "TGC": "C",
"GAT": "D", "GAC": "D",
"GAA": "E", "GAG": "E",
"TTT": "F", "TTC": "F",
"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
"CAT": "H", "CAC": "H",
"ATA": "I", "ATT": "I", "ATC": "I",
"AAA": "K", "AAG": "K",
"TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
"ATG": "M",
"AAT": "N", "AAC": "N",
"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
"CAA": "Q", "CAG": "Q",
"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
"TGG": "W",
"TAT": "Y", "TAC": "Y",
"TAA": "_", "TAG": "_", "TGA": "_"
```

### open reading frames (ORF)
These frames are read by ribosomes.
open reading frame is a portion of DNA molecule that when translated to an amino acid contains no "stop"codon.

Sturcture of open Reading frame
START - Codon - Stop
Codons give mRNA
mRNA gives functional protein

for eg:-
given DNA sequence:
AGTATAAGGTGCTGGCACATCAGTATGGCTTACTGAGCGTGGGCGCGAAG

output for amino acids for entire seq:
['S', 'I', 'R', 'C', 'W', 'H', 'I', 'S', 'M', 'A', 'Y', '_', 'A', 'W', 'A', 'R']

output ORF:
['S', 'I', 'R', 'C', 'W', 'H', 'I', 'S', 'M', 'A', 'Y', '_', 'A', 'W', 'A', 'R']
['V', '_', 'G', 'A', 'G', 'T', 'S', 'V', 'W', 'L', 'T', 'E', 'R', 'G', 'R', 'E']
['Y', 'K', 'V', 'L', 'A', 'H', 'Q', 'Y', 'G', 'L', 'L', 'S', 'V', 'G', 'A', 'K']
['L', 'R', 'A', 'H', 'A', 'Q', '_', 'A', 'I', 'L', 'M', 'C', 'Q', 'H', 'L', 'I']
['F', 'A', 'P', 'T', 'L', 'S', 'K', 'P', 'Y', '_', 'C', 'A', 'S', 'T', 'L', 'Y']
['S', 'R', 'P', 'R', 'S', 'V', 'S', 'H', 'T', 'D', 'V', 'P', 'A', 'P', 'Y', 'T']

First row continue with the 3 letter set with the numbering starting from the first letter onwards 

the 'S' in the first row is the amino acid corresponding to the 3 letters set starting form the second letter in the DNA sequence "AGT"

Second row continue with the 3 letter set with the numbering starting from the Seond letter onwards

the 'V' in the first row second letter is the amino acid corresponding to the3 letters set starting form the second letter in the DNA sequence "GTA" 
