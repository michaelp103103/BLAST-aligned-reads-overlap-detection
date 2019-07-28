# BLAST_overlap_finder
A Cythonized Python program that identifies overlaps between DNA reads based on the BLAST alignment of the reads to their reference 
genomes. 

# General Usage
The Python code for the overlap finder must first be Cythonized using the following command:

```
python setup.py build_ext --inplace
```

Once the Python code is Cythonized, the program is ready to be run. The program accepts three arguments. 
`filepath1` should be set to the file path to the BLAST output file and `savepath1` to the desired file path to 
which the program output will be saved. The minimum overlap length (`overlap_length1`) is set to 50 bp by default and may be changed to 
any minimum overlap length greater than or equal to 1 bp. 

The following is a template that may be used for running the program. The template is also located in the file `find_overlaps_cy.py`.

```python
import find_overlaps_cy

find_overlaps_cy.find_overlaps(filepath1='/file/path/to/blast/output', overlap_length1 = 50, 
savepath1='/save/file/path')
```

The program recognizes BLAST output files of the following output format only:

```
-outfmt "6 qseqid sstart send evalue bitscore pident"
```

The program outputs a tab-delimited text file in the following format:

```
read_id1  read_id2  read_length1  read_length2  overlap_length  strand_orientation1  strand_orientation2
```

Where `read_id1` and `read_id2` are the identification numbers of two overlapping reads. `read_length1` and `read_length2` are the 
lengths of the two overlapping reads corresponding to `read_id1` and `read_id2`, respectively. `overlap_length` is the length of the 
overlap. `strand_orientation1` and `strand_orientation2` are the orientations of the reads corresponding to `read_id1` and `read_id2`, 
respectively. Strand orientation can be either `F`, for forward, or `R`, for reverse.

The output file does NOT contain any duplicate overlaps (i.e., if A overlaps with B, then this overlap will be recorded only once) and 
it does NOT contain record a read overlapping with itself.

# Authors
Michael Piskozub

# Acknowledgements
Thank you to Dr. Jaroslaw Zola and Vicky Zheng of SCoRe group at the University at Buffalo, as well as CLIMB UP of the University at 
Buffalo. 
