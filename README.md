# README
_This repository contains the files related to the publication (NOMBRE PUBLICACION)_
## To replicate the results

To replicate the results obtained, the following steps must be carried out.

1. MUMer 3.0 must be installed.
2. For each sequence, the following commands must be run in a Linux enviroment and from the directory containing all sequences (reference and strains):

`` nucmer --prefix=ref_qry ref.fasta qry.fasta ``  
`` show-snps -CHTlr -x 10 ref_qry.delta > ref_qry.snps ``

3. Then all the information is integrated into one file (``example_file.txt``) using the following commands:

``ls -1 *.snps > my_list``  
``foreach i ( `cat my_list` )``  
``cat $i | awk -v awkvar="$i" '{printf "%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s\n", $1, 9, $2, 9, $3, 9, $4, 9, $5, 9, $6, 9, $7, 9, $8, 9, $9, 9, $10, 9, $11, 9, $12, 9, $13, 9, $14, 9, awkvar}' > $i.out``  
``cat *.out >> example_file.txt``

4. ``visualization.py``, ``table.py`` and ``example_file.txt`` must be in the same directory. Run ``visualization.py`` to generate two files. The first one contains, for all the sequences, the nucleotides that are in the positions where a SNP was found. The second file contains the positions where the SNPs were found with respect to the reference genome.

## Considerations

The file ``visualization.py`` can be edited to change the names of input and output files. On line 3 the parameter ``context`` can be set to ``True`` if it is desired that the SNP file contains a 10-nucleotide upstream and downstream context for each SNP.

If a strain genome has an insertion, a gap at that position will be displayed for the rest of the sequences (reference included) that do not have that SNP.
