NCBI taxonomy tree made easy
=============================

This is a simple program that I use to query the NCBI taxonomy
tree. It requires the ETE python library (ete.cgenomics.org) to work.
Features are still very rudimentary. Please, refer/cite this
repository if you use the program.

Usage:
*********

update taxonomy tree (download and parse the latest NCBI taxonomy DB): 
-----------------------------------------------------------------------
  $ wget  ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

  $ tar zxf taxdump.tar.gz 

  $ python update_taxadb.py

get help:
------------
  $ python ./ncbi_taxa.py -h 

get name-to-taxid translation: 
------------------------------------
  $ python ./ncbi_taxa.py -n Bos taurus, Gallus gallus, Homo sapiens 

get NCBI topology from species names:
------------------------------------------------
  $ python ./ncbi_taxa.py -n Bos taurus, Gallus gallus, Homo sapiens -x

get NCBI lineage and info from species names: 
------------------------------------------------
  $ python ./ncbi_taxa.py -n Bos taurus, Gallus gallus, Homo sapiens -i

get the same as above using taxids: 
------------------------------------
  $ python ./ncbi_taxa.py -t 9913 9031 9606 -x

Future features: 
******************

* Support for incomplete taxa names


Contact: jhcepas[at]gmail.com
