kraken-biom
===========

Create BIOM-format tables (http://biom-format.org) from Kraken output 
(http://ccb.jhu.edu/software/kraken/).

Installation
------------

From PyPI:

.. code-block:: bash

    $ pip install kraken-biom

From GitHub:

.. code-block:: bash

    $ pip install git+http://github.com/smdabdoub/kraken-biom.git

From source:

.. code-block:: bash

    $ python setup.py install

From docker:

.. code-block:: bash

    $ git clone https://github.com/smdabdoub/kraken-biom.git && cd kraken-biom
    $ docker build . -t kraken_biom
    $ docker run -it --rm -v ${pwd}:/data kraken_biom


Citation
--------
kraken-biom does not yet have a published article, but it can be cited as:

    Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software].  
    Available at https://github.com/smdabdoub/kraken-biom.

Requirements
------------

- biom-format >= 2.1.5

Documentation
-------------

The program takes as input, one or more files output from the kraken-report
tool. Each file is parsed and the counts for each OTU (operational taxonomic
unit) are recorded, along with database ID (e.g. NCBI), and lineage. The
extracted data are then stored in a BIOM table where each count is linked
to the Sample and OTU it belongs to. Sample IDs are extracted from the input
filenames (everything up to the '.').

OTUs are defined by the --max and --min arguments. By default these are
set to Order and Species respectively. This means that counts assigned
directly to an Order, Family, or Genus are recorded under the associated
OTU ID, and counts assigned at or below the Species level are assigned to
the OTU ID for the species. Setting a minimum rank below Species is not yet
available.

The BIOM format currently has two major versions. Version 1.0 uses the 
JSON (JavaScript Object Notation) format as a base. Version 2.x uses the
HDF5 (Hierarchical Data Format v5) as a base. The output format can be
specified with the --fmt option. Note that a tab-separated (tsv) output
format is also available. The resulting file will not contain most of the
metadata, but can be opened by spreadsheet programs.

Version 2 of the BIOM format is used by default for output, but requires the
Python library 'h5py'. If the library is not installed, kraken-biom will 
automatically switch to using version 1.0. Note that the output can 
optionally be compressed with gzip (--gzip) for version 1.0 and TSV files. 
Version 2 files are automatically compressed.

Currently the taxonomy for each OTU ID is stored as row metadata in the BIOM
table using the standard seven-level QIIME format: k__K; p__P; ... s__S. If
you would like another format supported, please file an issue or send a pull
request (note the contribution guidelines).
::

    usage: kraken-biom [-h] [--max {D,P,C,O,F,G,S}] [--min {D,P,C,O,F,G,S}]
                          [-o OUTPUT_FP] [--fmt {hdf5,json,tsv}] [--gzip]
                          [--version] [-v]
                          kraken_reports [kraken_reports ...]

Usage examples
--------------

1. Basic usage with default parameters::

    $ kraken-biom S1.txt S2.txt

  This produces a compressed BIOM 2.1 file: table.biom

2. BIOM v1.0 output::

    $ kraken-biom S1.txt S2.txt --fmt json

  Produces a BIOM 1.0 file: table.biom

3. Compressed TSV output::

    $ kraken-biom S1.txt S2.txt --fmt tsv --gzip -o table.tsv

  Produces a TSV file: table.tsv.gz

4. Change the max and min OTU levels to Class and Genus::

    $ kraken-biom S1.txt S2.txt --max C --min G

5. Basic usage with default parameters and metadata::

    $ kraken-biom S1.txt S2.txt -m metadata.tsv
This produces a compressed BIOM 2.1 file: table.biom

Program arguments
-----------------

positional arguments::

    kraken_reports        Results files from the kraken-report tool.

optional arguments::
    
      -h, --help            show this help message and exit
      --max {D,P,C,O,F,G,S}
                            Assigned reads will be recorded only if they are at or
                            below max rank. Default: O.
      --min {D,P,C,O,F,G,S}
                            Reads assigned at and below min rank will be recorded
                            as being assigned to the min rank level. Default: S.
      -o OUTPUT_FP, --output_fp OUTPUT_FP
                            Path to the BIOM-format file. By default, the table
                            will be in the HDF5 BIOM 2.x format. Users can output
                            to a different format using the --fmt option. The
                            output can also be gzipped using the --gzip option.
                            Default path is: ./table.biom
     -m METADATA, --metadata METADATA
                            Path to the sample metadata file. This should be in
                            TSV format. The first column should be the sample id.
                            This is the same name as the input files. If no
                            metadata is given, basic metadata is added to help
                            when importing the biom file on sites like phinch
                            (http://phinch.org/index.html).

      --fmt {hdf5,json,tsv}
                            Set the output format of the BIOM table. Default is
                            HDF5.
      --gzip                Compress the output BIOM table with gzip. HDF5 BIOM
                            (v2.x) files are internally compressed by default, so
                            this option is not needed when specifying --fmt hdf5.
      --version             show program's version number and exit
      -v, --verbose         Prints status messages during program execution.
