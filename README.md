# PPCAS-MEL
Memory Efficient Probabilistic Pairwise model for Consistency-based multiple alignment in Apache Spark - SUPE2018

Prerequisites
--------------
PPCAS compilation requires the following tools installed on your system ``make``, ``gcc-c++`` and ``Glib 2.14 or superior``. 

The execution requires a ``Hadoop`` and ``Spark`` infrastructure with the environment variables correctly set and its ``path``. Also a ``Python`` installation with ``Numpy`` is needed.


Compile 
--------
Download the git repository on your computer.
    
Make sure you have installed the required dependencies listed above. 
When done, move in the project root folder and enter the following commands:     
    
    $ make
    

The shared library will be automatically generated.


Usage
--------
It is included a script named ``run`` which executes PPCAS with the required and optional parameters.

Required parameters:

    $ -i [sequence_file]
    $ -m [master_spark_ip]
    
Optional parameters:

    $ -b (library bound, default: none)
    $ -h (HDFS path, default: /user/root)
    $ -o (output, collection of files)
    $ -p (number of partitions, default: number of sequences)
    $ -t (T-Coffee output, single file)
    

Example
--------

There are input sequences in the examples folder.

``BB11001.tfa`` a small dataset from ``BAliBASE``.

``ghf13_*`` being  ``*`` the number of sequences obtained from ``HomFam`` dataset.

Calculate the library with a maximum of 500 constraints (BB11001.tfa contains 1540 constraints):

    $ ./run.sh -i examples/BB11001.tfa -m 192.168.101.51 -b 500
    
