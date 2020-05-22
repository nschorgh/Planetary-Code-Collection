#!/bin/bash

# run mars_mapi2p


# name of executable
PROG=./mars_mapi2p
PROG=./a.out


# run up to 4 jobs in parallel with integers 1...N as command line arguments,
# where N is the number of lines in input file mapgrid.slp
parallel --jobs 4  $PROG mapgrid.slp ::: $(seq 1 2)
# (in this example N is only 2, but it could be much larger)



# On a cluster use a job scheduler to provide the input index.
