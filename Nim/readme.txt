This directory contains the later and much faster 'nim' implementation of the algorithm described in my thesis.
Specifically, the program in 'circsearch.nim' computes all surfaces with Delta up to a given bound.
The file 'orbit.nim' contains code for finding computing the set of planes in an orbit intersecting the fundamental domain.
The file 'eisenstein.nim' implements the eisenstein integers with operator overloading for convenience, along with matrices.

If you're using windows, you can run the included executables without installing anything:

circsearch_no_orientation.exe:
    Computes a list of surfaces without orientation.
    usage:
        circsearch_no_orientation.exe bound
    The argument bound determines the maximum magnitude of Delta searched up to, i.e. 0 > Delta > -bound
    Output is 'outfile.csv', which you can open with Microsoft Excel or equivalent

circsearch_orientation.exe:
    Computes a list of surfaces without orientation. Slower than circsearch_no_orientation.exe.
    usage:
        circsearch_orientation.exe bound
    The argument bound determines the maximum magnitude of Delta searched up to, i.e. 0 > Delta > -bound
    Output is 'outfile.csv', which you can open with Microsoft Excel or equivalent

Try running 'circsearch_no_orientation.exe 100' from command prompt

If you are on another system and you have nim installed, then you can build these executables like so:

    nim c circsearch.nim

    OR

    nim c -d:compute_orientation circsearch.nim

The output in both cases will be called circsearch.exe
