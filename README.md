# multilouvain
Matlab/C++ implementation of Traag Louvain methods appliable to different quality functions. This code has some modifications with respect to the Traag's one to make available calls to member function that have been moved from private to public, as well as other minor things.

To compile the code:

    mkdir build
    cmake -DMATLAB_SUPPORT=True ..
    make

This will make the multilouvain.mex64 mex file that you can use as wrapper to run the Louvain algorithm on square adjacency matrices.
