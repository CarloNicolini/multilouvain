# multilouvain
Matlab/C++ implementation of Traag Louvain methods appliable to different quality functions. This code has some modifications with respect to the Traag's one to make available calls to member function that have been moved from private to public, as well as other minor things.
You can find the original Vincent Traag code here:

    https://github.com/vtraag/louvain-igraph

To clone this repository you must

    git clone https://github.com/carlonicolini/multilouvain
    git submodule update --init --recursive

The last command is because it contains Eigen as a submodule

To compile the code:

    cd multilouvain
    mkdir build
    cmake -DMATLAB_SUPPORT=True ..
    make

This will make the `multilouvain.mex64` mex file that you can use as wrapper to run the Louvain algorithm on square adjacency matrices.

You can also compile the code for Octave, but first because of some incompatibilities between Octave and Matlab mex files, you have to clean the repo:

    git clean -df .
    mkdir build
    cmake -DOCTAVE_SUPPORT=True ..
    make

This will make the `multilouvain.mex` file to be used within Octave. To compile it you need the following packages

    sudo apt-get install liboctave-dev

# Compiling with no admin privileges

You don't need to install igraph to compile multilouvain. It is possible to compile igraph and install it locally, as shown in these lines:

1. Download igraph and compile it
    
    cd
    wget igraph-0.7.1.tar.gz
    tar zxvf igraph-0.7.1.tar.gz
    cd igraph-0.7.1
    ./configure --prefix=$HOME/igraph-0.7.1
    make install

Now go to your `multilouvain` folder and run cmake with the `IGRAPH_INCLUDES` and `IGRAPH_LIBRARIES` options, set to the `include` and `lib` folders of your local igraph installation.

    cd multilouvain
    mkdir build
    cd build
    cmake -DIGRAPH_INCLUDES=$HOME/igraph-0.7.1/include -DIGRAPH_LIBRARIES=$HOME/igraph-0.7.1/igraph/lib ..
    make
