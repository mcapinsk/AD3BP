# AD3BP
The proof of Arnold diffusion in the full three body problem. The proof is based on the CAPD library:

CAPD main webpage: http://capd.ii.uj.edu.pl

Full documentation: http://capd.ii.uj.edu.pl/html/

CAPD requirements: http://capd.ii.uj.edu.pl/html/capd_requirements.html

The proof also requires the OpenMP:

https://www.openmp.org

Please make sure you have it installed.

## Quick guide on how to build the CAPD library required for the proof

Clone the repository:

    git clone https://github.com/mcapinsk/AD3BP.git
    
Enter the repository, create the build folder, configure the library and then build:

    cd AD3BP
    mkdir build
    cd build
    cmake ..
    make

The above commands will build only the CAPD library.

For detailed decription on how to build the library see

http://capd.ii.uj.edu.pl/html/capd_compilation.html

## Building and running the Adnold diffusion proof

To compile and run the proof of Arnold diffusion execute:

    cd ../AD3BP
    make
    ./AD3BP

The proof will take roughly 70 days on a single thread. The computation runs on all threads available on the machine. (It is a good idea to run the proof on some unused office desktop computer or on a cluster.) The results for succesive computations of the proof will be stored in result files, as outlined below.

## Results

During the computation the results are stored in the folder:

AD3BP/AD3BP/results

Once the proof is conducted successfuly the file:

AD3BP/AD3BP/results/0_final_result.txt

will store the final result. Should the proof fail due to errors or due to required bounds not being satisfied, this will be reported in files

AD3BP/AD3BP/results/errors.txt
AD3BP/AD3BP/results/failure_.txt

(After a succesful run these files should be empty.)

## Authorship

Only the contents of the folder AD3BP/AD3BP constitute the proof Arnold diffusion in the 3bp, and the files included there have been written by Maciej J. Capinski. The remaining files are part of the CAPD library. These have been created by the CAPD Group.


