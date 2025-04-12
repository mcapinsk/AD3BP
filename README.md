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

For detailed description on how to build the library see

http://capd.ii.uj.edu.pl/html/capd_compilation.html

## Building and running the Arnold diffusion proof

To compile and run the full proof of Arnold diffusion execute:

    cd ../AD3BP
    make
    ./AD3BP

Instead of the above,

    cd ../AD3BP
    make
    ./AD3BP 12345

will validate just a single connecting sequence, which should take under 30 seconds. Instead of 12345 one can choose any number between 0 and 89999.

The full proof will take roughly 30 days on a single thread. The computation runs on all threads available on the machine. (It is a good idea to run the proof on some unused office desktop computer or on a cluster.) The results for successive computations of the proof will be stored in result files, as outlined below.

## Results

During the computation the results are stored in the folder:

AD3BP/AD3BP/results

Once the proof is conducted successfully the file:

AD3BP/AD3BP/results/0_final_result.txt

will store the final result. Should the proof fail due to errors or due to required bounds not being satisfied, this will be reported in files

AD3BP/AD3BP/results/errors.txt
AD3BP/AD3BP/results/failure_.txt

(After a successfull run these files will be empty.)

## Suggested order of reading the files

- VectorField.h This file contains routines for reading the vector field into the CAPD computer assisted proof. The formulae themselves are in the folder AD3BP/1_VectorField. They were written into the text files by using a simple Mathematica routine in 
AD3BP/1_VectorField/MathematicaFormulaeComputation.nb

- readSections.h/cpp 
This reads the sequence of points q[i] at which the sections are positioned, the sequence of matrices A[i] which define the sections, and the matrices B[i] which give the local coordinates on the sections. These are stored in the folder:
AD3BP/0_sections

- localCoordinateChange.h/cpp 
This file contains classes for passing to and from loacl coordinates on the surfaces of sections. These changes of coordinates are explained in section "5 Surfaces of sections and local coordinates" from the paper.

- localMap.h/cpp
This contains a class for the computation of the section-to-section map along the flow, expressed in the appropriate local coordinates.

- main.cpp
Here we validate all the connecting sequences and compute the global bound (66) on the change of I from the paper. 

- constants.h
This file contains the mass parameter "mu" of the Neptune-Triton system so that it can be accessed by various parts of the program. 

## Authorship

Only the contents of the folder AD3BP/AD3BP constitute the proof of Arnold diffusion in the 3bp, and the files included there have been written by Maciej J. Capinski. The remaining files are part of the CAPD library. These have been created by the CAPD Group.


