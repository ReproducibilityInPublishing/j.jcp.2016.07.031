# A fast accurate approximation method with multigrid solver for two-dimensional fractional sub-diffusion equation

"A fast accurate approximation method with multigrid solver for two-dimensional
fractional sub-diffusion equation" implements an approximate solution method
for the fractional sub-diffusion equation. At each solution step, a large
dimensional array has to be inverted to get the next step. For each inversion,
a multigrid method is used to reduce the dimensionality of the matrix in
question. Then an approximate inversion equation is applied when the
dimensionality becomes small enough. After inversion, the solution is
interpolated back to the full dimensionality. 

## Build Instructions

### Requirements
Instructions were tested using Docker version 18.06.0-ce, build 0ffa825, on Ubuntu 16.04.5 LTS.

### Building with Docker
    docker build -t ${DOCKER_IMAGE_NAME} .

## Run Instructions

### Running with Docker
To start a container for the Docker image:

    docker run -it --rm -v $(pwd):/Scratch ${DOCKER_IMAGE_NAME}

#### Run Everything
Within the Docker container, to run everything, computational scripts for
experiments and visualization scripts, run

    ./run.sh

Please be aware of computational efforts for the scripts, see...

See sections below provide for details about the individual steps.

### Building Code
Within the Docker container, run

    ./build.sh

After building, the resulting artifacts are binaries:
    example1
    example1_BDADI
    example2
    example2_BFSMGM
    example3

#### Running Computational Scripts
Within the Docker container, run

    ./computation.sh

Output will be tables:
    table1.csv
    table2.csv
    table3.csv
    table4.csv

Expected tables are in directory `expected_tables/`.

## Reproduction Notes
We kept track of our progress and issues inside `notes.txt`. We also have an
jupyter notebook showing this progress over time `ReproducibilityPlot.ipynb`.

## Acknowledgements
We want acknowledge the authors for their fine work on this experiment. We
succeeded with this project where many others had failed. The authors should be
commended on putting together high quality work.
