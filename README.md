# A fast accurate approximation method with multigrid solver for two-dimensional fractional sub-diffusion equation

The authors of [A fast accurate approximation method with multigrid solver for two-dimensional fractional sub-diffusion equation](https://dx.doi.org/10.1016/j.jcp.2016.07.031) Xue-lei Lin, Xin Lu, Micheal K. Ng, and Hai-Wei Sun were kind enough to share the code they used for their computational experiments with us. If you use this code please be sure to cite their work. You can use the bibtex citation in `cite.bib` or you can use the citation below:

Xue-lei Lin, Xin Lu, Micheal K. Ng, and Hai-Wei Sun "A fast accurate approximation method with multigrid solver for two-dimensional fractional sub-diffusion equation", Journal of Computational Physics, 323, 204-218, 2016

## Docker

To build the docker image, go into the root directory and run:

docker build -t multigrid:fresh .   # This makes the image with the name multigrid, with version fresh

To run the experiment scripts in the docker image, run:

docker run -it -v $(pwd):/Scratch multigrid:fresh

This runs the same image (multigrid:fresh) from before, while mapping current
directory to /Scratch inside the spun up docker container. The spun up
container starts you off in the bash shell. While inside the docker container,
just run:

bash run_experiments.sh

This will build all the C/C++ code and then run the different scripts
(table*.sh) that creates the tables within corresponding table*.csv files. The
run_experiments.sh script also runs the check.sh that checks the created
table*.csv files have error values that match the expected values.
