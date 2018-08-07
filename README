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
