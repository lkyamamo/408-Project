preface:
some of the code written here was done by chatGPT but all classes and respective methods in addition to how these classes initeract with the external libraries were done by myself. In instances where code or ideas were adapted from other sources, they are cited (most notably the mpi source usage in the DockerFile). The current state of this library is to showcase the multiple methods that can be used for numerical solutions including known discrete definitions, method of lines from py-pde, and finite element method using FEniCS. FEniCS was not used in the end, but is still incorporated into class methods because this is an ongoing project. 

setup:

In order for this library and example files to work you need to have docker installed on your computer. 
Instructions for this can be found here https://www.docker.com/get-started/. Come back after you have successfully installed it and can run "docker --version"

1. get your terminal into the same directory as this file and run the following line
    "docker build -t {username}/408-project:latest ." 
    This will build your environment to include the necessary external libraries
    NOTE: this build does take a substantial amount of time (~15 min)

2. run 
    "docker images"
    and save the image id of the image that has the tag you specified before

3. run 
    "docker run -it {image_id}"
    This will take you into the docker container and allow you to run the scripts from the command line

There are example files that you can run outright that will export different gif's of potential scalar fields.