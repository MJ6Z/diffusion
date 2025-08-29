# Nuclear Design Diffusion Model
A part of my A-level computer sciene NEA project.

The idea of the program is to simulate the core of a nuclear reactor with a low level (but not no) realism using differential equations. This repo was orginally created as a source for code that simulates divergence-gradient diffision of one variable over a hexagonal grid.

The program allows the user to change a variety of factors in [params.json](/params.json) aswell assign certain hexagons special properties. These being nuclear fuel rods, control rods, or coolant channels as well as neutron source locations.

The model is heavily dependant on the [mathplot](https://github.com/sebsjames/mathplot) library installed in-tree with the code.


## dependancies
This project requires an installation of Mathplot's dependancies.

Extract from [Mathplot/README.md]()
```bash
# Install dependencies for building graph1.cpp and (almost) all the other examples (assuming Debian-like OS)
sudo apt install build-essential cmake git wget \
                 nlohmann-json3-dev librapidxml-dev \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev \
                 libglfw3-dev libfreetype-dev libarmadillo-dev libhdf5-dev
```

## Installation
Installation guide for Debian/Ubuntu Linux. It has not been tested on any other linux distros. 

(these insturctions include installation of mathplot in-tree)
```bash
#cloning my code
git clone https://github.com/MJ6Z/diffusion
#cloning mathplot in tree.
cd diffusion/
git clone --recurse-submodules git@github.com:sebsjames/mathplot
#building the makefile.
cmake -B build/
#compliling my code.
cd build/
make
```

## Running the program
From your install directory (not your build directory) run:
```bash
./build/model params.json
```

## changing the programs parameters
The parameters the model runs off are all based in [params.json](/params.json) they have brief descriptions.

## disclaimer
This is an student project, and is not intented for real world use or development. But you're still welcome to have a look at what I've been getting up to!

### README.md to be updated in the future.
