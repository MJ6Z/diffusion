# Nuclear Design Diffusion Model
A part of my A-level computer sciene NEA project.

The idea of the program is to simulate the core of a nuclear reactor with a low level (but not no) realism using differential equations. This repo was orginally created as a source for code that simulates divergence-gradient diffision of one variable over a hexagonal grid.

The program allows the user to change a variety of factors in [params.json](/params.json) aswell assign certain hexagons special properties. These being nuclear fuel rods, control rods, or coolant channels as well as neutron source locations.

The model is heavily dependant on the [morphologica](https://github.com/ABRG-Models/morphologica) library installed in-tree with the code.


## dependancies
This project requires an installation of Morphologica's dependancies.

Extract from [morphologica/README.md]()
```bash
sudo apt install build-essential cmake git wget  \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev \
                 libglfw3-dev libfreetype-dev libarmadillo-dev libhdf5-dev
```

## Installation
Installation guide for Debian/Ubuntu Linux. It has not been tested on any other linux distros. Morphalogica is supported on Linux, Mac and Windows - but this project is not.

(these insturctions include installation morphalogica in-tree)
```bash
#cloning my code
git clone https://github.com/MJ6Z/diffusion
#cloning morphalogica in tree.
cd diffusion/
git clone https://github.com/ABRG-Models/morphologica
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
This is an A-level CS project, and is not intented for real world use or development, I am not profficent enough with c++ to guarantee the safety of the program. But you're still welcome to have a look at what I've been getting up to!

### README.md to be updated in the future.
### Code docs will be available after this project is finished.
