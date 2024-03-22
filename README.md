# Nuclear Design Diffusion Model
A part of my A-level computer sciene NEA project. This code simulates the diffusion of neutrons over a hexagonal grid.

The model is heavily dependant on the [morphalogica](https://github.com/ABRG-Models/morphologica) library installed in-tree with the code.



## dependancies
Dependancies should be flagged up when required when building the program. Once development is finished I'll add a comprehensive list of dependancies.

## Installation
Installation guide for Debian/Ubuntu Linux. This project is not yet supported for other operating systems. It has not been tested on any other linux distros.
(this includes installing morphalogica in-tree)
```bash
#cloning my code
git clone https://github.com/MJ6Z/diffusion

#cloning and building morphalogica in tree.
cd diffusion/
git clone https://github.com/ABRG-Models/morphologica
cmake -B build/
cd build/
make
```

## Running the program
From your install directory (not your build directory) run:
```bash
./build/model params.json
```

## changing the programs parameters
The parameters the model runs off are all based in [params.json](/params.json)

## disclaimer
This is an A-level CS project, and is not intented for real world use or development, I am not profficent enough with c++ to guarantee the safety of the program. But you're still welcome to have a look at what I've been getting up to!

### README.md to be updated in the future
