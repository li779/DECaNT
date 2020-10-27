DECaNT
=========================================
### Project Contributors
Y. C. Li, A. H. Davoody, A. J. Gabourie, S.W. Belling, and I. Knezevic

Official implementation of [DECaNT: Simulation Tool for Diffusion of Excitons in Carbon Nanotube Films](https://arxiv.org/abs/2010.11992). This project grew out of the work by former group members Amirhossein Davoody and Alexander J Gabourie. The initial project can be found here: [Mesh](https://github.com/amirhosseindavoody/carbon_nanotube_mesh), [Monte Carlo](https://github.com/amirhosseindavoody/cnt_film_monte_carlo).

<p align="center"><img src="graphs/Figure6_simulation_schematic.png" width="400px"></p>

Dependency
-------------
This project is intented to run on linux or linux-like operating systems. The following external libraries are required:
   - Python 3
   - Armadillo
   - BulletPhysics
### Armadillo
Armadillo is a linear algebra library. Before installing armadillo, we need to make sure BLAS and LAPACK are installed. Installation:

    $ sudo apt install libopenblas-dev liblapack-dev

Then we are ready to install Armadillo. There are two approaches: either go to the [website](http://arma.sourceforge.net/download.html) to obtain the tar ball or, directly using package management:

    $ sudo apt-get install libarmadillo-dev
    
### BulletPhysics
BulletPhysics is an open-source real-time physics simulation library. We use this library to generate carbon nanotube meshes with different morphologies. Installation requires the OpenGL library to render the scene. Detailed instructions can be found at the [old repo](https://github.com/amirhosseindavoody/carbon_nanotube_mesh/wiki)
   
Mesh Generation
----------------
### Code Structure
Properties of carbon nanotubes, such as length, segment length, tube chirality, etc., are stored in arrays that are properties of the class "cnt_mesh". The entire mesh of carbon nanotubes is an instance of the class "cnt_mesh". The methods defined in this class construct carbon nanotubes out of segments, apply constraints, and drop the nanotube from a specific height. Main.cpp will act as an upper-level structure to set up and run the BulletPhysics scene and also call functions in cnt_mesh to generate new tubes. The output file will record the position, orientation and chirality of each carbon nanotube segments. More detailed explanation can be found [Here](https://github.com/amirhosseindavoody/carbon_nanotube_mesh).

After generating the film, we use a python script to increase the density of gridpoints along the length of each tube using interpolation. Once the interpolation is finished, the script sets the gridpoints along each tube to act as scattering sites in the Monte Carlo simulation.
### Run example
First, navigate to the folder containing the project and use the makefile to generate the executable program.

    $ cd mesh
    $ make main
    
Then, modify input.json as needed, and run main.exe

    $ main.exe input.json

Run the python script.

    $ cd python_script
    $ python3 create_fine_mesh.py
    

Monte Carlo Simulation
----------------
### Code Structure
The most important code is contained in the monte_carlo folder. The code is divided into three major class structures: "simulation", "exciton", and "scattering site". Each class has its own properties and methods. The simulation object is used to control things like boundary conditions, simulation domain trimming, and the construction of exciton and scattering site objects. The output of the simulation is a file containing exciton displacement and position at each time step.

Other folders contain complementary code that is used to calculate exciton bandstructure and transfer rates. These will be called by the Monte Carlo simulation code as needed to generate the appropriate scattering tables.

Controlling the Monte Carlo simulation is best done using the input file.

### Run example
First, navigate to the folder containing the project and use the makefile to generate the executable program.

    $ cd montecarlo
    $ make main
    
Then, modify input.json as needed, and run main.exe

    $ main.exe input.json

