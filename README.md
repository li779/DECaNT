DECaNT
=========================================
Project Contributers: Y. C. Li, A. H. Davoody, S.W. Belling, A. J. Gabourie, and I. Knezevic
Official implementation of [DECaNT: Simulation Tool for Diffusion of Excitons in Carbon Nanotube Films](https://arxiv.org/pdf/1703.05192.pdf). 

<img src="assets/discogan.png" width="600px">

Prerequisites
-------------
   - Python 3
   - Armadillo
   - BulletPhysics

Mesh Generation
----------------
### CelebA
Download CelebA dataset using

    $ python ./datasets/download.py celebA 

(Currently, the link for downloading [CelebA](http://mmlab.ie.cuhk.edu.hk/projects/CelebA.html) dataset is not available).

To train gender conversion,

    $ python ./discogan/image_translation.py --task_name='celebA' --style_A='Male'

To train hair color conversion 

    $ python ./discogan/image_translation.py --task_name='celebA' --style_A='Blond_Hair' --style_B='Black_Hair' --constraint='Male'

### Handbags / Shoes
Download Edges2Handbags dataset using 

    $ python ./datasets/download.py edges2handbags

Download Edges2Shoes dataset using 

    $ python ./datasets/download.py edges2shoes

To train Edges2Handbags,

    $ python ./discogan/image_translation.py --task_name='edges2handbags'

To train Edges2Shoes,

    $ python ./discogan/image_translation.py --task_name='edges2shoes' 

To train Handbags2Shoes,

    $ python ./discogan/image_translation.py --task_name='Handbags2Shoes' --starting_rate=0.5

### Facescrub
Download Facescrub dataset using 

    $ python ./datasets/download.py facescrub

To train gender conversion,

    $ python ./discogan/image_translation.py --task_name='facescrub'

### Car, Face
Download [3D car dataset](http://www.scottreed.info/files/nips2015-analogy-data.tar.gz) used in [Deep Visual Analogy-Making]( http://www-personal.umich.edu/~reedscot/nips2015.pdf), and [3D face dataset](http://faces.cs.unibas.ch/bfm/main.php?nav=1-2&id=downloads) into ./datasets folder and extract them.

To train Car2Car translation,

    $ python ./discogan/angle_pairing.py --task_name='car2car' 

To train Car2Face translation,

    $ python ./discogan/angle_pairing.py --task_name='car2face'

Run script.sh in order to train a model using other datasaet, after uncommenting corresponding line.

Monte Carlo Simulation
----------------
### Code Structure
Main Code is in monte_carlo folder. The whole simulation is divided into three major object: Simulation itself, exciton and scattering sites. Each object has its class with detailed attributions and methods. Simulation object is the top-layer client that will create exciton and scattering sites based on simulation metrics. The output file contains displacement and position of excitons at each time step.
There are complementary code in exciton transfer and helper folder. In exciton transfer folder, codes are used in calculating Carbon Nanotube bandstructure and resonant exciton transfering. Functions will be called by Monte Carlo Simulation code to generate scattering table.
