# Kelvin-Helmholtz Simulations

This project is a collection of scripts used to run numerical simulations of Kelvin-Helmholtz instability. It uses [Dedalus](https://dedalus-project.org), an open source tool for spectrally solving differential equations.

To learn more about Kelvin-Helmholtz instability and how these scripts can be used to learn more about it, see the included paper about the topic [here](./paper/kelvin-helmholtz-instability.pdf).

## Getting Started
These simulations were run on Ubuntu 20.04 using Python 3.8, but newer versions of Ubuntu should work just as well. To set up an environment for running the simulations, do the following:
1. Make sure you have the Open MPI development package installed
    ```
    sudo apt install libopenmpi-dev
    ```
1. Create a virtual environment for running the simulations and install dependencies:
    ```
    python3 -m venv env
    source env/bin/activate
    pip install -r requirements.txt
    ```
    This should install everything you need to run the simulation and image processing scripts.
1. Before running any simulations, you'll need to activate the virtual environment:
    ```
    source env/bin/activate
    ```

## Scripts
All simulation and processing scripts can be found in the `scripts/` folder.
- **K-H-instability.py** numerically solves the equations of motion for Kelvin-Helmholtz instability with the given parameters, using free-slip boundaries on the top and bottom, and periodic boundaries on the sides.
- **interpolate_img_sets.py** interpolates the data produced in the simulation, producing high-resolution images.
- **make_movie.py** accepts a path to a list of frames and produces an animation of the fluid motion from the frames.
- **plot_slices.py** takes the data produced by `K-H-instability.py` and plots frames of the fluid motion without interpolation.

To see the available options for each script, run `python3 <script-name>.py --help`. These scripts are intended to be run with MPI, as shown in `run_sim.sh`.