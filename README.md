# Greenhill Observatory Scheduler

The UTAS Greenhill Observatory (UTGO) scheduler dashboard is a manual tool for planning a night of observing. It uses the Astroplan package to calculate the optimal observing times, given a list of targets and constraints. It will either read this from the 'UTGO_Targets' Google sheet, or a .CSV file placed in the folder. It then outputs a schedule and sky chart, hosted on a local server. Note, this tool is not for driving the telescope at the current time.

### Installation

It is recommended that the user downloads and installs Anaconda or Miniconda to manage their
python distributions. First, download/clone the Github repository to your local machine. You will need to change into the downloaded directory, using:

#### cd [insert path to folder here]

Then, you will need to install the necessary requirements to run the dashboard. To do this, go to
your terminal (or Anaconda Prompt on Windows) and create a new conda environment (for example
‘obs’):

#### conda create -n obs python=3.9.16

You will then need to activate the environment:

#### conda activate obs

Then, assuming that ‘pip’ installed into this environment during the process, type:

#### pip install -r requirements.txt

### How to Use

To run the dashboard, activate
your conda environment, change to the UTGO_Scheduler directory and type:

#### python scheduler.py

This will launch the dashboard on your local host. Go to a web browser and type the address given in
the command line.

For more detailed instructions and background, please see the 'UTGO_Scheduler.pdf' file in the 'Documentation folder. 
