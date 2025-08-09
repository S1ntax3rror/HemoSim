# Steps to be taken:

## Initial setup --> Running simulation
- create a base setup
- allready specify the simulation time
- modify the run.sh file to copy all nescessary files (ensure to not copy dcd files)
- configure the Observer.py as needed
- use the build_setups_from_main_simulation.py script to generate multiple initial conditions
- use the HemoSetupGenerator.py script to generate all required permuatations
- use the run_multiple_observers.sh script to launch all Observers required
- monitor the Observers via the observer.log file

## Monitoring:
- use nvidia-smi to check GPU utilization
- if less than expected stop simulations and go back to building a solid initial setup
- use squeue -u [username] and make sure that all observers are running properly

## Patching after setup:
- The deploy_observer.py script allows to deploy scripty into all folders automatically.


## TODO

- update NPT inp file to also vary sim lengths !!
- check sim container / NPT / 10ns / 0000 -- (step 4 crashed for some reason)


- verify NVT runs
- setup NVT env
- start test run

- verify test run
- start main run