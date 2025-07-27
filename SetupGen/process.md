# Steps to be taken:

- create a base setup
- allready specify the simulation time
- modify the run.sh file to copy all nescessary files (ensure to not copy dcd files)
- configure the Observer.py as needed
- use the generate_branches script to generate multiple initial conditions
- use the HemoSetupGenerator.py script to generate all required permuatations
- use the run_multiple_observers.sh script to launch all Observers required
- monitor the Observers via the observer.log file
