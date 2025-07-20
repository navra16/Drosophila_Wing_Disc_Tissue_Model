# DSP_Wing_Disc_Eversion
Tissue model. The initial condition is an open system

Case 1. 3rd Instar Larval Shape

      Starting from an initial flat cylindrical triangulated model representing a piece of tissue, a selected region can undergo changes in mechanical properties leading to deformation of the tissue into a dome shaped starting structure for the end of this stage which will form the initial conditions for Case 2.  

Case 2. Eversion

      Starting from an initial domed triangulated model representing a piece of tissue, selected regions can undergo changes in mechanical properties leading to deformation and formation of a bilayer then extension to form the fully formed fly wing. 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

The general flow of simulation steps can be found at "void System::solveSystem()" in System.cu.

For particular functions such as linear spring, please refer to LinearSpring.cu and LinearSpring.h files. The same applies to bending spring and area springs.

For edge-swap algorithm and general data structure manipulation functions, please refer to Edgeswap_test.cpp and Edgeswap_test.h.

To change the name of saved animation and data output, please refer to Storage.cpp and Storage.h.

To change the simulation job title and to some extend simulation time step size, please refer to SBATCH.sh.

Initial data structure (built via MATLAB functions) is located in Data_Structure.xml.

Overall flow of the simulation steps:

I. Initialization of global parameters and data structures.

II. Run a predetermined number of relaxation steps of the model system to attain quasi-steady state.

III. Start the actual simulation:

      1. Update parameters (if necessary).
      2. Run a predetermined number of relaxation steps or dynamical number of relaxation steps depending on simulation types (molecular dynamics vs energy minimization).
      3. Run edge-swap algorithm (if applicable).
      4. Repeat (1-3) for a number of times, and then test for growth (if applicable).
      5. Repeat (1-4) until simulation terminates.

To run the simulation on UCR HPCC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load extra; module load GCCcore/4.9.3; module load cuda/9.1
4. type: make (Or make -j N, N can be 2,3,4,....,or 12. But this should only be done if you are using an interactive gpu session. See UCR HPCC website for detail)
5. type: sbatch -p gpu --gres=gpu:1 --time=144:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh 

###################################################



Current working simulation steps (modified from above):
To run simulation on UCR HPCC cluster:
1. *either clone or upload file onto HPCC user storage.*
2. *change directory to where the "SBATCH_try_this_one_if_the_original_does_not_work.sh" is located using*:
      a. cd [*folder name with "SBATCH_try_this_one_if_the_original_does_not_work.sh" file*]
3. module load singularity
4. module load centos
5. centos.sh
6. module load extra; module load GCCcore/4.9.3; module load cuda/9.1;
7. module load cmake
8. make
9. *After compilation completes, enter*: exit;
10. sbatch -p gpu --gres=gpu:1 --time=X:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh (*X here is the number of hours you want to keep the simulation running*)

