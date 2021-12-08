# dynCogPD

## Description
The main objective of this project is to investigate the time-varying modulations in dynamic brain networks states through 'dynamic-source-connectivity' approach combined with a 'dimensionality reduction method' (Independent Component Analysis: ICA) using HD-EEG signals.
In addition, an approach based on micro-state metrics is implemented in order to quantify differences between two groups of subjects (i.e., controls vs patients).

Using this code, you will be able to:

- Compute Sources (weighted minimum norm estimate: wMNE)\
- Compute dynamic Functional Connectivity dFC (Phase-Locking Value: PLV with sliding window)\
- Compute dynamic states (Independent Component Analysis: ICA-JADE)\
- Compute statistics to extract significant 'group-level' states\
- Compute microstates parameters (backfitting code) to extract significant differences between two groups (i.e., Controls vs Patients)\


## Structure

This project is structured into the following folders:

- Inputs  : includes necessary input variables (to be used in the code, automatically loaded in the code)
- Results : includes two subfolders : conn-PLV (to save connectivity results) and state-ICA (to save ICA results), each divided into HC (control group) and PD (patient group)
- Code    : includes two subfolders : Function folder (that contains necessary functions used at each step) and Main Scripts folder (that contains the two main scripts to run called : 'run_pipeline.m' and 'run_microstats.m')


## TODO

1- Download and unzip the project folder 'DynCogPD' (preserve folders and subfolders structure)

2- Do NOT manually add folders to matlab path (to prevent complications between some matlab and toolbox functions). Necessary inputs and codes to be used are automatically added to Matlab path through code.

3- Download Toolbox folder from [here](https://github.com/judytabbal/dynCogPD/releases/tag/v1), unzip it, and add it to the path of the project: dynCogPD\Toolbox

4- Download the realigned MRI structure and the original NII file of icbm MRI from [here](https://github.com/judytabbal/dynCogPD/releases/tag/v1) and add it to the path of the project inputs: dynCogPD\Inputs\icbm

5- Open 'run_pipeline.m' script in DynCogPD\Code\Main Scripts folder.

6- Run the previous script (run_pipeline.m) for each group seperately and progressively (controls then patients for example).
	It includes 6 sections, to run Section by Section:

	Section 1: (MANIPULATE appropriately Section 1.1. for paths definition)
                   Define General Variables to be used and automatically add necessary paths and variables.

	Section 2: (MANIPULATE appropriately Section 2.1. for database configuration)
                   Create data structure (fieldTrip format) to be used.

	Section 3: (RUN AS IT IS)
                   Compute HeadModel using OpenMEEG (common for all subjects, using template MRI)

	Section 4: (MANIPULATE appropriately Section 4.1. for source and dFC configuration)
                   Compute sources wMNE (filters output) and connectivity PLV with sliding window (cmat output) for the group of subjects.

	Section 5: (MANIPULATE appropriately Section 5.1. for ICA configuration)
                   Compute ICA results using JADE method. + Plot option

	Section 6: (MANIPULATE appropriately Section 6.1. for null distribution configuration)
                   Build null distribution based on sign flip. + Plot option

7- Open 'run_microstats.m' script in DynCogPD\Code\Main Scripts folder.

8- Run the previous script (run_microstats.m) once to extract microstates parameters relative to both groups.
	It inludes 5 sections, to run Section by Section:

	Section 1: (MANIPULATE appropriately Section 1 variables)
                   Define Global Variables to be used: including band of interest, paths, number of subjects per group, number of states per group.

	Section 2: (RUN AS IT IS)
                   Load already saved results from the first script (PLVs, ICAs, Perms) for both groups

	Section 3: (RUN AS IT IS)
                   Extract automatically significant states for both groups

	Section 4: (RUN AS IT IS)
                   Apply Backfitting approach for both groups to assign for each subject and at each temporal window the 'most similar' group state

	Section 5: (RUN AS IT IS)
                   Extract Microstates parameters for both groups

	Final Output: variables 'microparams_HC' and 'microparams_PD' structures containing 5 main parameters:
		(1) fraction coverage time
		(2) frequency of occurrence
		(3) average lifespan or duration
		(4) Global Explained Variance (GEV)
		(5) Transition between states (symmetric 'TRsym0' and non-symmetric 'TR' versions)


## Notes

Note A. The user can refer to the script called code_for_inputs.m to see the computation of some saved variables in Inputs Folder:
- mri_realign (realigned template MRI ICBM)
- scout_mni, scout_scs (position/orientation of destrieux centroids in mni/scs coordinates)
- subgrid (source grid format used)
- elec_BS_mm (EEG electrodes fieldTrip format)
- indata.mat and nb_trials_persub.mat (input preprocessed EEG trials to be segmented and converted into fieldTrip data structure in 'run_pipeline.m' script.\

Note B. indata.mat and nb_trials_persub.mat saved variables serve only as examples to understand correct structure fields and dimensions. However, the values are random and do not refer to the real EEG data values treated. The user should create his own.\

Note C. The user can visualize results of brain networks states with temporal evolution with/without null distribution: uncomment the code lines at Sections 6.4 / 5.4 of run_pipeline.m.
