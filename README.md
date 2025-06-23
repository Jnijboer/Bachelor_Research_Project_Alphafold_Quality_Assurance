
##### Quality assessment of AlphaFold structures through molecular dynamics

This project is about using Molecular Dynamics (MD) simulations to determine the accuracy of AlphaFold predictions and this repository contains all programs and files used for the analysis. 

*Group: Molecular Dynamics, University of Groningen.*
*Supervisor: Dr. Tsjerk A. Wassenaar*

Since the release of AlphaFold, it has become much easier to obtain protein structures. Instead of having to perform time-consuming experimental methods such as X-ray crystallography, cryo-EM or NMR, AlphaFold can now use Artificial Intelligence (AI) to make predictions about the structure. It does this with the help of the Protein Data Bank (PDB), which already contains many experimentally determined structures.Alphafold2 returns 5 predicted structures from one input sequence, but it does not indicate which of the 5 is the most accurate or most relevant in biological context. And when a protein can have multiple conformations, it makes the choice even more difficult. What has been done in this project is an attempt at using MD simulations to find out what features of the protein are a sign of an accurate prediction. This is based on the assumption that inaccurate proteins will show more 'stress' in the form of movements and internal force profiles. Therefore, this analysis focuses on visualising these patterns. The main test subjects that have been used in this analysis are the ALK kinase domain and the closely related ROS1 kinase domain.  

# Aims of this project
- processing binary trr files to read out the data of interest
- performing RMSD clustering analysis
- performing movement and force analysis 
- visualising the movements and changes in force on the kinase structures

# Contents
- trr processing program
- timer code
- automation code for the simulation 
- clustering analysis program
- movement/force analysis program
- vector placement program



