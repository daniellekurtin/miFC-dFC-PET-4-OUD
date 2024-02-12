## miFC-dFC-PET-4-OUD ##
# Multimodal differences in brain function between controls and people with Opiate Use Disorder. #
This study investigates whether there are generalisable differences in connectivity among reward, attention, and cognitive networks in methadone-dependent (MD) people with Opioid Use Disorder (OUD) and healthy controls. We assessed the generalisability of between-group differences through analysis of two tasks - a Cue Reactivity task, with drug-related stimuli, and the Monetary Incentive Delay (MID) task, which has non-drug rewards. This work is a part of the Neural Correlates of Reward and Emotion (NCORE) study: https://www.imperial.ac.uk/brain-sciences/research/psychiatry/ncore/

Please reach out to danielle.kurtin18@imperial.ac.uk for any questions or feedback about this project.

## Repository contents ##
The code in this repository is to make the figures in the paper and run the following analysis: 
<ol>
<li> Compute pairwise mutual information functional connectivity (miFC) </li> 
<li> Assess group differences in miFC and their relationship to functional networks </li>   
<li> Characterise whether differences in miFC relate to μ-Opiod Receptor (MOR) or Dopamine D2 Receptor (DRD2) availability </li> 
<li> Compute a timeseries of brain states via: </li>
  <ol>
  <li> A novel dynamical functional connectivity analysis, HomeBrew State Dynamics </li> 
  <li> Leading Eigenvector Dynamic Analysis </li> 
  </ol>
<li> Compute the following metrics of brain state dynamics: </li> 
  <ol>
  <li> State lifetime </li>
  <li> State probability </li>
  <li> Lempel Ziv Complexity of state timeseries </li>
  <li> Block Decomposition Methods of Complexity of state timeseries </li>
  <li> 0-4th order transition entropy </li>
  </ol>
<li> Assess the effect of group on brain state dynamics </li> 
</ol>

## Methods figure ##

![Methods](https://github.com/daniellekurtin/miFC-dFC-PET-4-OUD/assets/45391054/03f1df69-0a92-4259-a134-0905f384b794)

## Folder description ##
### AtlasAndOtherInputs: ### 
PET images and text files with names and coordinates for ROIs. Some of the PET images are too large to upload. However, all are available at Justine Hansen's repository: https://github.com/netneurolab/hansen_receptors/tree/main. All cortical regions are defined using the Yeo/Schaefer 200-region parcellation, with all volumetric atlas, region names, and coordinates found in this repo: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal. Finally, the dataframe used as an input for the miFC and dFC analysis (df.mat) is too big to upload. Should you want this data, please reach out to danielle.kurtin18@imperial.ac.uk.

### Functions: ### 
All functions necessary to run the analysis. The subfolder Functions/complexity has all complexity metrics, which are used to compute information theoretic metrics of brain state dynamics. The Functions/mi folder contains scripts necessary to compute mutual information functional connectivity based off Peng et al 2015, which can be found here: https://uk.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation 

### LEiDA: ### 
A misleading folder name, as this contains code for classic LEiDA (created by Joana Cabral, sourced here: https://github.com/juanitacabral/LEiDA) and my HomeBrewStateDynamics. This contains code to run the analysis and make the figs. Some scripts are task-specific; however, they are all relatively interchangeable with some small tweaks in naming. 

### MI: ### 
All code to compute miFC, relate to PET data, and make figures.





