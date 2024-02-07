# miFC-dFC-PET-4-OUD
## Multimodal differences in brain function between controls and people with Opiate Use Disorder. ##
This study investigates whether there are generalisable differences in connectivity among reward, attention, and cognitive networks in methadone-dependent (MD) people with Opioid Use Disorder (OUD) and healthy controls. We assess the generalisability of between-group differences through analysis of two tasks - a Cue Reactivity task, with drug-related stimuli, and the Monetary Incentive Delay (MID) task, which has non-drug rewards. This study is a part of the Neural Correlates of Reward and Emotion (NCORE) study: https://www.imperial.ac.uk/brain-sciences/research/psychiatry/ncore/

The code in this repository is to make the figures in the paper and run the following analysis: 
<ol>
<li> Compute pairwise mutual information functional connectivity (miFC) </li> 
<li> Assess group differences in miFC and their relationship to functional networks </li>   
<li> Characterise whether differences in miFC relate to μ-Opiod Receptor (MOR) or Dopamine D2 Receptor (DRD2) availability </li> 
<li> Compute a timeseries of brain states via: </li>
  <ol>
  <li> A novel dynamical functional connectivity analysis </li> 
  <li> Leading Eigenvector Dynamic Analysis </li> 
  </ol>ol>
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

A visual schematic of this pipeline is below:

![Methods](https://github.com/daniellekurtin/miFC-dFC-PET-4-OUD/assets/45391054/03f1df69-0a92-4259-a134-0905f384b794)




A note on dependencies and availability: 
Some of the PET images are too large to upload. However, all are available at Justine Hansen's repository: https://github.com/netneurolab/hansen_receptors/tree/main. 
Moreover, the dataframe used as an input for the miFC and dFC analysis is too big to upload. Should you want this data, please reach out to danielle.kurtin18@imperial.ac.uk. 
