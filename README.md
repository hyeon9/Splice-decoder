# What is the Splice-decoder?
* Splice decoder provides functional annotation for your differential splicing events (DESs) or any splicing events of interest
* The functional annotation contains NMD probability, functional domain alteration (such as DNA binding, motif, regions, protein domain, and so on), CDS/UTR alterations, and Effect Score
* You can use the effect score to prioritize your DSEs
* Now we only support hg38 (hg19 will be supported soon)

# Workflow overview
<img width="1401" alt="image" src="https://github.com/user-attachments/assets/bf9fcb4d-d8ae-499c-8a1c-25eb232f520d" />

* It consist of five  main steps
  1. Processing input (Generate All Possible Splicing Cases): This step makes proper format of input data from the output of event-based splicing tools
  2. Mapping DSEs and ORFs (Map Splicing Cases): This step explores the given transcriptome (.GTF) to find Ref-TX (Reference transcript, it contains perfectly matched exon structure for the given DSE) and assign the best three open reading frames (ORFs)
  3. Simulation (Simulate Splicing Events): Based on the Ref-TXs and their ORFs, this step perform simulation of alternative splicing (e.g., if the Ref-TX has exon inclusion (EI) form, this step makes a simulated transcript (Sim-TX) with exon skipped (ES) form)
  4. Annotation (Functional Annotation): Based on the Uniprot DB, SpliceDecoder assigns known functional domains and estimates functional changes between Ref-TX and Sim-TX
  5. Scoring (DSEs with Effect Score): SpliceDecoder assigns an effect score to each DSE based on multiple biological factors, enabling prioritization of your DESs

<br>

# Install & Usage
## Quick start (`For HPC users`)
* You can run SplicDecoder without install

      wget https://github.com/hyeon9/Splice-decoder/archive/refs/heads/main.zip

      OR

      cp -r /projects/anczukow-lab/kangh/Tool/Splice-decoder_git ./

      vi paths.config  # Update several variables
  
      sbatch Main.sh Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all ${Your_config}  # You can use all or specifiy certain function e.g., Make_input
  
* If you need some information to update `paths.config` file, [you can find it](#build-your-configuration-and-input-file):
  
<br>

## Quick start (`For Non-HPC users`)
* Splice-decoder can be downloaded from https://github.com/hyeon9/Splice-decoder/
* Before run the install script, user should install mamba or conda (we strongly recommend using mamba)
* If you are using a mamba, you can run this commend
  
      bash install.sh
  
* If you are using a conda, you can run this commend
  
      bash install_conda.sh

* Verify if splice-decoder is running properly

      mamba(or conda) activate splice-decoder
      bash Main.sh paths.config all

* If you want to run certain step, specify certain step

      mamba(or conda) activate splice-decoder
      bash Main.sh paths.config {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring}

* If you use SLURM, modifying configure file and using this command

      sbatch Main.sh paths.config {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all}

* All your output will be saved to `${input}/result`

<br>

## Build your configuration and input file
- You should specify your paths for variables `Main` (Splice-decoder install directory), `conda` (Conda path), `input` (It points your working directory: all processed files and outputs will be created here) in the `paths.config` file
- If you are using your own conda, you can find it using the command below and update `conda` variable in the `config file`

      conda activate splice-decoder
  
      conda info | grep "active env location"
  
- You should check whether your gtf file has `geneID` or `geneSymbol`. Then set the `geneID_type` in the paths.config file (doesn't need it anymore)
- You should specifiy the `species` as either `human` or `mouse` and specify the `seq_type` as either `SR` or `LR`
- If you have your own TPM matrix (sample * transcript) you should specify the `tpm` as `Y_own`. Otherwise, set the `tpm` as `Y` in the `paths.config` file
      
<br>

## Overview for output files

* You can ckech Summary HTML in `${SpliceDecoder_folder}/DS_input/figure`
  
![image](https://github.com/user-attachments/assets/fdfff5a3-b923-45d7-b5aa-4dfcd932c767)

* After completing the entire process successfully, you can find several output files in `${SpliceDecoder_folder}/DS_input/result`
  * `*Main_output.txt`
  * `*NMD_check.txt`
  * `*Domain_integrity_indi.txt`
  * `Whole_DS_score_Whole.txt`

* First of all, check `Whole_DS_score_Whole.txt`, which contains [You can find more details here](#description-for-output):
  * `Effect_Score`
  * `Domain_change_rate`
  * `Probability_of_NMD`
  * `DOA_direction`
  

* The `Domain_integrity_indi.txt` includes:
  * `DS-TX pair ID`
  * `ORF priority`
  * `Domain information`
  * `Domain change ratio` (|Sim domain - Ref domain| / Ref domain) per domain block
  * `Change direction` (1: Gain of Function, -1: Loss of Function)
  
  The individual domain block change ratio measures the relative change size normalized by the reference domain block size.


## Description for Output
*Example of important values of `Whole_DS_score_Whole.txt`*
![image](https://github.com/user-attachments/assets/957b665a-829e-4885-afed-ea02a7a9cf8e)

### Key Metrics

- **`LongID`**: DS event ID
- **`gene`**: Gene Symbol
- **`Reference_transcript`**: Matched Transcript (==Ref_TX)
- **`Simulated_event`**: Simulated event
- **`Effect_Score`**: **A score to prioritize your DS events**
- **`Domain_change_rate`**: Average rate of domain changes in Sim-TX compared to Ref-TX
- **`Probability_of_NMD`**: NMD **(-1)**, PTC removal **(1)**, No NMD related event **(0)**
- **`DOA_direction`**: GoD (Gain of Domain), LoD (Loss of Domain), NMD, no_change, CDS_alts, and UTR_alts)
- **`Delta_PSI`**: PSI difference (group2 - group1), it came from rMATS
- **`Transcript_usage`**: Proportion of expression of reference transcript for each gene

---
### Supplementary Metrics
- `ORF`: Used ORF (This file only contains ORF with the highest potential)
- `AUG (Ref-Sim)`: Start codon position on the Ref TX and Sim TX (Ref-Sim)
- `Stop`: Stop codon position on the Ref TX and Sim TX (Ref-Sim)
- `Delta_Amino_acid`: Amino acid length difference (Ref TX - Sim TX)
- `5'UTR_difference`: 5' UTR length difference (Ref TX - Sim TX)
- `3'UTR_difference`: 3' UTR length difference (Ref TX - Sim TX)
- `Domain_integrity`: (Sim_domain_length / Ref_domain_length) * 100
- `Length_of_simulated_tx_domain`: Total domain length of Sim TX
- `Length_of_referece_tx_domain`: Total domain length of Ref TX
- `rMATS_FDR(-log10)`: -Log10 scale FDR, it came from rMATS

## Visualize your DS simulation
* Based on your Main_output file, you can pcik ceratin DS event to visualize it using this code

      cd ${SpliceDecoder_folder}
  
      python code/02-3_v3_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5
![image](https://github.com/user-attachments/assets/d700e8c1-efb6-40a9-bb12-b6576b955ef1)


* If you want to `remove some information` in figure space, using `ri` option
  
      python code/02-3_v3_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5 -ri region
![image](https://github.com/user-attachments/assets/f5866dfd-2ca5-4bf5-a40d-54ad3bf6e827)



* All figures will be saved at `${input}/figure/consequence/`
