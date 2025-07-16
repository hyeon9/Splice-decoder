## What is the Splice-decoder?
<img width="1290" alt="image" src="https://github.com/user-attachments/assets/7206811b-971a-45be-9221-3b4f5d7387f9" />

* Splice decoder provides functional annotation for your differential splicing events (DESs)
* The functional annotation contains NMD probability, functional domain alteration (such as DNA binding, motif, regions, protein domain, and so on), CDS/UTR alterations, and effect score
* You can use the effect score to prioritize your DSEs
* Currently, we only support hg38 (hg19 will be supported soon)

## Workflow overview
<img width="1401" alt="image" src="https://github.com/user-attachments/assets/bf9fcb4d-d8ae-499c-8a1c-25eb232f520d" />

1. Generate All Possible Splicing Cases (Processing input) : This step makes proper format of input data from the output of event-based splicing tools
2. Map Splicing Cases (Mapping DSEs and ORFs): This step explores the given transcriptome (.GTF) to find Ref-TX (Reference transcript, it contains perfectly matched exon structure for the given DSE) and assign the best three open reading frames (ORFs)
3. Simulate Splicing Events (Simulation): Based on the Ref-TXs and their ORFs, this step perform simulation of alternative splicing (e.g., if the Ref-TX has exon inclusion (EI) form, this step makes a simulated transcript (Sim-TX) with exon skipped (ES) form)
4. Functional Annotation (Annotation): Based on the Uniprot DB, SpliceDecoder assigns known functional domains and estimates functional changes between Ref-TX and Sim-TX
5. DSEs with Effect Score (Scoring): SpliceDecoder assigns an effect score to each DSE based on multiple biological factors, enabling prioritization of your DESs

<br>

## Quick start (conda is required)
* Splice-decoder can be downloaded from https://github.com/hyeon9/Splice-decoder/
* To install SpliceDecoder, you should specify proper yml file for your OS (`For_Mac_user.yml` or `For_Linux_user.yml`)
  
      cd ./Splice-decoder && bash install.sh ${yml}

* To perform a test run, you can use the provided toy_data
* To to this you should build a toy configuration file through an interactive way [You can find more details here](#guide-for-making-config-file):

      cd code/
      bash Make_config.sh toy

* Alternatively, you can create a configuration file using your own data to analyze your own data

      cd code/
      bash Make_config.sh

* If you successfully created `Your.config`, now you can run SpliceDecoder
* The steps are intended to be executed in order, so it is recommended to use `all`

      bash Main.sh all ${Your.config}

* If needed, you can run a specific step by selecting one of the following: `Make_input`, `DS_mapping`, `ORF_mapping`, `Simulation` and `Scoring`

      bash Main.sh {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring} ${Your.config}

* If you use SLURM, use this command to submit your job

      sbatch Main.sh {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all} ${Your.config}

* All your output will be saved to `${input}/result`

<br>

## Guide for making config file
* Make_config.sh will ask..

      ? Specify your config file name (e.g. HGjob)
      > You just need to specify your config file

      ? Enter the path of SpliceDecoder (e.g. /User/usr/Tool/Splice-decoder-main/)
      > You just need to specify the install path of SpliceDecoder
  
      ? Enter your working directory (e.g. /User/usr/Tool/Splice-decoder-main/project1)
      > You just need to specify your new working directory
  
      ? Enter your rMATS output path (e.g. /User/usr/Tool/Splice-decoder-main/toy_data)
      > You just need to specify the rMATS output path
  
      ? Enter your GTF file that you used in rMATS with its full path (e.g. /User/usr/Tool/Splice-decoder-main/toy_data/toy.gtf or /User/usr/Tool/Splice-decoder-main/toy_data/*.gtf)
      > You just need to specify the full path + GTFfile

      ? Do you want to calculate the effect score? [yes/no]
      > Simply type yes or no. If you type "yes", SpliceDecoder will ask TPM matrix or bamfile path to calculate the effect score

      ? Enter your TPM matrix with full path (e.g. /User/usr/Tool/Splice-decoder-main/toy_data/tpm.tsv or N)
      > Specify the full path to your TPM matrix, or enter 'N' if you don’t have one
  
      ? Enter your bamlist which should contains bamfile with their full path in each line (e.g. /User/usr/Tool/Splice-decoder-main/toy_data/bam_list.txt or N)
      > If you don’t have a TPM matrix, specify the full path to your BAM list file, or enter 'N'
  
      ? Enter a species of your data (e.g. human or mouse)
      > You just need to specify the species of your data
  
      ? Enter a sequencing type of your data (e.g. SR (short-read) or LR (long-read) )
      > You just need to specify the sequencing method of your data
  
      ? Enter a FDR cut off for your rMATS (float [0-1], default 0.05)
      > Specify rMATS FDR cut off
  
      ? Enter a |dPSI| cut off for your rMATS (float [0-1], default 0.1)
      > Specify rMATS FDR cut off
  
      ? Enter a number of cpu in splice-decoder job (int [0-?])
      > Specify a number of cpu will be used in your job

* You can reuse a pre-existing config file by copying it:

      cp ${existing_config} project2.config

* Then, update the following fields in the new config: `input`, `Your_GTF`, and `Your_rMATS`

      
  

<br>

## Summary stats for DES mapping and annotation
* You can ckech Summary HTML in `${SpliceDecoder_folder}/DS_input/figure`
  
![image](https://github.com/user-attachments/assets/fdfff5a3-b923-45d7-b5aa-4dfcd932c767)

* After running the entire process successfully, you can find several output files in your `${SpliceDecoder_folder}/DS_input/result`
  * `*Main_output.txt`
  * `*NMD_check.txt`
  * `*Domain_integrity_indi.txt`
  * `Whole_DS_score_Whole.txt`

* First of all, check `Whole_DS_score_Whole.txt`, which contains [You can find more details here](#descriptions-for-outputs):
  * `Effect_Score`
  * `Domain_change_rate`
  * `Probability_of_NMD`
  * `DOA_direction`
  

* The `Domain_integrity_indi.txt` includes:
  * `DS-TX pair ID`
  * `ORF priority`
  * `Domain information`
  * `Domain change ratio`
  * `Change direction`

<br>

## Descriptions for Outputs
*Example of important values of `Whole_DS_score_Whole.txt`*
<img width="1090" alt="image" src="https://github.com/user-attachments/assets/ef03e748-6fb3-4d07-9211-9a81b19e430b" />

### Key Metrics

- **`LongID`**: DS event ID
- **`gene`**: Gene Symbol
- **`Reference_transcript`**: Matched Transcript (==Ref_TX)
- **`Simulated_event`**: Simulated event (ES = Exon skipping, EI = Exon inclusion, SI = Skipped intron, RI = Retained intron, Can A3/5SS = canonical 3/5' splice site, Alt A3/5SS = alternative 3/5' splice site)
- **`Effect_Score`**: **A score to prioritize your DS events**
- **`Domain_change_rate`**: Average rate of domain changes in Sim-TX compared to Ref-TX
- **`Probability_of_NMD`**: NMD **(-1)**, PTC removal **(1)**, No NMD related event **(0)**
- **`Functional_class`**: It contains the following functional classes: GoD (Gain of Domain), LoD (Loss of Domain), NMD, CDS_alts, and UTR_alts
  
### Supplementary Metrics
- `Delta_PSI`: PSI difference (group2 - group1)
- `Transcript_usage`: Proportion of expression of reference transcript for each gene
- `ORF`: Used ORF (This file only contains pORF1 which has the highest coding potential)
- `AUG (Ref-Sim)`: Start codon position on the Ref TX and Sim TX (Ref-Sim)
- `Stop_codon (Ref-Sim)`: Stop codon position on the Ref TX and Sim TX (Ref-Sim)
- `CDS_difference`: Coding sequence length difference (Ref TX - Sim TX)
- `5'UTR_difference`: 5' UTR length difference (Ref TX - Sim TX)
- `3'UTR_difference`: 3' UTR length difference (Ref TX - Sim TX)
- `Domain_integrity`: (Sim_domain_length / Ref_domain_length) * 100
- `Length_of_simulated_tx_domain`: Total domain length of Sim TX
- `Length_of_referece_tx_domain`: Total domain length of Ref TX
- `rMATS_FDR(-log10)`: -Log10 scale FDR, it came from rMATS

---
*Example of important values of `Domain_integrity_indi.txt`*
![image](https://github.com/user-attachments/assets/266bcdec-d4aa-4e5b-a0bb-3bf0db247a14)

### Key Metrics
- **`DS-TX pair ID`**: It contains, in order Long_ID, Ref-TX ID, and simulated event type
- **`ORF priority`**: A priority of the used reading frame in simulation
- **`Domain information`**: A name of altered domain by the simulated alternative splicing event
- **`Functional_change_ratio (∆L)`**: A difference of functional change ratio for simulated alternative splicing
- **`Change direction`**: It indicates whether the altered domain is a gain (1) or a loss of domain (-1)

<br>

## Visualize your DS simulation
* Based on your Main_output file, you can pcik ceratin DS event to visualize it using this code

      conda activate splice-decoder
      python code/02-3_v3_Draw_consequence.py --input ${working directory} --splicing_event RI --gene ENSMUSG00000027470.9 --sim_splicing_event RI --transcript ENSMUST00000028970.7
      python code/02-3_v3_Draw_consequence.py -h  # You can get more details
![RI;ENSMUSG00000027470 9;chr2;+;152919325;152919453;152919454;152920285;152920286;152920438_ENSMUST00000028970 7_splicing_map](https://github.com/user-attachments/assets/9a9e69b1-5ae7-4efb-9228-47c829a0ff40)


* If you want to `remove some information` in figure space, using `ri` option (all categories should be separated by space)

      python code/02-3_v3_Draw_consequence.py --input ${working directory} --splicing_event A3SS --gene ENSMUSG00000028864.7 --sim_splicing_event Ori_A3SS --transcript ENSMUST00000195957.4 -ri proteome chain
![image](https://github.com/user-attachments/assets/3306e051-8fa7-47c3-ab25-13c6df061da1)



* All figures will be saved at `${input}/figure/consequence/`

<br>

## Create a 3D Protein structure based on simulated 
* You can use `Make_aa_fa.py` to extract amino acid sequences from your interesting targets
* !! This function requires the `Effect_score.tsv`, Toy data is not eligible for this function
* You can find the ${input} variable in your `.config` file

      conda activate splice-decoder
      python code/Make_aa_fa.py -i ${input} -r human -t ENST00000438015.6 -e ES

* You can copy and paste the amino acid sequences into the Alphafold server (https://alphafoldserver.com) as input
