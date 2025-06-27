# What is the Splice-decoder?
![image](https://github.com/user-attachments/assets/fe9dde61-83e3-4208-9fb2-d3887266b004)


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

# Install & Usage
## Quick start (Before install SpliceDecoder, user should install conda first)
* Splice-decoder can be downloaded from https://github.com/hyeon9/Splice-decoder/
* To install SpliceDecoder, you should specify propoer yml file for your OS (`For_Mac_user.yml` or `For_Linux_user.yml`)
  
      bash install_conda.sh ${yml}

* To verify if splice-decoder is installed properly, you can do a test run with the toy_data (subset of rMATS and GTF file)

      conda activate splice-decoder

* To run SpliceDecoder with the toy_data build your configuration file through an interactive way

      cd code/
      bash Make_config.sh

* If you successfully created `Your.config`, now you can run SpliceDecoder
* The steps are intended to be executed in order, so it is recommended to use 'all')

      bash Main.sh ${Your.config} all

* If needed, you can run a specific step by selecting one of the following: `Make_input`, `DS_mapping`, `ORF_mapping`, `Simulation` and `Scoring`

      bash Main.sh ${Your.config} {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring}

* If you use SLURM, modifying configure file and using this command

      sbatch Main.sh ${Your.config} {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all}

* All your output will be saved to `${input}/result`

<br>

## For JAX users
* You can run SplicDecoder without install

      cp -r /projects/anczukow-lab/kangh/Tool/Splice-decoder_git ./

      vi paths.config  # Update several variables
  
      sbatch Main.sh Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all ${Your_config}
  
<br>

# Overview for output files

## Summary stats for DES mapping and annotation
* You can ckech Summary HTML in `${SpliceDecoder_folder}/DS_input/figure`
  
![image](https://github.com/user-attachments/assets/fdfff5a3-b923-45d7-b5aa-4dfcd932c767)

* After running the entire process successfully, you can find several output files in your `${SpliceDecoder_folder}/DS_input/result`
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
  * `Domain change ratio`
  * `Change direction`

<br>

## Descriptions for Outputs
*Example of important values of `Whole_DS_score_Whole.txt`*
![image](https://github.com/user-attachments/assets/da4dca2d-7c94-4277-b155-8261bd57a8ad)

### Key Metrics

- **`LongID`**: DS event ID
- **`gene`**: Gene Symbol
- **`Reference_transcript`**: Matched Transcript (==Ref_TX)
- **`Simulated_event`**: Simulated event
- **`Effect_Score`**: **A score to prioritize your DS events**
- **`∆L (Functional_change_ratio)`**: Average rate of domain changes in Sim-TX compared to Ref-TX
- **`Probability_of_NMD`**: NMD **(-1)**, PTC removal **(1)**, No NMD related event **(0)**
- **`DOA_direction`**: It contains the following functional classes: GoD (Gain of Domain), LoD (Loss of Domain), NMD, no_change, CDS_alts, and UTR_alts
  
### Supplementary Metrics
- `Delta_PSI`: PSI difference (group2 - group1), it came from rMATS
- `Transcript_usage`: Proportion of expression of reference transcript for each gene
- `ORF`: Used ORF (This file only contains ORF with the highest potential)
- `AUG (Ref-Sim)`: Start codon position on the Ref TX and Sim TX (Ref-Sim)
- `Stop_codon (Ref-Sim)`: Stop codon position on the Ref TX and Sim TX (Ref-Sim)
- `Delta_Amino_acid`: Amino acid length difference (Ref TX - Sim TX)
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

      cd ${SpliceDecoder_folder}
  
      python code/02-3_v3_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5
![image](https://github.com/user-attachments/assets/d700e8c1-efb6-40a9-bb12-b6576b955ef1)


* If you want to `remove some information` in figure space, using `ri` option
  
      python code/02-3_v3_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5 -ri region
![image](https://github.com/user-attachments/assets/f5866dfd-2ca5-4bf5-a40d-54ad3bf6e827)



* All figures will be saved at `${input}/figure/consequence/`

<br>

## Create 3D Protein structure based on simulated 
* You can use `Make_aa_fa.py` to extract amino acid sequences from your interesting targets
* You can find the ${input} variable in your `.config` file

      python code/Make_aa_fa.py -i ${input} -r human -t ENST00000438015.6 -e ES

* You can copy and paste the amino acid sequences into the Alphafold server (https://alphafoldserver.com) as input
