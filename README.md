# TO DO
- [] Make a detailed dscription for output files
- [X] Make HTML summary using sort of figures
- [] Make a description for vis function


# What is the main purpose of Splice-decoder?
* Splice decoder can add biological information to your differential splicing or your interesting splicing target list
* This information contains NMD probability, functional domain alteration (such as DNA binding, motif, regions, protein domain, and so on) and UTR alterations as their target transcript
* You can use this additaional information to prioritize your differential splicing events in your analysis
* Now we only support hg38 (hg19 will be supported soon)

# Workflow overview
![image](https://github.com/user-attachments/assets/57495564-6112-4127-a4fb-3a79192af674)

* It consist of four main parts
  1. Make_input: This part makes proper format of input data from the output of event-based splicing tools
  2. DS_mapping and ORF_mapping: This part explores the given transcriptome (GTF file) to find Reference TX (it contains perfectly matched exon structure for certain differential splicing event) and define Open-reading frame using CPAT2
  3. Simulation: Using discovered Reference TXs (Ref-TXs) and their ORF, this part simulates their countuer-part splicing event (e.g., If the Ref-TX has exon inclusion form, this part makes exon skipped form)
  4. Scoring: Using transcript usage, splicing likelihood, and ORF prioirty, this part calculate effect score to prioritize each differential splicing - transcript pair

<br>

# Install & Usage
## Quick start (`For HPC users`)
* You can run SplicDecoder without install

      wget https://github.com/hyeon9/Splice-decoder/archive/refs/heads/main.zip
      sbatch Main.sh paths.config {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all}
  
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
**Make your input directory**
- You should prepare rMATS output file, gtf file (which was used in rMATS)
- In this example, we used "DS_input" as a name of input directory
  
      cd ${SpliceDecoder_folder}
      mkdir DS_input/

- Then, make symbolic link of your rmat output and GTF file. They should have fixed name (e.g., rmat.csv and main.gtf).
- If you don't have rmat.csv file,
[you can make it](#post-processing-for-rmats-output-files)

      cd DS_input
      ln -s ${Your_rMATS} ./rmat.csv
      ln -s ${Your_GTF} ./main.gtf
  
- If you want to use tpm normalized counts in the effect score calculation step, you should set the tpm as Y of paths.config file (It IS HIGHLY RECOMMENDED)
- If you don't have TPM matrix,
[you can make it](#post-processing-for-rmats-output-files)

      cd ${SpliceDecoder_folder}
      cd DS_input
      mkdir tpm
      cd tpm
      ln -s ${Your_TPM} ./matrix.tpm

- If you used long-read RNA seq, you should check whether your gtf file has `geneID` or `geneSymbol`. Then set the geneID_type and seq_type of paths.config file
- You should put your path for variables `Main` (Splice-decoder install directory), `conda` (conda path), `input` (${SpliceDecoder_folder}/DS_input) of the paths.config file
- You can find your conda path

      conda activate splice-decoder
      conda info | grep "active env location"

- If you find your conda path, just modify your paths.config

      cd ${SpliceDecoder_folder}
      vi paths.config

<br>

## Post Processing for rMATS output files
* Splice-decoder use rMATS JECE outputs, usnig this commend splice-decoder make proper input format from your rMATS output path
  
      python ${SpliceDecoder_folder}/code/NEW_make_input_from_rmats.py ${Your_rMATS} ${SpliceDecoder_folder}/DS_input/

* You can make TPM matrix by using BAM files (

      bash ${SpliceDecoder_folder}/code/stringtie.sh ${SpliceDecoder_folder}/main.gtf ${bam_list} ${SpliceDecoder_folder}/DS_input/
      
* Your gtf file should be named as `main.gtf` using this

      ln -s ${Your_gtf} ${input}/main.gtf
  
* Splice-decoder also offers gtf processing and TPM calculate fucntion, each script requires own config file

      bash gtf_proc.sh ${config}
      bash stringtie.sh ${config}

<br>

## Description for output files

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
- **`DOA_direction`**: GoD (Gain of Domain), LoD (Loss of Domain), NMD, no_change, and other_retions_diff (e.g., UTRs and unannot CDS)
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
      python code/02-3_v2_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5
![image](https://github.com/user-attachments/assets/d700e8c1-efb6-40a9-bb12-b6576b955ef1)


* If you want to `remove some information` in figure space, using `ri` option
  
      python code/02-3_v2_Draw_consequence.py --input `pwd`/DS_input/ --splicing_event CA --gene SPG7 --sim_splicing_event EI --transcript ENST00000561911.5 -ri region
![image](https://github.com/user-attachments/assets/f5866dfd-2ca5-4bf5-a40d-54ad3bf6e827)



* All figures will be saved at ${input}/figure/consequence/
