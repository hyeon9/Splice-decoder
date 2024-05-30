# What is the main purpose of Splice-decoder?
* Splice decoder can add biological information to your differential splicing or your interesting splicing target list.
* This information contains NMD probability, functional domain alteration (such as DNA binding, motif, regions, protein domain, and so on) and UTR alterations as their target transcript.
* Finally you can use this additaional information to prioritize your differential splicing events
* Now we only support hg38 (hg19 will be supported soon)

## Workflow overview
![image](https://github.com/hyeon9/Splice-decoder/assets/51947181/37692184-60f3-48d8-8e91-bc03f5596d69)


## Install & Usage
**Quick start:**  
* Splice-decoder can be downloaded from https://github.com/hyeon9/Splice-decoder/
* Before run the install script, user should install mamba or conda (we strongly recommend using mamba)
* If you are using a mamba, run this commend
  
       bash install.sh
  
* If you are using a conda, run this commend
  
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
  
## Build your own configuration file to run Splice-decoder
- You should prepare rMATS output file, gtf file (which was used in rMATS), and RNAseq bam file
- If you want to use tpm normalized counts in the effect score calculation step, you should set the tpm as Y of paths.config file
- If you used long-read RNA seq, you should check whether your gtf file has geneID or geneSymbol. Then set the geneID_type and seq_type of paths.config file
- You should set your Main (where Splice-decoder main directory), conda (where conda directory), input, and rMATS_path of the paths.config file
- All results were saved in ${input}/result
- Your gtf file should be named as "main.gtf", using this

       ln -s ${Your_gtf} ${input}/main.gtf

## Input Format
* Splice-decoder use rMATS JECE outputs, usnig this commend splice-decoder make proper input format from your rMATS output path
  
       bash Main.sh Make_input

* If you want to modify significant level, you can change FDR and dPSI var
* Splice-decoder also offers gtf processing and TPM calculate fucntion, each script requires own config file

       bash gtf_proc.sh ${config}
       bash stringtie.sh ${config}

## Simulation Output files
* Simulation analysis makes three output files (Main_output.txt, NMD_check.txt, and Domain_integrity_indi.txt)
* `Main_output.txt` contains overall information (e.g., NMD probability, average domain change ratio, UTR/CDS alteration, start/stop codon positions, and total domain length of each differential splicing (DS) and transcript (TX) pair
* `NMD_check.txt` contains DS-TX pair ID, ORF, NMD score (if it is larger than 50, it is considered an NMD), total domain length, and TX type (Reference or Simulation)
* `Domain_integrity_indi.txt` contains DS-TX pair ID, ORF, domain information, and domain change ratio (|Sim domain - Ref domain| / Ref domain) per domain block - The individual domain block change ratio measures relative change size by normalized reference domain block size

## Description for Main_output
- LongID: DS event ID
- Target_TX: Matched Transcript (==Ref_TX)
- occurred_event: Simulated event
- ORF_priority: pORF1 has the highest coding potential
- Start: [Ref TX start codon-Sim TX start codon]
- Stop: [Ref TX stop codon-Sim TX stop codon]
- 5'UTR: 5' UTR length difference (Ref TX - Sim TX)
- dAA: Amino acid length difference (Ref TX - Sim TX)
- 3'UTR: 3' UTR length difference (Ref TX - Sim TX)
- Domain_integrity: (Sim_domain_length / Ref_domain_length) * 100
- Domain_change_ratio: Average domain change ratio (average of |Ref TX - Sim TX| / Ref TX for individual domain)
- Ref_domain_length: Total domain length of Ref TX
- Sim_domain_length: Total domain length of Sim TX
- pNMD: -1 = NMD, 1 = PTC remove, 0 = No NMD related event
