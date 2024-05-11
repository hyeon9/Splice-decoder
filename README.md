# What is main purpose of Splice-decoder?
* Splice decoder can add biological information to your differential splicing or your interesting splicing target list.
* This information contains NMD probability, functional domain alteration (such as DNA binding, motif, regions, protein domain, and so on) and UTR alterations as their target transcript.
* Finally you can use this additaional information to prioritize your differential splicing events 

## Install & Usage
**Quick start:**  
* Splice-decoder can be downloaded from https://github.com/hyeon9/Splice-decoder/

       bash install.sh

* Verify if splice-decoder is running:

       python ${splice-decoder}/code/00-1_add_exon_n.py --h

* If you use SLURM, modifying configure file and using this command

       bash ${splice-decoder} 
  
## Make your own configuration file to run Splice-decoder
- Using pre-exist paths.config file
- You only modify your Main path (where Splice-decoder main folder) and conda path (where conda folder, you can find "which conda")
  
- If you finish the configuration, run bash script
- Proc.sh > Mapp.sh > Sim_splicing.sh
- All results were saved in ${input}/result

## INPUTS Format
* Splice-decoder requires your 

* Splice-decoder use rMATS JECE outputs
  
       python make_input_from_rmats.py ${rMATS_outdir}

* Splice-decoder also offers gtf processing and TPM calculate fucntion

       bash gtf_proc.sh ${config}
  


## How to USE


## Description for Main_output
LongID: DS event ID
Target_TX: Matched Transcript (==Ref_TX)
occurred_event: Simulated event
ORF_priority: pORF1 has the highest coding potential
Start: [Ref TX start codon-Sim TX start codon]
Stop: [Ref TX stop codon-Sim TX stop codon]
5'UTR: 5' UTR length difference (Ref TX - Sim TX)
dAA: Amino acid length difference (Ref TX - Sim TX)
3'UTR: 3' UTR length difference (Ref TX - Sim TX)
Domain_integrity: (Sim_domain_length / Ref_domain_length) * 100
Domain_change_ratio: Average domain change ratio (average of |Ref TX - Sim TX| / Ref TX for individual domain)
Sim_domain_length: Total domain length of Sim TX
Ref_domain_length: Total domain length of Ref TX
pNMD: -1 = NMD, 1 = PTC remove, 0 = No NMD related event
