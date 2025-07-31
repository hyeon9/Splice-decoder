## Quick start (conda is required)
* transcript-toolkit automatically be installed with Splice-decoder
* You can create a configuration file using your own data to analyze your own data

      cd code/transcript-toolkit/
      bash Make_config.sh

* If you successfully created `Your.config`, now you can run transcript-toolkit
* The steps are intended to be executed in order, so it is recommended to use `all`

      bash transcript-toolkit.sh ${Your.config}

## Guide for making config file
* Make_config.sh will ask..

      ? Specify your config file name (e.g. HGjob)
      > You just need to specify your config file

      ? Enter the path of transcript-toolkit (e.g. /User/usr/Tool/transcript-toolkit-main/)
      > You just need to specify the install path of transcript-toolkit

      ? Enter your working directory (e.g. /User/usr/Tool/transcript-toolkit-main/project1)
      > You just need to specify your new working directory

      ? Enter your query list with its path (e.g. /User/usr/Tool/transcript-toolkit-main/query.tsv)
      > You just need to specify the query list path

      ? Enter your GTF file that you used in transcriptome assembly with its full path (e.g. /User/usr/Tool/transcript-toolkit-main/assem.gtf)
      > You just need to specify the full path + GTFfile

      ? Enter your data species (e.g. human or mouse)
      > You just need to specify the species of your data

      ? Enter your data sequencing type? [LR (Long-read) or SR (Short-read)]
      > You just need to specify the sequencing method of your data

      ? Enter a number of cpu in transcript-toolkit job (int [0-?])
      > Specify a number of cpu will be used in your job

 
## Descriptions for Outputs
*Example of important values of `Main_output.txt`*

### Key Metrics

- **`Comparison`**: Canonical vs Query transcript ID
- **`Reference_transcript`**: Canonical transcript ID
- **`Query_transcript`**: Query transcript ID
- **`Domain_change_rate`**: Average rate of domain changes in Sim-TX compared to Ref-TX
- **`Probability_of_NMD`**: NMD **(-1)**, PTC removal **(1)**, No NMD related event **(0)**
- **`Functional_class`**: It contains the following functional classes: GoD (Gain of Domain), LoD (Loss of Domain), NMD, CDS_alts, and UTR_alts
  
### Supplementary Metrics
- `Transcript_usage`: Proportion of expression of reference transcript for each gene
- `ORF`: Used ORF (This file only contains pORF1 which has the highest coding potential)
- `AUG (Ref-Sim)`: Start codon position on the Ref TX and Sim TX (Ref-Sim)
- `Stop (Ref-Sim)`: Stop codon position on the Ref TX and Sim TX (Ref-Sim)
- `5'UTR_difference (nt)`: 5' UTR length difference (Ref TX - Sim TX)
- `CDS_difference (nt)`: Amino acid length difference (Ref TX - Sim TX)
- `3'UTR_difference (nt)`: 3' UTR length difference (Ref TX - Sim TX)
- `Domain_integrity`: (Sim_domain_length / Ref_domain_length) * 100
- `Length_of_simulated_tx_domain`: Total domain length of Sim TX
- `Length_of_referece_tx_domain`: Total domain length of Ref TX


## Visualize your DS simulation
* Based on your Main_output file, you can pick ceratin comparison to visualize it using this code

      conda activate splice-decoder
      python code/transcript-toolkit/02_Vis.py --input ${working directory} --cano_tx ${Canonical_transcript} --query_list ${Your_query_list}

## Create a 3D Protein structure based on simulated 
* You can use `Make_aa_fa.py` to extract amino acid sequences from your interesting targets
* !! This function requires the `Main_output.txt`.
* You can find the `${input}` and `${Main}` in your `.config` file

      conda activate splice-decoder
      python code/transcript-toolkit/Make_aa.py -i ${input} -r human -rt ${ref_ENST}.6 -qt ${query_ENST} -d ${Main}

* You can copy and paste the amino acid sequences into the Alphafold server (https://alphafoldserver.com) as input
