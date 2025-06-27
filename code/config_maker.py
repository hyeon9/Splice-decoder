#!/usr/bin/env python3
import questionary
import os
import sys

config = {}

output_prefix = questionary.text("Specify your config file name").ask()
config['Main'] = questionary.text("Enter the path of SpliceDecoder").ask()+"/"
config['code'] = os.path.join(config['Main'], "code/")
config['cpatdb'] = os.path.join(config['Main'], "dat/")
config['input'] = questionary.text("Enter your working directory").ask()+"/"
config['Your_rMATS'] = questionary.text("Enter your rMATS output path").ask()+"/"
config['Your_GTF'] = questionary.text("Enter your GTF file that you used in rMATS with its full path").ask()
config['Your_TPM'] = questionary.text("Enter your TPM matrix with full path").ask()
config['bam_list'] = questionary.text("Enter your bamlist which should contains bamfile with their full path in each line").ask()
config['species'] = questionary.text("Enter a species of your data").ask()
config['seq_type'] = questionary.text("Enter a sequencing type of your data e.g., SR (short-read) or LR (long-read)").ask()
config['FDR'] = questionary.text("Enter a FDR cut off for your rMATS").ask()
config['PSI'] = questionary.text("Enter a |dPSI| cut off for your rMATS").ask()
config['njobs'] = questionary.text("Enter a number of cpu in splice-decoder job").ask()
config['tpm'] = questionary.text("Specify your TPM option e.g., Y (Calculate TPM by using SpliceDecoder), Y_own (Using your own TPM matrix), and N (Not using TPM)").ask()

config['conda'] = sys.argv[1]
config['cpat'] = os.path.join(config['conda'], "bin/cpat")
config['bedtools'] = os.path.join(config['conda'], "bin/bedtools")
with open(f"{config['Main']}{output_prefix}.config", "w") as f:
    for key, value in config.items():
        f.write(f'{key}="{value}"\n')

print(f"{config['Main']}{output_prefix}.config has been created! \nYou can use `bash Main.sh {output_prefix}.config $Function` or `sbatch Main.sh {output_prefix}.config $Function` in your terminal.")
print(f"You can get more information for the Function variable of Main.sh by using bash `Main.sh -h`")
