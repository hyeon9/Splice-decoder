#!/usr/bin/env python3
# %%
from pdf2image import convert_from_path
import os
import sys
import argparse

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
    Description
    #########################################################
    This script makes summary stats of mapping process and merged input file.
    OUTPUT: mapping.stats, mat_tx_numbers.pdf, mapping_rate.pdf, fig_input.txt
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


def pdfs_to_html_with_titles(pdf_paths, titles, output_dir, txt_file, output_html):
    """
    Args:
        pdf_paths (list): PDF file path
        titles (list): title for each pdf
        output_html (str): out path
    """
    if len(pdf_paths) != len(titles):
        raise ValueError("pdf_paths should have same length of titles)

    # Make html 
    with open(f"{output_dir}{output_html}", 'w') as f:
        f.write('<html><body>\n')
        target_size = (400, 400)
        ## Main box >> it has several divs
        f.write('<div style="display: flex; align-items: center; margin-bottom: 20px;">\n')


        # PDF processing
        for idx, pdf_path in enumerate(pdf_paths[:2]):
            images = convert_from_path(pdf_path)
            image = images[0]
            image.thumbnail(target_size)
            img_save = f"{output_dir}/pdf_{idx+1}_page_{idx+1}.png"
            image.save(img_save, 'PNG')
            img_filename = f"pdf_{idx+1}_page_{idx+1}.png"
            f.write('<div style="text-align: center; margin-right: 40px;">\n')
            f.write(f'<h2>{titles[idx]}</h2>\n')
            f.write(f'<img src="{img_filename}" alt="PDF {idx+1} page {idx+1}"><br>\n')
            f.write('</div>\n')

        # Add summary table
        f.write('<div style="text-align: left; margin-left: 40px; margin-right: 40px;">\n') 
        f.write('<h2>Mapping Summary</h2>\n')
        f.write('<table border="1" style="width: 100%; margin: 0;">\n')

        with open(txt_file, 'r') as file:
            lines = file.readlines()

        for line in lines:
            f.write("      <tr>\n")
            columns = line.strip().split()

            for column in columns:
                f.write(f"        <td>{column}</td>\n")
            f.write("      </tr>\n")

        f.write('</table>\n')
        f.write('</div>\n')
        
        # 2nd Flexbox
        for idx, pdf_path in enumerate(pdf_paths[2:]):
            images = convert_from_path(pdf_path)
            image = images[0]
            image.thumbnail(target_size)
            img_save = f"{output_dir}/pdf_{idx+3}_page_{idx+3}.png"
            image.save(img_save, 'PNG')
            img_filename = f"pdf_{idx+3}_page_{idx+3}.png"
            f.write('<div style="text-align: center; margin-right: 40px;">\n')
            f.write(f'<h2>{titles[idx+2]}</h2>\n')
            f.write(f'<img src="{img_filename}" alt="PDF {idx+3} page {idx+3}"><br>\n')
            f.write('</div>\n')

        f.write('</div>\n')
        f.write('</body></html>\n')

    print(f"HTML file created: {output_dir}{output_html}")


args, parser = parse_args(sys.argv[1:])    
input = args.input
dir = f"{input}"
pdf_paths = [f'{dir}/figure/mapping_rate.pdf', f'{dir}/figure/mat_tx_numbers.pdf',\
             f'{dir}/figure/splicing_categories_stacked_plot.pdf', f'{dir}/figure/merged_stacked_plot.pdf']
titles = ['Mapping results', 'Average number of Reference TX', "Distribution of Functional DS", "Merged Distribution of Functional DS"]
pdfs_to_html_with_titles(pdf_paths, titles, f'{dir}/figure/', f'{dir}/mapping.stats', 'Summary.html')
# %%
