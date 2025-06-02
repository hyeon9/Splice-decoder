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
    여러 개의 PDF 파일을 각 페이지별 이미지로 변환하고,
    각 이미지에 알맞은 제목을 추가하여 HTML 파일 생성.

    Args:
        pdf_paths (list): PDF 파일 경로 리스트.
        titles (list): 각 PDF에 해당하는 제목 리스트.
        output_html (str): 생성될 HTML 파일 경로.
    """
    if len(pdf_paths) != len(titles):
        raise ValueError("PDF 경로 리스트와 제목 리스트의 길이가 같아야 합니다.")

    # HTML 파일 생성
    with open(f"{output_dir}{output_html}", 'w') as f:
        f.write('<html><body>\n')
        target_size = (400, 400)
        ## Main box >> it has several divs
        f.write('<div style="display: flex; align-items: center; margin-bottom: 20px;">\n')


        # 각 PDF 파일에 대해 처리
        for idx, pdf_path in enumerate(pdf_paths[:2]):
            images = convert_from_path(pdf_path)  # PDF 파일을 이미지로 변환
            image = images[0]
            image.thumbnail(target_size)
            img_save = f"{output_dir}/pdf_{idx+1}_page_{idx+1}.png"
            image.save(img_save, 'PNG')
            img_filename = f"pdf_{idx+1}_page_{idx+1}.png"
            f.write('<div style="text-align: center; margin-right: 40px;">\n')  # 마진 추가
            f.write(f'<h2>{titles[idx]}</h2>\n')  # 제목 추가
            f.write(f'<img src="{img_filename}" alt="PDF {idx+1} page {idx+1}"><br>\n')
            f.write('</div>\n')  # 첫 번째 이미지 div 끝

        # 테이블 추가
        f.write('<div style="text-align: left; margin-left: 40px; margin-right: 40px;">\n')  # 테이블 div
        f.write('<h2>Mapping Summary</h2>\n')
        f.write('<table border="1" style="width: 100%; margin: 0;">\n')  # 테이블 스타일

        with open(txt_file, 'r') as file:
            # 텍스트 파일을 줄 단위로 읽기
            lines = file.readlines()

        for line in lines:
            f.write("      <tr>\n")  # 새로운 행 생성
            # 탭 또는 콤마로 구분된 데이터를 셀로 변환
            columns = line.strip().split()  # 기본 구분자는 공백

            for column in columns:
                f.write(f"        <td>{column}</td>\n")  # 각 데이터를 셀로 변환
            f.write("      </tr>\n")  # 행 종료

        f.write('</table>\n')  # 테이블 닫기
        f.write('</div>\n')  # 테이블 div 끝
        # f.write('</div>\n')  # 첫번째 Flexbox div 끝
        
        # ## 두번째 Flexbox
        # f.write('<div style="display: flex; align-items: center; margin-bottom: 20px;">\n')
        for idx, pdf_path in enumerate(pdf_paths[2:]):
            images = convert_from_path(pdf_path)  # PDF 파일을 이미지로 변환
            image = images[0]
            image.thumbnail(target_size)
            img_save = f"{output_dir}/pdf_{idx+3}_page_{idx+3}.png"
            image.save(img_save, 'PNG')
            img_filename = f"pdf_{idx+3}_page_{idx+3}.png"
            f.write('<div style="text-align: center; margin-right: 40px;">\n')  # 마진 추가
            f.write(f'<h2>{titles[idx+2]}</h2>\n')  # 제목 추가
            f.write(f'<img src="{img_filename}" alt="PDF {idx+3} page {idx+3}"><br>\n')
            f.write('</div>\n')  # 첫 번째 이미지 div 끝

        f.write('</div>\n')  # 두번째 Flexbox div 끝
        f.write('</body></html>\n')

    print(f"HTML file created: {output_dir}{output_html}")


# 사용 예시
args, parser = parse_args(sys.argv[1:])    
input = args.input
dir = f"{input}"
dir2 = f"{input}/result/"
pdf_paths = [f'{dir}/figure/mapping_rate.pdf', f'{dir}/figure/mat_tx_numbers.pdf',\
             f'{dir2}/figure/splicing_categories_stacked_plot.pdf', f'{dir2}/figure/merged_stacked_plot.pdf']  # PDF 파일 경로 리스트
titles = ['Mapping results', 'Average number of Reference TX', "Distribution of Functional DS", "Merged Distribution of Functional DS"]  # 각 PDF 파일의 제목 리스트
pdfs_to_html_with_titles(pdf_paths, titles, f'{dir}/figure/', f'{dir}/mapping.stats', 'Summary.html')
# %%
