# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 17:04:18 2022

@author: Abdulkadir

98 genom analiz tablosu temel başlıklar: 
    coverage,               genome_results.txt (mean coverageData)          +++
    fastq file size,        fastq_size.txt
    duplication rate,       fastp.json                                      +++
    read count,             fastp.json (before filtering total reads)       +++
    filtered read count,    fastp.json (after filtering total reads)        +++
    read mean length,       qualimap.pdf (read mean length)                 ++        
    total bases,            fastp.json (before and after filtering )        +++
    q20,                    fastp.json (before and after filtering?)        +++
    q30,                    fastp.json (before and after filtering?)        +++
    variant count,          variant_count                                   +++
    elapsed time,           ??
    ..

"""
# TO-DO:
    # Parametrelere ulaş.

import os
import json
import pandas as pd

path = "C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastp"
dir_list = os.listdir(path)

print(dir_list)


fastq_size_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastq_size.txt", "rt")
variant_count_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//variant_count.txt", "rt")

for i in dir_list:
    genome_results_txt = open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{i}/genome_results.txt", "rt")
    
    lines_gr = genome_results_txt.readlines()
    lines_fs = fastq_size_txt.readlines()
    lines_vc = variant_count_txt.readlines()
    
    filtered_read_count = lines_gr[19][22:37]
    mean_coverage = lines_gr[71][24:33] # space leri cikar.
    genome_results_txt.close()

#%%
import json
import os

path = "C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastp"
dir_list = os.listdir(path)

fastq_size = []

fastq_size_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastq_size.txt", "rt")
lines_fastq_sizes = fastq_size_txt.readlines()
for line in lines_fastq_sizes:
    if ('SAA11A2' in line) and ("R1" in line):
        fastq_size.append(line + '\n')
        print(fastq_size)
            
        
        

fastq_size_txt.close()

#%%

fastq_size_txt.close()
variant_count_txt.close()

#%%

import numpy as np
import pandas as pd
import os
import json
import PyPDF2 

fastq_size_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastq_size.txt", "rt")
# variant_count_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//variant_count.txt", "rt")
# genome_results_txt = open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{i}/genome_results.txt", "rt")
# fastp_json = open(f"C:/Users/Abdulkadir/Desktop/Genome_Project_Documents/haplo/fastp/{j}/fastp.json")    

path = "C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastp"
IDs = os.listdir(path)

Samples = [] 
Coverages = [] 
Read_count = [] 
Filtered_read_count = [] 
Total_base = [] 
Filtered_total_base = [] 
q20 = [] 
Filtered_q20 = [] 
q30 = []
Filtered_q30 = [] 
Duplication = [] 
Variants = [] 
Read_mean_length = [] 



for ID in IDs:
    
    genome_results_txt = open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{ID}/genome_results.txt", "rt")
    
    lines_gen_res = genome_results_txt.readlines()
    
    mean_coverage = lines_gen_res[71][24:32]

    file = open(f"C:/Users/Abdulkadir/Desktop/Genome_Project_Documents/haplo/fastp/{ID}/fastp.json")
    fastp_json = json.load(file)
    
    before_read_count = fastp_json["summary"]["before_filtering"]["total_reads"]
    before_total_bases = fastp_json["summary"]["before_filtering"]["total_bases"]
    before_q20 = fastp_json["summary"]["before_filtering"]["q20_bases"]
    before_q30 = fastp_json["summary"]["before_filtering"]["q30_bases"]
    after_read_count = fastp_json["summary"]["after_filtering"]["total_reads"]
    after_total_bases = fastp_json["summary"]["after_filtering"]["total_bases"]
    after_q20 = fastp_json["summary"]["after_filtering"]["q20_bases"]
    after_q30 = fastp_json["summary"]["after_filtering"]["q30_bases"]
    duplication = fastp_json["duplication"]["rate"]
    
    qualimap_txt = open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{ID}/qualimap_{ID}.txt", "rt")
    
    lines_qualimap = qualimap_txt.readlines()
    
    read_mean_length = lines_qualimap[27][11:17]
    print(read_mean_length)
    
    Samples.append(ID)
    Coverages.append(float(mean_coverage))
    Read_count.append(before_read_count)
    Filtered_read_count.append(after_read_count)
    Total_base.append(before_total_bases)
    Filtered_total_base.append(after_total_bases)
    q20.append(before_q20)
    Filtered_q20.append(after_q20)
    q30.append(before_q30)
    Filtered_q30.append(after_q30)
    Duplication.append(100 * duplication)
    Read_mean_length.append(float(read_mean_length.replace(",",".")))
    


variant_count_txt = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//variant_count.txt", "rt")

lines_var_count = variant_count_txt.readlines()

for i in lines_var_count:
    variant = (i[-8:-1]) # istediğim degerler i nin icinde.
    Variants.append(int(variant))
    

# fastq_size_txt.close()
variant_count_txt.close()
genome_results_txt.close()
file.close()
qualimap_txt.close()

master_table = pd.DataFrame()
master_table['Samples'] = Samples
master_table['Coverages'] = Coverages
master_table['Variants'] = Variants
master_table['Read Mean Length'] = Read_mean_length
master_table['Read Count'] = Read_count
master_table['Filtered Read Count'] = Filtered_read_count
master_table['Total Base'] = Total_base
master_table['Filtered Total Base'] = Filtered_total_base
master_table['q20'] = q20
master_table['Filtered q20'] = Filtered_q20
master_table['q30'] = q30
master_table['Filtered q30'] = Filtered_q30
master_table['Duplication'] = Duplication
master_table.set_index('Samples', inplace=True)

master_table.to_csv('C://Users//Abdulkadir//Desktop//Genome_Project_Documents//master_table.csv')

#%% PDF'leri .TXT dosyasına dönüştürdük.

import PyPDF2

for ID in IDs:    
    qualimap_pdf = open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{ID}/qualimap_{ID}.pdf", "rb")
    qualimap_pdf_read = PyPDF2.PdfFileReader(qualimap_pdf)
    pageObj = qualimap_pdf_read.getPage(2)
    text=pageObj.extractText()
    file1=open(f"C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_{ID}//qualimap_{ID}.txt","a")
    file1.writelines(text)
#%%

import PyPDF2
qualimap_pdf = open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_SAA11T4/qualimap_SAA11T4.pdf", "rb")
qualimap_pdf_read = PyPDF2.PdfFileReader(qualimap_pdf)
pageObj = qualimap_pdf_read.getPage(2)
text=pageObj.extractText()
file1=open("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//bam_qc//qualimap_SAA11T4/qualimap_SAA11T4.txt","a")
file1.writelines(text)

#%%
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

filtered_master_table = master_table[master_table['Coverages'] > 30]

# Histograms

Variant_his = master_table.hist(column = 'Variants')
plt.xlabel("Variant sayıları (milyon)")
plt.ylabel("Sayım")
plt.title("Varyant Dağılımı")


Filtered_Variant_his = filtered_master_table.hist(column = 'Variants')
plt.xlabel("Variant sayıları (milyon)")
plt.ylabel("Sayım")
plt.title("Filtrelenmiş veride varyant dağılımı")


Coverages_his = master_table.hist(column = 'Coverages')
plt.xlabel("Coverage değerleri")
plt.ylabel("Sayım")
plt.title("Coverage dağılımı")


Filtered_Coverages_his = filtered_master_table.hist(column = 'Coverages')
plt.xlabel("Coverage değerleri")
plt.ylabel("Sayım")
plt.title("Filtrelenmiş veride coverage dağılımı")


Coverages_his = master_table.hist(column = 'Duplication')
plt.xlabel("Duplication değerleri (%)")
plt.ylabel("Sayım")
plt.title("Duplication dağılımı")


Filtered_Coverages_his = filtered_master_table.hist(column = 'Duplication')
plt.xlabel("Duplication değerleri (%)")
plt.ylabel("Sayım")
plt.title("Filtrelenmiş veride duplication dağılımı")


Coverages_his = master_table.hist(column = 'Read Count')
plt.xlabel("Okuma sayısı")
plt.ylabel("Sayım")
plt.title("Okuma sayıları dağılımı")


Filtered_Coverages_his = filtered_master_table.hist(column = 'Read Count')
plt.xlabel("Okuma sayısı değerleri")
plt.ylabel("Sayım")
plt.title("Filtrelenmiş veride okuma sayıları dağılımı")


# Correlation Matrix

cormat = master_table.corr()
round(cormat,2)
#cormat.to_csv('C://Users//Abdulkadir//Desktop//Genome_Project_Documents//cormat.csv')

filtered_cormat = filtered_master_table.corr()
round(filtered_cormat,2)
#filtered_cormat.to_csv('C://Users//Abdulkadir//Desktop//Genome_Project_Documents//filtered_cormat.csv')

# Scatterplots

sns.lmplot(x='Coverages', y='Duplication', data=master_table)
sns.lmplot(x='Coverages', y='Duplication', data=filtered_master_table)

sns.lmplot(x='Coverages', y='Read Count', data = master_table)
sns.lmplot(x='Coverages', y='Read Count', data = filtered_master_table)

sns.lmplot(x='Coverages', y='Variants', data=master_table)
sns.lmplot(x='Coverages', y='Variants', data=filtered_master_table)

sns.lmplot(x='Variants', y='Read Count', data=master_table)
sns.lmplot(x='Variants', y='Read Count', data=filtered_master_table)


#%%
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

master_table_pd = pd.read_csv("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//master_table_fastq.csv")

cormat = master_table_pd.corr()
round(cormat,2)

filtered_master_table_pd = master_table_pd[master_table_pd['Coverages'] > 30]

Filtered_Total_fastq_his = filtered_master_table_pd.hist(column = 'Total_fastq_size')
plt.xlabel("Fastq boyutları")
plt.ylabel("Sayım")
plt.title("Filtrelenmiş veride fastq boyut dağılımı")


cormat.to_csv('C://Users//Abdulkadir//Desktop//Genome_Project_Documents//filtered_cormat_pd.csv')


sns.lmplot(x='Coverages', y='Total_fastq_size', data=filtered_master_table_pd)

