import pandas as pd
import os
from ftplib import FTP

# import study ids to downlaod
study_ids = open('MetDMN\other_studies.txt', 'r').read().splitlines()

print(study_ids)

# ftp download 
ftp_url = 'ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/'
ftp = FTP('ftp.ebi.ac.uk') 
ftp.login()

for study_id in study_ids:
    ftp.cwd('/pub/databases/metabolights/studies/public/' + study_id + '/')

    study_url = ftp_url + study_id
    sample_file_url = 's_' + study_id + '.txt'

    # can be multiple mafs
    maf_file_url = 'm_' + study_id + '*.txt'

    # make directory
    if not os.path.exists('MetDMN/Studies/' + study_id):
        os.makedirs('MetDMN/Studies/' + study_id)

    # download sample file
    sampfile = open('MetDMN/Studies/' + study_id + '/' + sample_file_url, 'wb')
    print(sample_file_url)
    ftp.retrbinary('RETR ' + sample_file_url, sampfile.write)
    sampfile.close()

    # download maf files
    maf_files = ftp.nlst('*maf.tsv')
    for file in maf_files:
        maf_file = open('MetDMN/Studies/' + study_id + '/' + file, 'wb')
        print(file)
        ftp.retrbinary('RETR ' + file, maf_file.write)
        maf_file.close()

ftp.quit()
