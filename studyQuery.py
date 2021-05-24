""" Query metabolomics workbench for raw MS data

    Fulgens Consulting, LLC """

import shutil, os, zipfile
import urllib.request as request
from contextlib import closing

studyID = 'ST000934'
baseURL = 'ftp://www.metabolomicsworkbench.org/Studies/' 

def main(studyID, baseURL):
    if studyID+'.zip' not in os.listdir('data'):

        fout = 'data/{}'.format(studyID+'.zip')
        with closing(request.urlopen(baseURL+studyID+'.zip')) as r:
            with open(fout, 'wb') as f:
                shutil.copyfileobj(r, f)

        with zipfile.ZipFile(fout, 'r') as zip_ref:
            zip_ref.extractall('data/{}'.format(studyID))

if __name__ == '__main__':
    main(studyID, baseURL)

    