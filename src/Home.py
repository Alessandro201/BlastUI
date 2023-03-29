# Description: This file is the entry point of the streamlit application
# It is called by run_app.py

import sys
from pathlib import Path

# Needed to search for scripts in the parent folder
sys.path.append(str(Path(__file__).parent))

import streamlit as st
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_downloader import BlastDownloader, DownloadError
from scripts.utils import *
import ftplib


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=Path(resource_path('.'), 'icon.png').read_bytes())

    st.title('BlastUI')
    st.markdown("""
    From the [blast website](https://blast.ncbi.nlm.nih.gov/Blast.cgi): 
    > The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. 
    The program compares nucleotide or protein sequences to sequence databases and calculates the 
    statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships 
    between sequences as well as help identify members of gene families. 

    This WebApp is a simple interface to run blast locally against your own sequences. You can then filter the 
    results and save them in a file.
    """)

    ##### 1) Install BLAST #####
    st.markdown("""
    ## Install BLAST
    Here you can download the latest version of BLAST but if you have already installed it you can skip this step.
    You can check below if the app was able to find it or not.
    """)

    with st.spinner('Checking if blast is installed...'):
        blast_exec = get_programs_path()

        if not blast_exec:
            st.warning(f'Blast was not found in your computer. Click the download button to get the latest version.')

        if blast_exec['blastn'].parent == Path('Binaries/bin'):
            blast_version = check_blast_version(blast_exec['blastn'])
            st.success(f'Blast {blast_version} was found in your computer, but if you want you can click the download '
                       f'button to get the latest version.')

        else:
            blast_version = check_blast_version(blast_exec['blastn'])
            st.success(f'Blast {blast_version} was found in your computer, but if you want you can click the download '
                       f'button to get the latest version.')

    if st.button('Download blast'):
        pbar = st.progress(0)

        download_folder = Path('./Binaries/')
        shutil.rmtree(download_folder, ignore_errors=True)
        download_folder.mkdir(parents=True, exist_ok=True)

        try:
            BlastDownloader(pbar=pbar)
        except DownloadError:
            st.error('The downloaded BLAST was corrupted (hashes do not match). Try again.')
            st.stop()
        except ftplib.all_errors:
            st.error('There was an error during the download. Try again.')

        blast_exec = get_programs_path()
        blast_version = check_blast_version(blast_exec['blastn'])
        st.success(f'Blast {blast_version} was downloaded successfully!')

    st.session_state['blast_exec'] = blast_exec

    if st.session_state['blast_exec'] is not None:
        ##### 2) Choose the database #####
        st.header('Choose the database')
        st.write('The next step is to create your genome database from the fasta files. You can manage your '
                 'databases by going to the corresponding page.')
        if st.button('Go to the database page'):
            switch_page('Manage Genome Databases')


if __name__ == '__main__':
    main()
