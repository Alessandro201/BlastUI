# Description: This file is the entry point of the streamlit application
# It is called by run_app.py

import sys
from io import BytesIO
from pathlib import Path
import shutil

# Needed to search for scripts in the parent folder
sys.path.append(str(Path(__file__).parent))

import streamlit as st
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_downloader import BlastDownloader
from scripts import utils
import ftputil.error


def set_not_check_blast():
    st.session_state['check_blast'] = False


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=BytesIO(Path(utils.resource_path('.'), 'icon.png').read_bytes()))

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

    # This variable is needed to avoid checking if blast is installed after the user clicks the
    # 'Go to the database page' button. This is due to how streamlit works as if the user clicks
    # on the button the page will be reloaded and this script executed again before actually doing
    # what comes in the if statement after the button.
    if 'check_blast' not in st.session_state:
        st.session_state['check_blast'] = True

    if 'blast_exec' not in st.session_state:
        st.session_state['blast_exec'] = None

    if 'blast_exec_version' not in st.session_state:
        st.session_state['blast_exec_version'] = None

    if st.session_state['check_blast']:
        with st.spinner('Checking if blast is installed...'):
            blast_exec: dict | None = utils.get_programs_path()
            blast_version: str | None = None

            if not blast_exec:
                st.warning(f'Blast was not found in your computer. '
                           f'Click the download button to get the latest version.')

            else:
                blast_version = utils.check_blast_version(blast_exec['blastn'])
                st.session_state['blast_exec_version'] = blast_version
                st.session_state['blast_exec'] = blast_exec

        if blast_version:
            st.success(f'Using BLAST {blast_version} but if you want you can '
                       f'click the download button to get the latest version.')

    if st.button('Download blast', on_click=set_not_check_blast):
        pbar = st.progress(0)

        download_folder = Path('./Binaries/')
        shutil.rmtree(download_folder, ignore_errors=True)
        download_folder.mkdir(parents=True, exist_ok=True)

        try:
            BlastDownloader(pbar=pbar)
        except ValueError as e:
            st.error(e)
            st.stop()

        blast_exec = utils.get_programs_path()
        blast_version = utils.check_blast_version(blast_exec['blastn'])
        st.success(f'Blast {blast_version} was downloaded successfully!')

        st.session_state['blast_exec_version'] = blast_version
        st.session_state['blast_exec'] = blast_exec

    # Reset the check_blast variable so that the next time this page is loaded it will check if blast is installed
    # unless the user clicks the 'Go to the database page' button
    st.session_state['check_blast'] = True

    if st.session_state['blast_exec'] is not None:
        ##### 2) Choose the database #####
        st.header('Choose the database')
        st.write('The next step is to create your genome database from the fasta files. You can manage your '
                 'databases by going to the corresponding page.')

        if st.button('Go to "Manage Databases" page', on_click=set_not_check_blast):
            switch_page('Manage Databases')


if __name__ == '__main__':
    main()
