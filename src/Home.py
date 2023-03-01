# Description: This file is the entry point of streamlit for the application.
# It is called by the main.py file.
import shutil
import streamlit as st

from scripts.blast_downloader import BlastDownloader
from scripts.utils import *


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='🧬')

    st.title('BlastUI')
    st.write('This app will let you blast a query against your own custom sequences.')
    st.write('Here you can set up some options.')

    ##### 1) Install BLAST #####
    st.subheader("1) Install BLAST")
    st.write('Blast executables are required to run this app. If you already have them installed, you can skip '
             'this step. You can check below if the app was able to find them or not. '
             'In any case, you can download the latest version by clicking the corresponding button. '
             'If you do, the executable will be saved inside the Binaries folder of this app. ')

    col_path, col_binaries = st.columns([1, 1])
    with col_path:
        st.write('Blast executables found in this PC:')
        for program, path in check_blast_executables():
            if path:
                st.write(f'✅ **{program}**: *{path}*')
            else:
                st.write(f'❌ **{program}**')

    with col_binaries:
        st.write("Blast executables found in this app's folder:")
        for program, path in check_blast_executables_in_bin():
            if path:
                st.write(f'✅ **{program}**: *{path}*')
            else:
                st.write(f'❌ **{program}**')

    col_btn, col_checkbox, *cols = st.columns([1, 2, 3])

    # Make the force download appear only if the download button was clicked and the files already exists
    force_download = False
    if 'force_blast_download' in st.session_state:
        force_download = st.session_state['force_blast_download']

    if col_btn.button('Download blast'):
        pbar = st.progress(0)

        try:
            BlastDownloader(force_download=force_download, pbar=pbar)
            del st.session_state['force_blast_download']

        except shutil.Error as e:
            if 'already exists' in str(e):
                st.warning('A bin folder containing blast already exists!')

                # Initialize the session state variable, to show the checkbox
                st.session_state['force_blast_download'] = False
            else:
                raise e

        except OSError as e:
            if 'already downloaded' in str(e):
                st.warning('Blast was already downloaded! ')

                # Initialize the session state variable, to show the checkbox
                st.session_state['force_blast_download'] = False
            else:
                raise e

    if 'force_blast_download' in st.session_state:
        st.session_state['force_blast_download'] = col_checkbox.checkbox('Force download and '
                                                                         'overwrite existing files')

    # Choose which blast binaries to use and save it in the config file
    configs = read_configs()
    match configs['BLAST']['use_executables_in']:
        case 'BlastUI':
            i = 0
        case 'PATH':
            i = 1
        case _:
            i = 0

    match st.radio('Use the blast binaries:', ('Downloaded with BlastUI', 'In this PC'), index=i):
        case 'Downloaded with BlastUI':
            configs['BLAST']['use_executables_in'] = 'BlastUI'
        case 'In this PC':
            configs['BLAST']['use_executables_in'] = 'PATH'
    write_configs(configs)

    ##### 2) Choose the database #####
    st.subheader('2) Choose the database')
    st.write('The next step is to create you genome database from the fasta files. You can manage your'
             'databases by going to the corresponding page.')


if __name__ == '__main__':
    main()
