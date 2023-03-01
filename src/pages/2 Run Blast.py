from multiprocessing import cpu_count

import base64
import shlex
import pandas as pd
import shutil
import subprocess
import sys
import base64
from datetime import datetime
from math import ceil
from pathlib import PurePath

import streamlit as st
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid, DataReturnMode, ColumnsAutoSizeMode, GridUpdateMode, JsCode
from streamlit_option_menu import option_menu

from scripts.blast_response import *
from scripts.utils import *

from timebudget import timebudget


@st.cache_data(show_spinner=False)
def blast(query: str, blast_mode: str, db: str, threads: int, max_target_seq: int) -> BlastResponse:
    """
    This function runs the blast command and returns the results as a dictionary.
    """

    # Write query to file to be used by blast
    query_file = './Analysis/query.fasta'
    Path(query_file).write_text(query)

    today = datetime.today()
    today_time = today.strftime("%Y%m%d_%H%M%S")
    out_file = Path(f'./Analysis/{today_time}_results.json')

    exec_in = read_configs()['BLAST']['use_executables_in']

    match blast_mode:
        case 'tblastn':
            tblastn_path = get_program_path('tblastn', binaries_in=exec_in)

            cmd = shlex.split(f'"{tblastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                              f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 '
                              f'-max_target_seqs {max_target_seq}')

        case 'blastn':
            blastn_path = get_program_path('blastn', binaries_in=exec_in)
            cmd = shlex.split(f'"{blastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                              f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 '
                              f'-max_target_seqs {max_target_seq}')

        case 'blastp':
            blastn_path = get_program_path('blastp', binaries_in=exec_in)
            cmd = shlex.split(f'"{blastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                              f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 '
                              f'-max_target_seqs {max_target_seq}')

        case 'blastx':
            blastn_path = get_program_path('blastx', binaries_in=exec_in)
            cmd = shlex.split(f'"{blastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                              f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 '
                              f'-max_target_seqs {max_target_seq}')

        case 'tblastx':
            blastn_path = get_program_path('tblastx', binaries_in=exec_in)
            cmd = shlex.split(f'"{blastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                              f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 '
                              f'-max_target_seqs {max_target_seq}')

        case _:
            raise ValueError(f'blast_mode {blast_mode} not recognized')

    try:
        run_command(cmd)
    except CalledProcessError as e:
        stderr = e.stderr
        st.error(f'Error running blast: {stderr}')

        if 'BLAST Database error: No alias or index file found for nucleotide database' in stderr:
            st.info(f'It seems you were trying to do a ***{blast_mode.upper()}*** which requires a '
                    f'nucleotide database, but ***{Path(db).parent.name}*** is of proteins.')
        elif 'BLAST Database error: No alias or index file found for protein database' in stderr:
            st.info(f'It seems you were trying to do a ***{blast_mode.upper()}*** which requires a '
                    f'protein database, but ***{Path(db).parent.name}*** is of nucleotides.')
        else:
            raise e

        st.stop()

    return BlastResponse(out_file)


def choose_database(container=None):
    if not container:
        container = st

    Path('./BlastDatabases').mkdir(parents=True, exist_ok=True)
    dbs = [path.name for path in Path('./BlastDatabases').iterdir() if path.is_dir()]

    if dbs:

        previous_db = st.session_state.get('db', None)
        previous_db_index = dbs.index(previous_db.parent.name) if previous_db else 0

        db = container.selectbox('Select Blast Database', dbs, index=previous_db_index)
        st.session_state.db = Path('./BlastDatabases', db, 'blastdb')
    else:
        container.warning('No databases found. Please make one in the Manage Genome Databases section.')


def sidebar_options():
    """
    This function creates the sidebar with the options and returns a dictionary of them.
    """

    st.sidebar.title("Options")

    st.session_state['threads'] = st.sidebar.number_input('Threads to use: ',
                                                          min_value=1, max_value=cpu_count(),
                                                          value=round(cpu_count() / 2), step=1)

    st.sidebar.subheader('Utilities')

    if st.sidebar.button("Clear Cache"):
        # Clear values from *all* all in-memory and on-disk data caches:
        st.cache_data.clear()
        st.sidebar.info('Cache cleared')

    if st.sidebar.button("Clear session state"):
        for key in st.session_state.keys():
            del st.session_state[key]

    if st.sidebar.button("Clear analysis folder"):
        st.sidebar.warning('All the analysis files will be deleted!')

        def clear_analysis_folder():
            for file in Path('./Analysis').iterdir():
                file.unlink()
            st.session_state['analysis_folder_cleared'] = True

        st.sidebar.button('Confirm', on_click=clear_analysis_folder)

    if 'analysis_folder_cleared' in st.session_state:
        st.sidebar.success('Analysis folder cleared!')
        del st.session_state['analysis_folder_cleared']


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='🧬')

    sidebar_options()
    st.title('Blast queries against your local database!')

    ###### Query and blast mode ######
    query = st.text_area('Insert the query you want to blast: ', placeholder="Query...", height=210).strip()
    st.session_state.query = query

    st.session_state.blast_spinner_container = st.empty()

    ###### Blast button and spinner ######
    choose_database()

    blast_modes = ["blastn", "blastp", "blastx", 'tblastn', 'tblastx']
    st.session_state.blast_mode = st.radio('Blast mode', options=blast_modes, index=3)

    st.session_state.max_target_seq = st.number_input('Max sequences per query: ', min_value=1,
                                                      value=500, step=100)

    if st.button('Blast query'):
        st.session_state.switch_to_result_page = False

        if 'db' not in st.session_state:
            st.warning('No databases found. Please make one in the Manage Genome Databases section.')
            st.stop()

        if len(st.session_state.get('query', '')) == 0:
            st.warning('Please enter a query!')
            st.stop()

        with st.session_state.blast_spinner_container:
            with st.spinner(f"Running {st.session_state.blast_mode}..."):
                blast_response = blast(query=st.session_state['query'],
                                       blast_mode=st.session_state['blast_mode'],
                                       db=st.session_state['db'],
                                       threads=st.session_state['threads'],
                                       max_target_seq=st.session_state['max_target_seq'])

        st.session_state['BlastResponse']: BlastResponse = blast_response
        st.session_state['df'] = blast_response.reindexed_df
        st.session_state['switch_to_result_page'] = True

    if st.session_state.get('switch_to_result_page', False):
        blast_response = st.session_state.BlastResponse
        st.session_state.switch_to_result_page = False

        if not blast_response.messages:
            switch_page('Results')

        for message in blast_response.messages:
            st.warning(message)

        if st.button('Go to results'):
            switch_page('Results')


if __name__ == "__main__":
    main()
