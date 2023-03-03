from multiprocessing import cpu_count

import base64
import shlex
import pandas as pd
import shutil
import subprocess
import sys
import base64
from io import StringIO, BytesIO
from datetime import datetime
from math import ceil
from pathlib import PurePath

import streamlit as st
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid, DataReturnMode, ColumnsAutoSizeMode, GridUpdateMode, JsCode

from streamlit_option_menu import option_menu
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_response import *
from scripts.utils import *

from timebudget import timebudget


@st.cache_data(show_spinner=False)
def blast(query: str, blast_mode: str, db: str, threads: int, **kwargs) -> BlastResponse:
    """
    This function runs the blast command and returns the results as a dictionary.
    """

    evalue = kwargs.get('evalue', None)
    matrix = kwargs.get('matrix', None)
    word_size = kwargs.get('word_size', None)
    gap_penalty = kwargs.get('gap_penalty', None)
    max_target_seq = kwargs.get('max_target_seq', None)

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


def read_query(uploaded_files: BytesIO) -> str:
    """
    This function reads the query from the uploaded files and returns it as a string.

    :param uploaded_files: files uploaded with st.file_uploader to be read
    :return: the uploaded files read and joined as a single string
    """

    queries = list()
    for uploaded_file in uploaded_files:
        # To convert to a string based IO
        query_io = StringIO(uploaded_file.getvalue().decode("utf-8"))
        query = query_io.read().strip()
        queries.append(query)

    return '\n'.join(queries)


def set_advanced_options(container=None):
    if container is None:
        container = st

    options = dict()

    row1 = container.container()
    row1_col1, row1_col2 = row1.columns([1, 1])

    with row1_col1:
        gap_penalty = ['Existence: 11 Extension: 2', 'Existence: 10 Extension: 2', 'Existence: 9 Extension: 2',
                       'Existence: 8 Extension: 2', 'Existence: 7 Extension: 2', 'Existence: 6 Extension: 2',
                       'Existence: 13 Extension: 1', 'Existence: 12 Extension: 1', 'Existence: 11 Extension: 1',
                       'Existence: 10 Extension: 1', 'Existence: 9 Extension: 1']

        options['evalue'] = st.select_slider('E-value: ', options=[10 ** i for i in range(-100, 4)], value=10)

    with row1_col2:
        matrixes = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']

        options['matrix'] = st.selectbox('Matrix: ', options=matrixes, index=matrixes.index('BLOSUM62'))

    row2 = container.container()
    row2_col1, row2_col2 = row2.columns([1, 1])

    with row2_col1:
        options['word_size'] = st.number_input('Word size: ', min_value=2, value=6, step=1)

        options['gap_penalty'] = st.selectbox('Gap penalty: ', options=gap_penalty,
                                              index=gap_penalty.index('Existence: 11 Extension: 1'))

    with row2_col2:
        options['max_target_seq'] = st.number_input('Max sequences per query: ', min_value=1, value=500,
                                                    step=100)

    st.session_state.advanced_options = options


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

    ###### BLAST MODE ######
    blast_modes = ["BLASTN", "BLASTP", "BLASTX", 'TBLASTN', 'TBLASTX']
    icons = ['list-task', 'list-task', "list-task", 'list-task', 'list-task']
    st.session_state.blast_mode = option_menu('', options=blast_modes, icons=icons,
                                              menu_icon="gear", default_index=3, orientation="horizontal").lower()

    ###### QUERY ######
    query = st.text_area('Insert the queries: ', placeholder="Query...", height=200).strip()
    st.session_state.query = query

    uploaded_files = st.file_uploader("Alternatively upload queries in fasta format", type=["fasta", "faa"],
                                      accept_multiple_files=True)
    if uploaded_files:
        if st.session_state.query:
            st.warning('You have uploaded files and written a query. The query will be ignored.')
        st.session_state.query = read_query(uploaded_files)

    ###### BLAST OPTIONS ######
    choose_database()

    exp = st.expander('Advanced options')
    set_advanced_options(exp)

    ###### BLAST ######
    if st.button('Blast query'):
        st.session_state.switch_to_result_page = False

        if 'db' not in st.session_state:
            st.warning('No databases found. Please make one in the Manage Genome Databases section.')
            st.stop()

        if len(st.session_state.get('query', '')) == 0:
            st.warning('Please enter a query!')
            st.stop()

        with st.spinner(f"Running {st.session_state.blast_mode}..."):
            blast_response = blast(query=st.session_state['query'],
                                   blast_mode=st.session_state['blast_mode'],
                                   db=st.session_state['db'],
                                   threads=st.session_state['threads'],
                                   **st.session_state['advanced_options'])

        st.session_state['BlastResponse']: BlastResponse = blast_response
        st.session_state['df'] = blast_response.df
        st.session_state['switch_to_result_page'] = True

    if st.session_state.get('switch_to_result_page', False):
        blast_response = st.session_state.BlastResponse

        if not blast_response.messages:
            st.session_state.switch_to_result_page = False
            switch_page('Results')

        for message in blast_response.messages:
            st.warning(message)

        if st.button('Go to results'):
            st.session_state.switch_to_result_page = False
            switch_page('Results')


if __name__ == "__main__":
    main()
