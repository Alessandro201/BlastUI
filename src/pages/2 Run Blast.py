import pandas as pd
import shlex
import shutil
import subprocess
from datetime import datetime
from math import ceil
from pathlib import PurePath
from multiprocessing import cpu_count
import sys
import base64

import streamlit as st
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode, ColumnsAutoSizeMode

from scripts.blast import *
from scripts.blast_response import *
from scripts.utils import *


def run_command(command, error_description=''):
    """
    This function is a wrapper to run a command in the console.
    It prints the output as soon as it's given, prints the error and raises a subprocess.CalledProcessError exception
    to keep track of where the error was found and which command returned it.

    """

    if not error_description:
        error_description = f'ERROR FOUND'

    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1,
                          universal_newlines=True) as popen:

        # Save the errors in a variable to print them outside the with statement in case the return code is not zero.
        # It needs to be this way because inside the with statement stderr works but no return code is given,
        # and outside stderr does not work but there is a return code
        stderr = ''.join(popen.stderr)

    if popen.returncode != 0:
        print(f'\n{error_description}: ')
        for line in stderr:
            print(line, end='')
        print('')
        raise subprocess.CalledProcessError(popen.returncode, popen.args)


def tblastn(query_file: Union[str, Path], db: str, threads: int, max_target_seq: int) -> BlastResponse:
    today = datetime.today()
    today_time = today.strftime("%Y%m%d_%H%M%S")

    out_file = Path(f'./Analysis/{today_time}_results.json')

    exec_in = read_configs()['BLAST']['use_executables_in']
    tblastn_path = get_program_path('tblastn', binaries_in=exec_in)

    cmd = shlex.split(f'"{tblastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                      f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 -max_target_seqs {max_target_seq}')

    run_command(cmd, error_description='Error running tblastn')
    return BlastResponse(out_file)


def blastn(query_file: Union[str, Path], db: str, threads: int, max_target_seq: int) -> BlastResponse:
    today = datetime.today()
    today_time = today.strftime("%Y%m%d_%H%M%S")

    out_file = Path(f'./Analysis/{today_time}_results.json')

    exec_in = read_configs()['BLAST']['use_executables_in']
    blastn_path = get_program_path('blastn', binaries_in=exec_in)

    cmd = shlex.split(f'"{blastn_path}" -query "{query_file}" -db "{db}" -outfmt 15 '
                      f'-out "{out_file}" -num_threads {threads} -qcov_hsp_perc 10 -max_target_seqs {max_target_seq}')

    run_command(cmd, error_description='Error running blastn')
    return BlastResponse(out_file)


@st.cache_data
def blast(query: str, blast_mode: str, db: str, threads: int, max_target_seq: int):
    """
    This function runs the blast command and returns the results as a dictionary.
    """

    query_file = './Analysis/query.fasta'
    save_query(query, query_file)

    if blast_mode == 'tblastn':
        blast_response = tblastn(query_file, db, threads, max_target_seq)
    elif blast_mode == 'blastn':
        blast_response = blastn(query_file, db, threads, max_target_seq)
    else:
        raise ValueError(f'blast_mode {blast_mode} not recognized')

    return blast_response


@st.cache_data
def generate_df(data: dict, blast_mode: str) -> pd.DataFrame:
    """
    This function generate a pandas.DataFrame from the results of a blast analysis.

    :param data: Dictionary with the results of the blast analysis
                 run with parameter "-outfmt 15" and read by read_json()
    :param blast_mode: String with the blast mode used to run the analysis:
        - 'tblastn'
        - 'blastn'
    :return: pandas.DataFrame with the results of the blast analysis
    """

    if blast_mode == 'tblastn':
        df: pd.DataFrame = generate_tblastn_df(data)
    elif blast_mode == 'blastn':
        df: pd.DataFrame = generate_blastn_df(data)
    else:
        raise ValueError(f'blast_mode {blast_mode} not recognized')

    return df


@st.cache_data
def generate_alignments(data: list, indexes: list, blast_mode: str) -> list:
    alignments = list()

    for match in get_matches_by_index(data, indexes):
        if blast_mode == 'tblastn':
            alignment = generate_tblastn_alignment(match)
        elif blast_mode == 'blastn':
            alignment = generate_blastn_alignment(match)
        else:
            raise ValueError(f'blast_mode {blast_mode} not recognized')

        alignments.append((match["strain"], match["node"], alignment))

    return alignments


@st.cache_data
def convert_df_to_tsv(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    # return df.drop(['index']).to_csv(sep='\t', index_col=False).encode('utf-8')
    return df.to_csv(sep='\t', index=False).encode('utf-8')


def save_query(query, query_file):
    # Write query to file because tblastn does not accept it as a string
    Path(query_file).write_text(query)


@st.cache_data
def extract_indexes(selected: list) -> list:
    indexes = []
    for row in selected:
        nodeRowIndex = row['_selectedRowNodeInfo']['nodeRowIndex']
        match_index = row['id']
        indexes.append((nodeRowIndex, match_index))

    indexes = sorted(indexes, key=lambda x: x[0])
    return [index[1] for index in indexes]


def _download_button(object_to_download, download_filename):
    """
    Generates a link to download the given object_to_download.
    Params:
    ------
    object_to_download:  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv,
    Returns:
    -------
    (str): the anchor tag to download object_to_download
    """

    ext = PurePath(download_filename).suffix

    if isinstance(object_to_download, pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)
    elif isinstance(object_to_download, str):
        pass
    # Try JSON encode for everything else
    else:
        object_to_download = json.dumps(object_to_download)

    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    dl_link = f"""
    <html>
    <head>
    <title>Start Auto Download file</title>
    <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
    <script>
    $('<a href="data:text/{ext};base64,{b64}" download="{download_filename}">')[0].click()
    </script>
    </head>
    </html>
    """
    return dl_link


def download_all_alignments():
    st.info('Generating all the alignments...')
    blast_response = st.session_state['BlastResponse']
    grid_df = st.session_state['grid_df']
    indexes = grid_df['id']

    alignments = blast_response.alignments(indexes=indexes)
    alignments = '\n'.join(alignments)

    filename = "alignments.txt"
    components.html(
        _download_button(alignments, filename),
        height=None,
    )


def set_download_buttons(container=None):
    if not container:
        container = st

    col1, col2, col3, col4 = container.columns([1, 1, 1, 1])

    # Downloads the whole table
    with col1:
        df = st.session_state['grid_df']
        col1.download_button(
            label="Download table",
            data=convert_df_to_tsv(df),
            file_name='blast.tsv',
            mime='text/tsv')

    # Downloads only selected rows
    with col2:
        selected = st.session_state['selected']

        # Check that it's not empty
        selected_data = ''
        if selected:
            selected_df = pd.DataFrame(selected)
            selected_df.drop(columns=['_selectedRowNodeInfo', 'id'], inplace=True)
            selected_data = convert_df_to_tsv(selected_df)

        col2.download_button(
            label="Download selected rows",
            data=selected_data,
            file_name='blast_rows.tsv',
            mime='text/tsv')

    # Download all alignments
    with col3:
        col3.button('Download all alignments', on_click=download_all_alignments)

    # Download only selected alignments
    with col4:
        # Row alignments download button
        row_alignments = ''
        if 'row_alignments' in st.session_state:
            row_alignments = st.session_state['row_alignments']

        col4.download_button(
            label="Download selected alignments",
            data=row_alignments,
            file_name='alignments_rows.text',
            mime='text')


@st.cache_data
def get_aggrid_options(df: pd.DataFrame) -> dict:
    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_column('id', hide=True)

    if df.sseq.str.len().max() >= 1000:
        gb.configure_column('sseq', hide=True)

    gb.configure_column('sseq', hide=True)

    gb.configure_default_column(groupable=True,
                                value=True,
                                enableRowGroup=True,
                                aggFunc='avg',
                                editable=False)

    gb.configure_column('query_title', aggFunc='count')
    gb.configure_column('strain', aggFunc='count')
    gb.configure_column('node', aggFunc='count')

    gb.configure_selection('multiple',
                           use_checkbox=False,
                           groupSelectsChildren=True,
                           groupSelectsFiltered=True,
                           rowMultiSelectWithClick=True)

    gb.configure_grid_options()

    kwargs = {
        'gridOptions': gb.build(),
        'height': 800 if len(df.index) > 30 else None,
        'width': '100%',
        'data_return_mode': DataReturnMode.FILTERED_AND_SORTED,
        'update_mode': GridUpdateMode.GRID_CHANGED,
        'fit_columns_on_grid_load': False,
        'columns_auto_size_mode': ColumnsAutoSizeMode.FIT_CONTENTS,
        'theme': 'streamlit',
    }

    return kwargs


def sidebar_options() -> dict:
    """
    This function creates the sidebar with the options and returns a dictionary of them.
    """

    options = dict()
    st.sidebar.title("Options")
    Path('../BlastDatabases').mkdir(parents=True, exist_ok=True)
    dbs = [path.name for path in Path('./BlastDatabases').iterdir() if path.is_dir()]
    if dbs:
        db = st.sidebar.selectbox('Select Blast Database', dbs)
        options['db'] = Path('./BlastDatabases', db, 'blastdb')
    else:
        st.sidebar.info('No databases found. Please make one in the Manage Genome Databases section.')
    options['blast_mode'] = st.sidebar.radio('Select blast operation mode', ['tblastn', 'blastn'])

    options['max_target_seq'] = st.sidebar.number_input('Max sequences per query: ', min_value=1, value=500, step=100)
    options['threads'] = st.sidebar.number_input('Threads to use: ',
                                                 min_value=1, max_value=cpu_count(),
                                                 value=round(cpu_count() / 2), step=1)

    st.sidebar.subheader('Filter by:')
    options['perc_identity'] = st.sidebar.slider('identity percentage', 0.0, 100.0, 60.0, step=5.0)
    options['perc_alignment'] = st.sidebar.slider('query coverage percentage', 0.0, 100.0, 50.0, step=5.0)

    st.sidebar.subheader('Utilities')

    if st.sidebar.button("Clear Cache"):
        # Clear values from *all* all in-memory and on-disk data caches:
        st.cache_data.clear()
        st.sidebar.info('Cache cleared')

    if 'clear_analysis_btn_disabled' not in st.session_state:
        st.session_state['clear_analysis_btn_disabled'] = True

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

    st.session_state['options'] = options

    return options


def main():
    HEADER_TBLASTN = ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                      'identity', 'positive', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                      's_end', 'evalue', 'bit_score', 'sseq']

    HEADER_BLASTN = ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                     'identity', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                     's_end', 'evalue', 'bit_score', 'sseq']

    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='None')

    st.title('Blast queries against your local database!')

    options = sidebar_options()
    query = st.text_area('Enter your query', placeholder="Query...")

    if st.button("Blast query"):
        if 'db' not in options:
            st.warning('No databases found. Please make one in the Manage Genome Databases section.')
            st.stop()

        if not query:
            st.warning('Please enter a query')
            st.stop()

        with st.spinner('Blasting...'):
            st.session_state['BlastResponse']: BlastResponse = blast(query,
                                                                     blast_mode=options['blast_mode'],
                                                                     db=options['db'],
                                                                     threads=options['threads'],
                                                                     max_target_seq=options['max_target_seq'], )

        with st.spinner('Converting to dataframe...'):
            if options['blast_mode'] == 'tblastn':
                columns = HEADER_TBLASTN + ['id']
            elif options['blast_mode'] == 'blastn':
                columns = HEADER_BLASTN + ['id']
            st.session_state['df'] = st.session_state['BlastResponse'].to_dataframe(columns=columns)

    if 'df' in st.session_state:
        blast_response = st.session_state['BlastResponse']
        df = st.session_state['df']
        df = filter_results(df, options['perc_identity'], options['perc_alignment'])

        tab_table, tab_alignments = st.tabs(['Table', 'Alignments'])

        with tab_table:
            st.header("Blast analysis results")
            st.write("Click on a row to show the alignment, click again to hide it.")
            st.write("By clicking on a column you can sort it. When you hover on it, an icon will appear "
                     "that will enable you to filter or aggregate the data. For example you can group by Strain "
                     "to see how many matches each strain has, and click on the *\"count(Query)\"* column to sort it."
                     " It's very useful if you want to find gene duplicates.")



            aggrid_options = get_aggrid_options(df)
            grid = AgGrid(df, **aggrid_options)

            grid_df = grid['data']
            selected = grid['selected_rows']

            st.session_state['grid_df'] = grid_df
            st.session_state['selected'] = selected

            st.write(f"Found {grid_df.shape[0]} results")

            # Show alignments of selected rows
            if selected:
                row_alignments = list()
                indexes = extract_indexes(selected)
                alignments = blast_response.alignments(indexes=indexes)

                st.session_state['row_alignments'] = '\n'.join(row_alignments)

                for i, alignment in enumerate(alignments):
                    st.markdown(f"###### **{i})**")
                    st.code(alignment)

            # Initialize the download buttons
            download_container = st.container()
            set_download_buttons(container=download_container)

        with tab_alignments:
            st.header("Blast analysis alignments")

            if df.sseq.str.len().max() >= 1000:
                st.info("The alignments are too long to be displayed. If you want to check some, select the row from "
                        "the table and it will be shown underneath. If you want to see all the alignments "
                        "please download them by clicking the button in the Table tab. "
                        "Beware that it may take a few minutes.")
                st.stop()

            items_per_page = 50
            n_rows = grid_df.shape[0]
            n_pages = ceil(n_rows / items_per_page)  # Round UP
            selectbox_options = [f'{j * items_per_page + 1} - {min((j + 1) * items_per_page, n_rows)}' for j
                                 in range(n_pages)]
            page = st.selectbox('Page', selectbox_options, index=0)

            if not page:
                st.warning('No alignments to show')
                st.stop()

            start, end = page.split(' - ')
            start, end = int(start) - 1, int(end)
            indexes = grid_df['id'][start:end]
            alignments = blast_response.alignments(indexes=indexes)

            download_button = st.container()

            for i, alignment in enumerate(alignments):
                num_col, alignments_col = st.columns([1, 20])
                num_col.subheader(f"{i})")
                alignments_col.code(alignment)

            download_button._download_button(
                label="Download alignments in this page",
                data='\n'.join(alignments),
                file_name=f'alignments_{page}.txt',
                mime='text',
                key='alignments_download_button2')

            st.session_state['alignments'] = '\n'.join(alignments)


if __name__ == "__main__":
    main()
