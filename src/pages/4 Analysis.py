import base64
import json
import os
import sys
from io import BytesIO
from pathlib import Path, PurePath

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid, DataReturnMode, ColumnsAutoSizeMode
from streamlit_extras.no_default_selectbox import selectbox as ndf_selectbox
from streamlit_extras.switch_page_button import switch_page

# Needed to search for scripts in the parent folder when using PyInstaller
sys.path.append(str(Path(__file__).parent))

from scripts.blast_parser import load_analysis, BlastParser, EmptyCSVError
from scripts import utils
from scripts.analysis import find_strain_with_multiple_hits, find_alignments_with_stop_codons


def download_table_xlsx():
    grid_df: pd.DataFrame = st.session_state['grid_df']
    if grid_df.empty:
        return 'No rows to show'.encode('utf-8')

    # Add sseq column to grid_df from blast_parser.whole_df
    whole_df = st.session_state.blast_parser.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    filename = 'blast.xlsx'
    data = utils.generate_xlsx_table(df)
    components.html(html_download(data, filename), height=None, width=None)


def download_table_csv():
    grid_df: pd.DataFrame = st.session_state['grid_df']
    if grid_df.empty:
        return 'No rows to show'.encode('utf-8')

    # Add sseq column to grid_df from blast_parser.whole_df
    whole_df = st.session_state.blast_parser.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    table_data: bytes = df.to_csv(index=False).encode('utf-8')

    filename = 'blast.tsv'
    components.html(html_download(table_data, filename), height=None, width=None)


def download_hit_sequences():
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    grid_df: pd.DataFrame = st.session_state.grid_df
    whole_df: pd.DataFrame = st.session_state.blast_parser.whole_df

    if grid_df.empty:
        return 'No rows to show'

    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    df.insert(0, 'headers', df[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1))
    headers: list[str] = df['headers'].to_list()
    sequences: list[str] = df['sseq'].to_list()

    lines = list()
    for header, sequence in zip(headers, sequences):
        lines.append(header)
        # Split the sequence in lines of 60 characters
        lines.append('\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)]))

    lines = '\n'.join(lines).encode('utf-8')
    filename = 'hits_sequences.fasta'
    components.html(html_download(lines, filename), height=None, width=None)


def download_all_alignments():
    grid_df = st.session_state.grid_df
    blast_parser = st.session_state.blast_parser

    alignments = blast_parser.alignments(indexes=grid_df['id'])
    alignments = '\n\n\n\n'.join(alignments).encode('utf-8')

    filename = 'alignments.txt'
    components.html(html_download(alignments, filename), height=None, width=None)


def html_download(object_to_download: str | bytes, download_filename):
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

    match object_to_download:
        case bytes():
            bytes_object = object_to_download
        case str():
            bytes_object = object_to_download.encode()
        case _:
            # try to dump the object to json
            object_to_download = json.dumps(object_to_download)
            bytes_object = object_to_download.encode()

    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(bytes_object).decode()

    except AttributeError:
        b64 = base64.b64encode(object_to_download).decode()

    dl_link = f"""
        <script src="https://code.jquery.com/jquery-3.2.1.min.js"></script>
        <script>
        $('<a href="data:text/{ext};base64,{b64}" download="{download_filename}">')[0].click()
        </script>
        """
    return dl_link


def set_download_buttons(container=None):
    if not container:
        container = st

    col1, col2, col3, col4 = container.columns([1, 1, 1, 1])

    # Downloads the whole table
    with col1:
        st.button(
            label="Table as XLSX",
            on_click=download_table_xlsx,
            use_container_width=True)

    # Downloads only selected rows
    with col2:
        st.button(
            label='Table as CSV',
            on_click=download_table_csv,
            use_container_width=True)

    # Download all alignments
    with col3:
        st.button(
            label="FASTA (hit sequences)",
            on_click=download_hit_sequences,
            use_container_width=True)

    with col4:
        st.button(
            label="TEXT (all alignments)",
            on_click=download_all_alignments,
            use_container_width=True)


def load_aggrid_options(df: pd.DataFrame) -> AgGrid:
    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_side_bar()
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=100)
    gb.configure_selection('multiple',
                           use_checkbox=False,
                           header_checkbox=False,
                           groupSelectsChildren=True,
                           groupSelectsFiltered=True,
                           rowMultiSelectWithClick=True)

    gb.configure_default_column(groupable=True, value=True, enableRowGroup=True, aggFunc='avg')
    gb.configure_column("Row", valueGetter="node.rowIndex + 1", enableRowGroup=False, aggFunc='first', pinned='left')
    gb.configure_columns(['query_title', 'strain', 'node'], aggFunc='count')
    gb.configure_column('id', hide=True)

    # Suppress the virtualisation of the columns. In this way the autosize of the column works as stated in the docs
    # https://www.ag-grid.com/javascript-data-grid/column-sizing/#auto-size-columns
    other_options = {'suppressColumnVirtualisation': True}
    gb.configure_grid_options(**other_options)

    gridOptions = gb.build()

    # Move the row column to the first position
    row_columns_index = len(gridOptions['columnDefs']) - 1
    gridOptions['columnDefs'].insert(0, gridOptions['columnDefs'].pop(row_columns_index))

    if df.shape[0] > 25:
        height = {'height': 760}
    else:
        height = {'domLayout': 'autoHeight'}

    # Define custom CSS
    custom_css = {
        ".ag-row-hover": {
            "background-color": "#c84949 !important",
            "color": "white !important",
        },

        '.ag-theme-streamlit': {
            # '--ag-foreground-color': 'rgb(126, 46, 132)',
            '--ag-background-color': 'rgb(255, 255, 255)',
            '--ag-odd-row-background-color': '#f0f2f6',
            '--ag-header-foreground-color': '#2A3439',
            '--ag-header-background-color': '#D5E0F5',
            '--ag-header-column-resize-handle-color': 'rgb(126, 46, 132)',
            '--ag-selected-row-background-color': '#ffa0a0',
            '--ag-selected-row-color': '#ffffff',
            '--ag-font-size': '15px',
            '--ag-font-family': 'monospace',
        },

        '.ag-theme-streamlit-dark': {
            '--ag-foreground-color': '#CFCFCF',
            '--ag-background-color': '#1a1c21',
            '--ag-odd-row-background-color': '#292a34',
            '--ag-header-foreground-color': '#ebebeb',
            '--ag-header-background-color': '#515367',
            '--ag-header-column-resize-handle-color': '#fafafa',
            '--ag-selected-row-background-color': '#873838',
            '--ag-font-size': '15px',
            '--ag-font-family': 'monospace',
        }
    }

    kwargs = {
        'gridOptions': gridOptions,
        'width': '100%',
        **height,
        'data_return_mode': DataReturnMode.FILTERED_AND_SORTED,
        'update_on': ['modelUpdated'],
        'columns_auto_size_mode': ColumnsAutoSizeMode.FIT_CONTENTS,
        'fit_columns_on_grid_load': False,
        'theme': 'streamlit',
        'allow_unsafe_jscode': True,
        'custom_css': custom_css,
        'row_buffer': 20,
        'enable_quicksearch': False
    }

    return kwargs


def sidebar_options():
    """
    This function creates the sidebar with the options and returns a dictionary of them.
    """

    st.sidebar.title("Options")

    st.sidebar.subheader('Filter by:')
    st.session_state['perc_identity'] = st.sidebar.slider('identity percentage', 0.0, 100.0, 60.0, step=5.0)
    st.session_state['perc_alignment'] = st.sidebar.slider('query coverage percentage', 0.0, 100.0, 50.0, step=5.0)

    st.sidebar.subheader('Utilities')

    if st.sidebar.button("Clear Cache"):
        # Clear values from *all* all in-memory and on-disk data caches:
        st.cache_data.clear()
        st.sidebar.info('Cache cleared')

    if st.sidebar.button("Clear session state"):
        for key in st.session_state.keys():
            del st.session_state[key]

        st.session_state['perc_identity'] = 60.0
        st.session_state['perc_alignment'] = 50.0

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


def choose_analysis_to_load() -> BlastParser | None:
    Path('./Analysis/').mkdir(parents=True, exist_ok=True)
    analysis_outputs = Path('./Analysis/').glob('*.tsv')

    if not analysis_outputs:
        st.warning('Please run a BLAST search first in "Blast Query" page.')

        if st.button('Go to "Blast Query" page'):
            switch_page('Blast Query')

        st.stop()

    analysis_outputs = sorted(analysis_outputs, reverse=True)
    analysis_outputs = [file.name for file in analysis_outputs]

    # If there is a blast_parser in session_state, preselect the corresponding analysis
    if 'blast_parser' in st.session_state:
        blast_parser = st.session_state['blast_parser']
        file_loaded = blast_parser.file.name
        preselected_index = analysis_outputs.index(file_loaded)
    else:
        preselected_index = 0

    file_to_load = st.selectbox('Select which analysis to load. Default: last',
                                options=analysis_outputs,
                                index=preselected_index,
                                format_func=lambda x: os.path.splitext(x)[0])

    if not file_to_load:
        return None

    try:
        with st.spinner('Loading analysis...'):
            file_to_load = Path('./Analysis/') / file_to_load
            blast_parser = load_analysis(file_to_load)

    except EmptyCSVError:
        st.error('The selected file is empty! Probably the BLAST search did not finish correctly.')
        st.stop()

    return blast_parser


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=BytesIO(utils.resource_path('./icon.png').read_bytes()))

    st.title('Blast results!')
    sidebar_options()

    # Load analysis and show selectbox to choose it. By default, choose the analysis last loaded
    blast_parser = choose_analysis_to_load()
    if not blast_parser:
        st.warning('No analysis done yet. Run BLAST first!')
        st.stop()

    if blast_parser.df.empty:
        st.warning('No results found by BLAST!')
        st.stop()

    st.session_state['blast_parser'] = blast_parser
    df = blast_parser.filtered_df(st.session_state['perc_identity'], st.session_state['perc_alignment'])

    analysis_choices = ['Find strain with multiple hits for the same query',
                        'Find strain with stop codons inside the alignment']

    choice = ndf_selectbox('Select the analysis you want to perform:', options=analysis_choices)

    if choice == analysis_choices[0]:
        find_strain_with_multiple_hits.analyze(df, container=st)
    elif choice == analysis_choices[1]:
        find_alignments_with_stop_codons.analyze(df, container=st)


if __name__ == "__main__":
    main()
