import base64
import json
import os
import sys
from io import BytesIO
from math import ceil
from pathlib import Path, PurePath
from typing import Union

import pandas as pd
import streamlit as st
from st_aggrid import GridOptionsBuilder, AgGrid, DataReturnMode, ColumnsAutoSizeMode
from streamlit_extras.switch_page_button import switch_page
import streamlit.components.v1 as components

# Needed to search for scripts in the parent folder when using PyInstaller
sys.path.append(str(Path(__file__).parent))
from scripts import utils
from scripts.blast_parser import load_analysis, BlastParser, EmptyCSVError


def extract_indexes(selected: list):
    indexes = []
    for row in selected:
        nodeRowIndex = row['_selectedRowNodeInfo']['nodeRowIndex']
        match_index = row['id']
        indexes.append((nodeRowIndex, match_index))
    try:
        indexes = sorted(indexes, key=lambda x: x[0])
    except TypeError:
        # The grid was grouped by a column, so the nodeRowIndex is None, and you can't sort with None
        pass
    return zip(*indexes)


def download_table_xlsx():
    whole_df = st.session_state.blast_parser.whole_df
    grid_df: pd.DataFrame = st.session_state['grid_df']

    if grid_df.empty:
        st.warning('No rows to show')

    # Add sseq column to grid_df from blast_parser.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    data = utils.generate_xlsx_table(df)
    filename = 'blast.xlsx'

    # download_component_container points to an empty container at the end of the page that is used to
    # contain the component to download the file. Otherwise, the component would be rendered at the top shifting
    # all the page down (EVEN with height=0) and streamlit would then reload the page to reset it.
    with st.session_state['download_component_container']:
        components.html(html_download(data, filename), height=0)


def download_table_csv():
    whole_df = st.session_state.blast_parser.whole_df
    grid_df: pd.DataFrame = st.session_state['grid_df']

    if grid_df.empty:
        st.warning('No rows to show')

    # Add sseq column to grid_df from blast_parser.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    data: bytes = df.to_csv(index=False).encode('utf-8')
    filename = 'blast.tsv'

    # download_component_container points to an empty container at the end of the page that is used to
    # contain the component to download the file. Otherwise, the component would be rendered at the top shifting
    # all the page down (EVEN with height=0) and streamlit would then reload the page to reset it.
    with st.session_state['download_component_container']:
        components.html(html_download(data, filename), height=0)


def download_hit_sequences():
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    whole_df: pd.DataFrame = st.session_state.blast_parser.whole_df
    grid_df: pd.DataFrame = st.session_state.grid_df

    if grid_df.empty:
        st.warning('No rows to show')

    df = pd.merge(grid_df, whole_df[['id', 'sseq']], on=['id'], how='inner')
    df.insert(0, 'headers', df[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1))
    headers: list[str] = df['headers'].to_list()
    sequences: list[str] = df['sseq'].to_list()

    lines = list()
    for header, sequence in zip(headers, sequences):
        lines.append(header)
        # Split the sequence in lines of 60 characters
        lines.append('\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)]))

    data = '\n'.join(lines).encode('utf-8')
    filename = 'hits_sequences.fasta'

    # download_component_container points to an empty container at the end of the page that is used to
    # contain the component to download the file. Otherwise, the component would be rendered at the top shifting
    # all the page down (EVEN with height=0) and streamlit would then reload the page to reset it.
    with st.session_state['download_component_container']:
        components.html(html_download(data, filename), height=0)


def download_all_alignments():
    grid_df = st.session_state.grid_df
    blast_parser = st.session_state.blast_parser

    if grid_df.empty:
        st.warning('No alignments to download')

    alignments = blast_parser.alignments(indexes=grid_df['id'])
    data = '\n\n\n\n'.join(alignments).encode('utf-8')
    filename = 'alignments.txt'

    # download_component_container points to an empty container at the end of the page that is used to
    # contain the component to download the file. Otherwise, the component would be rendered at the top shifting
    # all the page down (EVEN with height=0) and streamlit would then reload the page to reset it.
    with st.session_state['download_component_container']:
        components.html(html_download(data, filename), height=0)


def html_download(object_to_download: Union[str, bytes], download_filename):
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

    # Suppress the virtualization of the columns. In this way the auto size of the column works as stated in the docs
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

    # Define custom CSS styling
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
        'row_buffer': 100,
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
        st.error('The selected file empty! '
                 "\\\nIt's likely that the BLAST search did not finish correctly.")
        st.stop()

    return blast_parser


def show_table():
    blast_parser = st.session_state['blast_parser']
    df = st.session_state['df']

    st.subheader('Download:')
    download_container = st.empty()
    set_download_buttons(download_container)

    if blast_parser.df.empty:
        st.warning('No results found by BLAST!')
        return

    if df.shape[0] > 50_000:
        st.warning(f'The table has {df.shape[0]} rows. The program may take a few seconds to '
                   f'load or may even crash.')

    aggrid_options: dict = load_aggrid_options(df)
    grid = AgGrid(df, **aggrid_options)

    grid_df: pd.DataFrame = grid['data']
    selected = grid['selected_rows']

    st.session_state['grid_df'] = grid_df
    st.session_state['selected'] = selected

    st.write(f"Found {grid_df.shape[0]} results")

    # Show alignments of selected rows
    if selected:
        # Extract indexes in the order shown in the grid
        row_indexes, indexes = extract_indexes(selected)
        row_indexes = row_indexes[:50]
        indexes = pd.Series(indexes[:50])

        alignments = blast_parser.alignments(indexes=indexes)

        if row_indexes:
            st.subheader(f"Showing alignments for the selected rows")

        whole_df = blast_parser.whole_df
        for row, selected_json_index, alignment in zip(row_indexes, indexes, alignments):
            num_col, alignments_col = st.columns([1, 10])

            try:
                num_col.subheader(f"{row + 1})")
            except TypeError:
                num_col.subheader(f"#)")

            if whole_df.iloc[[selected_json_index]].sseq.str.len().max() > 1000:
                alignments_col.info("The alignment is too long to be displayed. If you want to see even "
                                    "long alignments please download them all by clicking the button above "
                                    "the table.")
                continue

            alignments_col.code(alignment)


def show_alignments():
    blast_parser = st.session_state['blast_parser']
    if blast_parser.df.empty:
        st.warning('No results found by BLAST!')
        return

    grid_df = st.session_state['grid_df']
    if grid_df.empty:
        st.warning('No alignments to show')
        return

    st.write(f"Tip: The alignments are shown in the same order as the table. "
             f"The number on the left corresponds to the row in the table.")

    items_per_page = 50
    n_rows = grid_df.shape[0]
    n_pages = ceil(n_rows / items_per_page)  # Round UP
    selectbox_options = [f'{j * items_per_page + 1} - {min((j + 1) * items_per_page, n_rows)}' for j
                         in range(n_pages)]
    page = st.selectbox('Page', selectbox_options, index=0)

    start, end = page.split(' - ')
    start, end = int(start) - 1, int(end)
    indexes = grid_df['id'][start:end]

    alignments = blast_parser.alignments(indexes=indexes)

    whole_df = blast_parser.whole_df
    for i, index_alignment in enumerate(zip(indexes, alignments)):
        index, alignment = index_alignment
        num_col, alignments_col = st.columns([1, 10])

        num_col.subheader(f"{start + i + 1})")

        if whole_df.iloc[[index]].sseq.str.len().max() > 1000:
            alignments_col.info("The alignment is too long to be displayed. If you want to see even "
                                "long alignments please download them all by clicking the button above "
                                "the table.")
            continue

        alignments_col.code(alignment)


def show_graphic_summary():
    blast_parser = st.session_state['blast_parser']
    if blast_parser.df.empty:
        st.warning('No results found by BLAST!')
        return

    grid_df = st.session_state['grid_df']
    if grid_df.empty:
        st.warning('No alignments to show')
        return

    query_titles: pd.Series = blast_parser.whole_df['query_title'].unique()

    # Show query titles as options but return the query ids
    st.selectbox('Select the query:',
                 options=query_titles, key='query_selectbox')

    indexes = grid_df[grid_df['query_title'] == st.session_state['query_selectbox']]['id']

    if not list(indexes):
        st.info('No alignments to show')
        st.stop()

    max_hits = st.number_input('Max hits to show', min_value=1, max_value=3000, value=200, step=50)
    st.subheader(f'Distribution of the top {min(max_hits, len(indexes))} Hits')
    st.write('Sorted by Evalue and alignment percentage')

    st.bokeh_chart(blast_parser.plot_alignments_bokeh(indexes=indexes, max_hits=max_hits),
                   use_container_width=True)


def show_about():
    blast_parser = st.session_state['blast_parser']

    program = blast_parser.metadata['program']
    version = blast_parser.metadata['version']
    params = blast_parser.metadata['params']

    st.markdown(f"""
        ## Analysis made with {program.upper()} {version}
        """)

    if params:
        st.subheader('Parameters')
        params_text = ''

        for key, value in params.items():
            params_text += f'{key}: {value}\n'

        st.code(params_text)
    st.markdown(f"""---""")

    # try to look for the file containing the query
    query_file = str(blast_parser.file.with_suffix('.fasta')).replace('_results', '_query')
    if Path(query_file).exists():

        query_text = Path(query_file).read_text()
        query_seqs = query_text.split('>')[1:]

        st.download_button(
            label="Download query sequences as FASTA",
            data=query_text,
            file_name='query.fasta',
            mime='text/fasta')

        for index, query in enumerate(blast_parser.queries):
            header, seq = query_seqs[index].split('\n', maxsplit=1)
            seq = seq.replace('\n', '').strip()
            query_len = len(seq)

            st.markdown(f"""
            #### Query {index + 1}:
            **Query title**: {query['query_title']}\n
            **Query length**: {query_len}\n""")

            if query['hits'] == 0:
                st.warning('No hits found for this query')
            else:
                st.markdown(f"""
                **Hits**: {query['hits']}\n
                """)

            if query_len < 1000:
                # Split into lines of 60 characters
                seq = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])
                st.code(f">{header}\n{seq}")
            else:
                st.markdown('*The sequence is too long to be shown here. You can download it instead*')

    else:
        st.info(f'The file containing the query sequences *"{query_file}"* was not found.')


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

    st.session_state['blast_parser'] = blast_parser
    st.session_state['df'] = blast_parser.filtered_df(st.session_state['perc_identity'],
                                                      st.session_state['perc_alignment'])

    tabs = ['Table', 'Alignments', 'Graphic summary', 'About this analysis']
    col_table, col_alignments, col_graphic_summary, col_about = st.columns([1, 1, 1, 1])

    with col_table:
        if st.button(tabs[0], use_container_width=True):
            st.session_state['results_tab'] = tabs[0]

    with col_alignments:
        if st.button(tabs[1], use_container_width=True):
            st.session_state['results_tab'] = tabs[1]

    with col_graphic_summary:
        if st.button(tabs[2], use_container_width=True):
            st.session_state['results_tab'] = tabs[2]

    with col_about:
        if st.button(tabs[3], use_container_width=True):
            st.session_state['results_tab'] = tabs[3]

    if 'results_tab' not in st.session_state:
        st.session_state['results_tab'] = 'Table'

    if st.session_state['results_tab'] == 'Table':
        show_table()
    elif st.session_state['results_tab'] == 'Alignments':
        show_alignments()
    elif st.session_state['results_tab'] == 'Graphic summary':
        show_graphic_summary()
    elif st.session_state['results_tab'] == 'About this analysis':
        show_about()
    else:
        raise ValueError('Invalid tab selected')

    st.session_state['download_component_container'] = st.empty()


if __name__ == "__main__":
    main()
