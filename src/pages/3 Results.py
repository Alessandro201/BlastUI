import sys
from pathlib import Path, PurePath

# Needed to search for scripts in the parent folder when using PyInstaller
sys.path.append(str(Path(__file__).parent))

import pandas as pd
import os
import base64
import json
from math import ceil
from typing import Union

import streamlit as st
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid, DataReturnMode, ColumnsAutoSizeMode
from streamlit_extras.no_default_selectbox import selectbox as ndf_selectbox
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_response import load_analysis, BlastResponse
from scripts import utils
from scripts.utils import fragile
from scripts.analysis import find_strain_with_multiple_hits, find_alignments_with_stop_codons

from io import BytesIO


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
    grid_df: pd.DataFrame = st.session_state['grid_df']
    if grid_df.empty:
        return 'No rows to show'.encode('utf-8')

    # Add hseq column to grid_df from blast_response.whole_df
    whole_df = st.session_state.blast_response.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df = df.drop(columns=['id', 'query_id'])

    filename = 'blast.xlsx'
    data = utils.generate_xlsx_table(df)
    components.html(html_download(data, filename), height=None, width=None)


def download_table_csv():
    grid_df: pd.DataFrame = st.session_state['grid_df']
    if grid_df.empty:
        return 'No rows to show'.encode('utf-8')

    # Add hseq column to grid_df from blast_response.whole_df
    whole_df = st.session_state.blast_response.whole_df
    df = pd.merge(grid_df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df = df.drop(columns=['id', 'query_id'])

    table_data: bytes = df.to_csv(index=False).encode('utf-8')

    filename = 'blast.tsv'
    components.html(html_download(table_data, filename), height=None, width=None)


def download_hit_sequences():
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    grid_df: pd.DataFrame = st.session_state.grid_df
    whole_df: pd.DataFrame = st.session_state.blast_response.whole_df

    if grid_df.empty:
        return 'No rows to show'

    df = pd.merge(grid_df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df = df.drop(columns=['id'])

    df.insert(0, 'headers', df[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1))
    headers: list[str] = df['headers'].to_list()
    sequences: list[str] = df['hseq'].to_list()

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
    blast_response = st.session_state.blast_response

    alignments = blast_response.alignments(indexes=grid_df['id'])
    alignments = '\n\n\n\n'.join(alignments).encode('utf-8')

    filename = 'alignments.txt'
    components.html(html_download(alignments, filename), height=None, width=None)


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
        <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
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
    gb.configure_column('query_id', hide=True)

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


def show_blast_errors():
    blast_response = st.session_state.blast_response

    for message in blast_response.messages:
        st.warning(message)


def choose_analysis_to_load() -> BlastResponse:
    # .iterdir() returns a generator, and sorted() has the side effect of consuming the generator.
    Path('./Analysis/').mkdir(parents=True, exist_ok=True)
    analysis_outputs = Path('./Analysis/').iterdir()
    analysis_outputs = sorted(analysis_outputs, reverse=True)

    if not analysis_outputs:
        st.warning('Please run a BLAST search first in "Blast Query" page.')

        if st.button('Go to "Blast Query" page'):
            switch_page('Blast Query')

        st.stop()

    # Remove queries from results
    analysis_outputs = [file for file in analysis_outputs if '_query.fasta' not in str(file)]
    analysis_outputs = [file.name for file in analysis_outputs]

    # add '### LAST ###' as prefix to the first element
    analysis_outputs[0] = '### LAST ### ' + analysis_outputs[0]

    json_to_load = st.selectbox('Select which analysis to load. Default: last', options=analysis_outputs, index=0,
                                format_func=lambda x: os.path.splitext(x)[0])

    # Remove the prefix
    if '### LAST ### ' in json_to_load:
        json_to_load = json_to_load.replace('### LAST ### ', '')

    json_to_load = Path('./Analysis/') / json_to_load

    try:
        with st.spinner('Loading analysis...'):
            blast_response = load_analysis(json_to_load)

    except json.JSONDecodeError:
        st.error('The selected file is not a valid JSON file! Probably the BLAST search did not finish correctly.')
        st.stop()

    return blast_response


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=BytesIO(utils.resource_path('./icon.png').read_bytes()))

    st.title('Blast results!')
    sidebar_options()

    # Check that BLAST is installed
    if 'blast_exec' not in st.session_state:
        blast_exec = utils.get_programs_path()
        st.session_state['blast_exec'] = blast_exec

    if st.session_state['blast_exec'] is None:
        st.error('Could not find BLAST. Please download it in the home section.')
        if st.button('Go to home'):
            switch_page('Home')
        st.stop()

    blast_response = choose_analysis_to_load()
    st.session_state.blast_response = blast_response

    ###### Blast results ######
    if blast_response.df.empty:
        st.warning('No results found!')
        st.stop()

    df = blast_response.filtered_df(st.session_state['perc_identity'], st.session_state['perc_alignment'])

    tabs = st.tabs(['Table', 'Alignments', 'Graphic summary', 'Analysis', 'About this analysis'])
    tab_table, tab_alignments, tab_graphic_summary, tab_analysis, tab_about = tabs

    with fragile(tab_table):

        st.subheader('Download')
        download_container = st.empty()
        set_download_buttons(download_container)

        tab_description = st.expander('❓ How to work with the table', expanded=False)
        with tab_description:
            st.markdown("""
            ##### 1) Main features
            You can work with this table in a similar fashion as you would do with a spreadsheet.
            Click on a column to sort it. You can also group by and filter the columns.
            
            """)

            # streamlit.image() cannot read from the temporary folder created by PyInstaller (sys._MEIPASS)
            # given by resource_path(), so we need to give it the bytes of the image.
            col1, col2, col3 = st.columns([1, 1, 1])
            col1.image(utils.resource_path('./media/hover_menu.png').read_bytes(),
                       "By hovering on a header you can access the menu of the column")
            col2.image(utils.resource_path('./media/hover_menu_group_by.png').read_bytes(),
                       "You can group the column by its values")
            col3.image(utils.resource_path('./media/hover_menu_filter_by.png').read_bytes(),
                       'You can filter the column and choose which values to show')

            st.image(utils.resource_path('./media/vertical_tabs_aggrid.png').read_bytes(),
                     "Alternatively, you can use the vertical tabs to access the same options")

            st.markdown("""
            ##### 1.1) Example
            For example, you can group by *"strain"* to see how many hits each strain has, and click on the
            *"count(query_title)"* column to sort it. It's very useful if you want to quickly find multiple
            matches for the same strain, like in a search for a duplicate gene. It works best if you have only 
            one query, otherwise it gets a bit messy.
            """)

            st.image(utils.resource_path('./media/group_by_strain.png').read_bytes(), "Group by strain")
            st.image(utils.resource_path('./media/group_by_strain2.png').read_bytes(),
                     "Sort by count(query_title). You can see two matches of \"Query_1\" in the same strain \"SEB053\"")

            st.markdown("""
            ##### 2) Selecting rows
            You can select rows by clicking on them. To select multiple rows, press *ctrl* while clicking, or
            click on a row and while pressing *shift* click another one to select all the rows in between.
            
            ##### 3) Resetting the table
            To reset any filtering, sorting, aggregation or selection you can click on the menu of any column and 
            click on "Reset Columns".
            """)

        if df.shape[0] > 100_000:
            st.warning(f'The table has {df.shape[0]} rows. The program may take a few seconds to '
                       f'load or may even crash.')

        aggrid_options: dict = load_aggrid_options(df)
        grid = AgGrid(df, **aggrid_options)

        grid_df = grid['data']
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

            alignments = blast_response.alignments(indexes=indexes)

            if row_indexes:
                st.subheader(f"Showing alignments for {len(row_indexes)} selected rows")

            whole_df = blast_response.whole_df
            for row, index, alignment in zip(row_indexes, indexes, alignments):
                num_col, alignments_col = st.columns([1, 10])

                try:
                    num_col.subheader(f"{row + 1})")
                except TypeError:
                    num_col.subheader(f"#)")

                if whole_df.iloc[[index]].hseq.str.len().max() > 1000:
                    alignments_col.info("The alignment is too long to be displayed. If you want to see even "
                                        "long alignments please download them all by clicking the button above "
                                        "the table.")
                    continue

                alignments_col.code(alignment)

    with fragile(tab_alignments):

        if grid_df.empty:
            st.warning('No alignments to show')
            raise fragile.Break()  # exit with statement

        items_per_page = 50
        n_rows = grid_df.shape[0]
        n_pages = ceil(n_rows / items_per_page)  # Round UP
        selectbox_options = [f'{j * items_per_page + 1} - {min((j + 1) * items_per_page, n_rows)}' for j
                             in range(n_pages)]
        page = st.selectbox('Page', selectbox_options, index=0)

        start, end = page.split(' - ')
        start, end = int(start) - 1, int(end)
        indexes = grid_df['id'][start:end]

        alignments = blast_response.alignments(indexes=indexes)

        whole_df = blast_response.whole_df
        for i, index_alignment in enumerate(zip(indexes, alignments)):
            index, alignment = index_alignment
            num_col, alignments_col = st.columns([1, 10])

            num_col.subheader(f"{start + i + 1})")

            if whole_df.iloc[[index]].hseq.str.len().max() > 1000:
                alignments_col.info("The alignment is too long to be displayed. If you want to see even "
                                    "long alignments please download them all by clicking the button above "
                                    "the table.")
                continue

            alignments_col.code(alignment)

    with fragile(tab_graphic_summary):

        if grid_df.empty:
            st.warning('No alignments to show')
            raise fragile.Break()

        # [{'query_id1': str, 'query_title': str, 'query_len': int}, ...]
        queries: list[dict] = blast_response.metadata['queries']
        queries: dict = {query['query_id']: query['query_title'] for query in queries}

        # Show query titles as options but return the query ids
        st.selectbox('Select the query:',
                     options=queries.keys(), key='query_selectbox', format_func=lambda x: queries[x])

        indexes = grid_df[grid_df['query_id'] == st.session_state['query_selectbox']]['id']

        if not list(indexes):
            st.info('No alignments to show')
            st.stop()

        max_hits = st.number_input('Max hits to show', min_value=1, max_value=3000, value=200, step=50)
        st.subheader(f'Distribution of the top {min(max_hits, len(indexes))} Hits')
        st.write('Sorted by Evalue and alignment percentage')

        st.bokeh_chart(blast_response.plot_alignments_bokeh(indexes=indexes, max_hits=max_hits),
                       use_container_width=True)

        # st.pyplot(blast_response.plot_alignments_bokeh(indexes=indexes, max_hits=max_hits))

    with fragile(tab_analysis):
        if grid_df.empty:
            st.warning('No hits to analyze')
            raise fragile.Break()

        analysis_choices = ['Find strain with multiple hits for the same query',
                            'Find strain with stop codons inside the alignment']

        choice = ndf_selectbox('Select the analysis you want to perform:', options=analysis_choices)

        if choice == 'Find strain with multiple hits for the same query':
            find_strain_with_multiple_hits.analyze(grid_df, container=st)
        elif choice == 'Find strain with stop codons inside the alignment':
            find_alignments_with_stop_codons.analyze(grid_df, container=st)

    with tab_about:

        program = blast_response.program
        program_and_version = blast_response.metadata['version']
        metadata = blast_response.metadata
        params = blast_response.metadata['params']

        if program == 'blastn':
            results_params = f"""
            **Evalue**: {params['expect']}\n
            **Gap penalty**: Existence: {params['gap_open']} 
            Extension: {params['gap_extend']}\n
            **Reward/Penalty**: Match: {params['sc_match']} 
            Mismatch: {params['sc_mismatch']}\n
            **Filter**: {params['filter']}\n
            """

        elif program == 'blastp':
            results_params = f"""
            **Matrix**: {params['matrix']}\n
            **Evalue**: {params['expect']}\n
            **Gap penalty**: Existence: {params['gap_open']} 
            Extension: {params['gap_extend']}\n
            **Filter**: {params['filter']}\n
            **Composition-based statistics**: {params['cbs']}\n
            """

        elif program == 'blastx':
            results_params = f"""
            **Matrix**: {params['matrix']}\n
            **Evalue**: {params['expect']}\n
            **Gap penalty**: Existence: {params['gap_open']} 
            Extension: {params['gap_extend']}\n
            **Filter**: {params['filter']}\n
            **Composition-based statistics**: {params['cbs']}\n
            **Query Gencode**: {params['query_gencode']}\n
            """

        elif program == 'tblastn':
            results_params = f"""
            **Matrix**: {params['matrix']}\n
            **Evalue**: {params['expect']}\n
            **Gap penalty**: Existence: {params['gap_open']} 
            Extension: {params['gap_extend']}\n
            **Filter**: {params['filter']}\n
            **Database Gencode**: {params['db_gencode']}\n
            """

        elif program == 'tblastx':
            results_params = f"""
            
            **Matrix**: {params['matrix']}\n
            **Evalue**: {params['expect']}\n
            **Filter**: {params['filter']}\n
            **Query Gencode**: {params['query_gencode']}\n
            **Database Gencode**: {params['db_gencode']}\n
            
            """

        else:
            results_params = ''

        st.markdown(f"""
            ## Analysis made with {program_and_version}
        
            #### Reference: 
            {metadata['reference']} 
        
            #### Options:
            {results_params}
            
            ---
            """)

        # try to look for the file containing the query
        try:
            query_file = str(blast_response.json_file).replace('_results.json', '_query.fasta')
            query_text = Path(query_file).read_text()
            query_seqs = query_text.split('>')[1:]

            st.download_button(
                label="Download query sequences",
                data=query_text,
                file_name='query.fasta',
                mime='text/fasta')

        except FileNotFoundError:
            st.info('The file containing the query sequences was not found')

        for index, query in enumerate(blast_response.metadata['queries']):
            st.markdown(f"""
            #### Query {index + 1}:
            **Query id**: {query['query_id']}\n
            **Query title**: {query['query_title']}\n
            **Query length**: {query['query_len']}""")

            # if the file containing the query doesn't exist, an error is raised
            try:
                if query['query_len'] < 1000:
                    header, seq = query_seqs[index].split('\n', maxsplit=1)

                    seq = seq.replace('\n', '')
                    seq = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])

                    st.markdown(f"""**Sequence**:""")
                    st.code(f">{header}\n{seq}")
                else:
                    st.markdown('*The sequence is too long to be shown here. You can download it instead*')
            except (IndexError, UnboundLocalError):
                pass


if __name__ == "__main__":
    main()
