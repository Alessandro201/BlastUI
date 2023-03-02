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
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_response import *
from scripts.utils import *

from timebudget import timebudget
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


def set_download_buttons(container=None):
    if not container:
        container = st

    col1, col2, col3, *_ = container.columns([1, 1, 1, 1])

    # Downloads the whole table
    with col1:
        st.download_button(
            label="Table as XLSX",
            data=download_table_xlsx(),
            file_name='blast.xlsx',
            mime='text/xlsx',
            use_container_width=True)

    # Downloads only selected rows
    with col2:
        st.download_button(
            label='Table as CSV',
            data=download_table_csv(),
            file_name='blast.tsv',
            mime='text/tsv',
            use_container_width=True)

    # Download all alignments
    with col3:
        st.download_button(
            label="FASTA (hits sequences)",
            data=download_aligned_sequences(),
            file_name='blast.tsv',
            mime='text/tsv',
            use_container_width=True)


def download_table_xlsx():
    output = BytesIO()

    # By setting the 'engine' in the ExcelWriter constructor.
    writer = pd.ExcelWriter(output, engine="xlsxwriter")

    grid_df: pd.DataFrame = st.session_state['grid_df'].drop(columns=['id', 'query_id'])

    # Convert the dataframe to an XlsxWriter Excel object.
    grid_df.to_excel(writer, sheet_name='Sheet1', index=False, engine='xlsxwriter')

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    # Get the dimensions of the dataframe.
    (max_row, max_col) = grid_df.shape

    # Create a list of column headers, to use in add_table().
    column_settings = []
    for header in grid_df.columns:
        column_settings.append({'header': header})

    # Add the table.
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings, 'style': 'Table Style Medium 11'})

    # Make the columns wider for clarity.
    worksheet.set_column(0, max_col - 1, 12)

    # Autofit columns
    worksheet.autofit()

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()
    output.seek(0)
    return output.getvalue()


def download_table_csv():
    grid_df: pd.DataFrame = st.session_state['grid_df']
    table_data: bytes = grid_df.to_csv(index=False).encode()
    return table_data


def download_aligned_sequences():
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    df: pd.DataFrame = st.session_state['grid_df'][['strain', 'node', 'query_title', 'sseq']]
    df.insert(0, 'headers', df[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1))
    headers: list[str] = df['headers'].to_list()
    sequences: list[str] = st.session_state['grid_df']['sseq'].to_list()

    lines = list()
    for header, sequence in zip(headers, sequences):
        lines.append(header)
        # Split the sequence in lines of 60 characters
        lines.append('\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)]))

    return '\n'.join(lines)


def load_aggrid_options(df: pd.DataFrame) -> AgGrid:
    gb = GridOptionsBuilder.from_dataframe(df)

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

    if df.sseq.str.len().max() >= 200:
        gb.configure_column('sseq', hide=True)

    gb.configure_side_bar()

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

    gb.configure_grid_options()

    built = gb.build()

    # Move the row column to the first position
    row_columns_index = len(built['columnDefs']) - 1
    built['columnDefs'].insert(0, built['columnDefs'].pop(row_columns_index))

    kwargs = {
        'gridOptions': built,
        'height': 800,
        'width': '100%',
        'data_return_mode': DataReturnMode.FILTERED_AND_SORTED,
        'update_on': ['modelUpdated'],
        'columns_auto_size_mode': ColumnsAutoSizeMode.FIT_CONTENTS,
        'fit_columns_on_grid_load': False,
        'theme': 'streamlit',
        'allow_unsafe_jscode': True,
        'custom_css': custom_css,
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


def show_alignments(alignments, start, container=None):
    if container is None:
        container = st

    for i, alignment in enumerate(alignments):
        num_col, alignments_col = container.columns([1, 10])
        num_col.subheader(f"{start + i + 1})")
        alignments_col.code(alignment)


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='🧬')

    sidebar_options()

    st.title('Blast results!')

    if 'df' not in st.session_state or 'BlastResponse' not in st.session_state:
        st.warning('Please run a BLAST search first in "Run Blast" page.')

        if st.button('Go to "Run Blast" page'):
            switch_page('Run Blast')

        if st.button('Reaload last Blast search'):
            result_files = list()
            for file in Path('./Analysis/').iterdir():
                if 'query.fasta' not in file.name:
                    result_files.append(file)

            last_analysis = sorted(result_files, reverse=True)
            if len(last_analysis) == 0:
                st.warning('No previous analysis found!')
                st.stop()

            last_analysis = last_analysis[0]
            st.session_state['BlastResponse'] = BlastResponse(json_file=last_analysis)
            st.session_state['df'] = st.session_state['BlastResponse'].df
            st.experimental_rerun()

        st.stop()

    ###### Blast results ######
    if 'df' in st.session_state:
        blast_response: BlastResponse = st.session_state['BlastResponse']
        df = blast_response.filtered_df(st.session_state['perc_identity'], st.session_state['perc_alignment'])

        tabs = st.tabs(['Table', 'Alignments', 'Graphic summary', 'About this analysis'])
        tab_table, tab_alignments, tab_graphic_summary, tab_about = tabs

        with tab_table:

            st.subheader('Download')
            dowmload_container = st.empty()

            # tab_description, *_ = st.tabs(['How to work with the table'])

            tab_description = st.expander('How to work with the table', expanded=False)
            with tab_description:
                st.markdown("""
                ##### 1) Main features
                You can work with this table in a similar fashion as you would do with a spreadsheet.
                You can sort columns, and filter/aggregate the data by hovering on them and 
                clicking on the menu icon that appears. Alternatively, you can perform the same steps with the 
                vertical tabs on the right.
                
                
                ##### 1.1) Example
                For example, you can group by *"Strain"* to see how many hits each strain has, and click on the
                *"count(query_title)"* column to sort it. It's very useful if you want to quickly find multiple
                matches for the same strain, like in a search for a duplicate gene.
                
                ##### 2) Selecting rows
                You can select rows by clicking on them. To select multiple rows, press *ctrl* while clicking, or
                click on a row and while pressing *shift* click another one to select all the rows in between.
                To reset any filtering, sorting, aggregation or selection you can click again on the *reset table*
                button just below the table on the right.
                
                """)

            aggrid_options: dict = load_aggrid_options(df)
            grid = AgGrid(df, **aggrid_options)

            grid_df = grid['data']
            selected = grid['selected_rows']

            st.session_state['grid_df'] = grid_df
            st.session_state['selected'] = selected

            st.write(f"Found {grid_df.shape[0]} results")

            set_download_buttons(dowmload_container)

            # Show alignments of selected rows
            if selected:
                row_alignments = list()
                row_indexes, match_ids = extract_indexes(selected)
                row_indexes = row_indexes[:50]
                match_ids = match_ids[:50]
                alignments = blast_response.alignments(indexes=match_ids)

                st.session_state['row_alignments'] = '\n'.join(row_alignments)

                if row_indexes:
                    st.subheader(f"Showing alignments for {len(row_indexes)} selected rows")

                for row, alignment in zip(row_indexes, alignments):
                    num_col, alignments_col = st.columns([1, 10])

                    try:
                        num_col.subheader(f"{row + 1})")
                    except TypeError:
                        num_col.subheader(f"#)")
                    alignments_col.code(alignment)

        with fragile(tab_alignments):

            if df.sseq.str.len().max() >= 1000:
                st.info("The alignments are too long to be displayed. If you want to check some, select the row from "
                        "the table and it will be shown underneath. If you want to see all the alignments "
                        "please download them by clicking the button in the Table tab. "
                        "It may take a few second to download them.")
                raise fragile.Break()  # exit with statement

            items_per_page = 50
            n_rows = grid_df.shape[0]
            n_pages = ceil(n_rows / items_per_page)  # Round UP
            selectbox_options = [f'{j * items_per_page + 1} - {min((j + 1) * items_per_page, n_rows)}' for j
                                 in range(n_pages)]
            page = st.selectbox('Page', selectbox_options, index=0)

            if not page:
                st.warning('No alignments to show')
                raise fragile.Break()  # exit with statement

            start, end = page.split(' - ')
            start, end = int(start) - 1, int(end)
            indexes = grid_df['id'][start:end]
            alignments = blast_response.alignments(indexes=indexes)

            download_btn = st.container()

            show_alignments(alignments, start)

            download_btn.download_button(
                label="Download alignments in this page",
                data='\n'.join(alignments),
                file_name=f'alignments_{page}.txt',
                mime='text',
                key='alignments_download_button2')

        with tab_graphic_summary:

            # {query_id1: query_title1, ...}
            queries: dict = blast_response.queries

            # Show query titles as options but return the query ids
            st.selectbox('Select the query:',
                         options=queries.keys(), key='query_selectbox', format_func=lambda x: queries[x])

            indexes = grid_df[grid_df['query_id'] == st.session_state['query_selectbox']]['id']

            if not list(indexes):
                st.info('No alignments to show')
                st.stop()

            max_hits = 200
            st.subheader(f'Distribution of the top {min(max_hits, len(indexes))} Hits')

            st.bokeh_chart(blast_response.plot_alignments_bokeh(indexes=indexes, max_hits=max_hits),
                           use_container_width=True)

        with tab_about:
            st.markdown(f"""
            Analysis made with {blast_response.json['BlastOutput2'][0]['report']['version']}
            
            #### Reference: 
            {blast_response.json['BlastOutput2'][0]['report']['reference']} 
            
            #### Options:
            **Database**: {blast_response.json['BlastOutput2'][0]['report']['search_target']['db']}
            
            **Matrix**: {blast_response.json['BlastOutput2'][0]['report']['params']['matrix']}\n
            **Gap penalty**: Existence: {blast_response.json['BlastOutput2'][0]['report']['params']['gap_open']} 
            Extension: {blast_response.json['BlastOutput2'][0]['report']['params']['gap_extend']}\n
            **Filter**: {blast_response.json['BlastOutput2'][0]['report']['params']['filter']}\n
            **CBS**: {blast_response.json['BlastOutput2'][0]['report']['params']['cbs']}\n
            **DB Gencode**: {blast_response.json['BlastOutput2'][0]['report']['params']['db_gencode']}\n
            """)

            query_text = Path('Analysis/query.fasta').read_text()
            query_seqs = query_text.split('>')[1:]

            st.download_button(
                label="Download query sequences",
                data=query_text,
                file_name='query.fasta',
                mime='text/fasta')

            for index, report in enumerate(blast_response.reports):
                st.markdown(f"""
                #### Query {index + 1}:
                **Query id**: {report['results']['search']['query_id']}\n
                **Query title**: {report['results']['search']['query_title']}\n
                **Query length**: {report['results']['search']['query_len']}""")

                if report['results']['search']['query_len'] < 1000:
                    header, seq = query_seqs[index].split('\n', maxsplit=1)

                    seq = seq.replace('\n', '')
                    seq = '\n'.join([seq[i:i + 60] for i in range(0, len(seq), 60)])

                    st.markdown(f"""**Sequence**:""")
                    st.code(f">{header}\n{seq}")
                else:
                    st.markdown('*The sequence is too long to be shown here. You can download it instead*')


if __name__ == "__main__":
    main()
