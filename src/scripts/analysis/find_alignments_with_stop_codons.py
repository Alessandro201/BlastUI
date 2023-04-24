import random
from string import ascii_letters

import pandas as pd
import streamlit as st

from scripts.utils import generate_xlsx_table
from st_keyup import st_keyup

from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import numpy as np


def analyze(df, container):
    st.header('Hits with stop codons in the alignment')

    blast_parser = st.session_state.blast_parser
    whole_df = blast_parser.whole_df

    df_with_seqs: pd.DataFrame = __find_stop_codons_in_sequence(df, whole_df)

    if df_with_seqs.empty:
        container.info('No hits have stop codons inside the matched sequence.')
        return

    n_query = df_with_seqs['query_title'].unique().shape[0]
    n_matches = df_with_seqs.shape[0]
    container.write(f"Found {n_query} queries which matched {n_matches} sequences with stop codons inside.")
    container.write(f"Download all as:")

    __set_download_buttons(df_with_seqs, container)

    container.subheader('Search in the table:')
    # __set_download_buttons(df_with_seqs, container)
    query_col, btn_col = container.columns([7, 1])
    with query_col:
        filter_query = st_keyup('Search in the whole table', debounce=100)

    with btn_col:
        use_regex = st.checkbox('Use regex', key='use_regex', value=False)

    if filter_query:
        # query_title is the first column
        filtered_df = df_with_seqs[df_with_seqs['query_title'].str.contains(filter_query, case=False, regex=use_regex)]

        for column in df_with_seqs.columns[1:]:
            if is_string_dtype(df_with_seqs[column]):
                filtered_df_2 = df_with_seqs[
                    df_with_seqs[column].str.contains(filter_query, case=False, regex=use_regex)]

            elif is_numeric_dtype(df_with_seqs[column]):
                try:
                    float(filter_query)
                except ValueError:
                    continue

                # Search if the value is close to the filter_query, to avoid floating point errors
                filtered_df_2 = df_with_seqs[np.isclose(df_with_seqs[column].astype(float), float(filter_query))]
            else:
                continue

            filtered_df = pd.concat([filtered_df, filtered_df_2])

        if filtered_df.empty:
            container.info('The search did not match any results.')
            return

        __set_download_buttons(filtered_df, container)

        filtered_df = filtered_df.drop_duplicates()
        st.dataframe(filtered_df)
        st.write(f'Found {filtered_df.shape[0]} results.')
    else:
        st.dataframe(df_with_seqs)


@st.cache_data(show_spinner=False)
def __find_stop_codons_in_sequence(df, whole_df) -> pd.DataFrame:
    df_with_seqs = pd.merge(df, whole_df[['id', 'sseq']], on=['id'], how='inner')
    df_with_seqs = df_with_seqs[df_with_seqs['sseq'].str.count('\*') > 0]

    return df_with_seqs


@st.cache_data(show_spinner=False)
def __download_table_xlsx(df) -> bytes:
    # Add sseq column to grid_df from blast_parser.whole_df
    whole_df = st.session_state.blast_parser.whole_df
    df_with_seqs = pd.merge(df, whole_df[['id', 'sseq']], on=['id'], how='inner')

    return generate_xlsx_table(df_with_seqs)


@st.cache_data(show_spinner=False)
def __download_table_csv(df_with_seqs) -> bytes:
    table_data: bytes = df_with_seqs.to_csv(index=False).encode('utf-8')

    return table_data


@st.cache_data(show_spinner=False)
def __download_hit_sequences(df_with_seqs) -> bytes:
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    # insert without inplace
    headers = df_with_seqs[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1)
    headers: list[str] = headers
    sequences: list[str] = df_with_seqs['sseq'].to_list()

    lines = list()
    for header, sequence in zip(headers, sequences):
        lines.append(header)
        # Split the sequence in lines of 60 characters
        lines.append('\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)]))

    lines = '\n'.join(lines).encode('utf-8')
    return lines


@st.cache_data(show_spinner=False)
def __download_all_alignments(df_with_seqs) -> bytes:
    blast_parser = st.session_state.blast_parser

    alignments = blast_parser.alignments(indexes=df_with_seqs['id'])
    alignments = '\n\n\n\n'.join(alignments).encode('utf-8')
    return alignments


def __get_unique_keys(n=1) -> tuple:
    """Returns a tuple of n unique keys which are not already used by streamlit"""

    keys = set()
    while len(keys) < n:
        key = random.choices(ascii_letters, k=10)
        key = ''.join(key)
        if key not in st.session_state.keys():
            keys.add(key)

    return tuple(keys)


def __set_download_buttons(df, container=None):
    if not container:
        container = st

    col1, col2, col3, col4 = container.columns([1, 1, 1, 1])

    keys = __get_unique_keys(4)

    # Downloads the whole table
    with col1:
        st.download_button(
            label="Table as XLSX",
            file_name='multiple_hits.xlsx',
            data=__download_table_xlsx(df),
            mime='text/xlsx',
            use_container_width=True,
            key=keys[0])

    # Downloads only selected rows
    with col2:
        st.download_button(
            label='Table as CSV',
            file_name='multiple_hits.tsv',
            data=__download_table_csv(df),
            mime='text/tsv',
            use_container_width=True,
            key=keys[1])

    # Download all alignments
    with col3:
        st.download_button(
            label="FASTA (hit sequences)",
            file_name='multiple_hits_sequences.fasta',
            data=__download_hit_sequences(df),
            mime='text/fasta',
            use_container_width=True,
            key=keys[2])

    with col4:
        st.download_button(
            label="TEXT (all alignments)",
            file_name='multiple_hits_alignments.txt',
            data=__download_all_alignments(df),
            mime='text/txt',
            use_container_width=True,
            key=keys[3])
