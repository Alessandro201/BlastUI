import random
from string import ascii_letters

import pandas as pd
import streamlit as st

from scripts.utils import generate_xlsx_table


def analyze(df, container):
    st.header('Strain with multiple hits for the same query')

    blast_response = st.session_state.blast_response
    queries = blast_response.queries

    no_results = True
    for query in queries:
        query_title = query['query_title']
        query_id = query['query_id']

        temp_df = df[df['query_id'] == query_id]

        v = temp_df['strain'].value_counts()
        dup = temp_df[temp_df['strain'].isin(v.index[v.gt(1)])]

        if dup.empty:
            continue

        no_results = False

        container.subheader(query_title)

        dup.sort_values(by=['strain'], inplace=True)

        # The index is the row in grid_df shown in the table. The rows shown though add 1 to make it more readable.
        # We need to do the same here so the user can switch between the tables.
        dup.index = dup.index + 1
        __set_download_buttons(dup, container)

        container.dataframe(dup)

    if no_results:
        container.info('No queries have multiple hits for the same strain.')


def __download_table_xlsx(df) -> bytes:
    # Add hseq column to grid_df from blast_response.whole_df
    whole_df = st.session_state.blast_response.whole_df
    df_with_seqs = pd.merge(df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df_with_seqs = df_with_seqs.drop(columns=['id', 'query_id'])

    return generate_xlsx_table(df_with_seqs)


def __download_table_csv(df) -> bytes:
    # Add hseq column to grid_df from blast_response.whole_df
    whole_df = st.session_state.blast_response.whole_df
    df_with_seqs = pd.merge(df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df_with_seqs = df_with_seqs.drop(columns=['id', 'query_id'])

    table_data: bytes = df_with_seqs.to_csv(index=False).encode('utf-8')

    return table_data


def __download_hit_sequences(df) -> bytes:
    def get_header(strain, node, query_title):
        return f">{strain}_NODE_{node};{query_title}"

    whole_df: pd.DataFrame = st.session_state.blast_response.whole_df

    df_with_seqs = pd.merge(df, whole_df[['id', 'hseq']], on=['id'], how='inner')
    df_with_seqs = df_with_seqs.drop(columns=['id'])

    df_with_seqs.insert(0, 'headers',
                        df_with_seqs[['strain', 'node', 'query_title']].apply(lambda x: get_header(*x), axis=1))
    headers: list[str] = df_with_seqs['headers'].to_list()
    sequences: list[str] = df_with_seqs['hseq'].to_list()

    lines = list()
    for header, sequence in zip(headers, sequences):
        lines.append(header)
        # Split the sequence in lines of 60 characters
        lines.append('\n'.join([sequence[i:i + 60] for i in range(0, len(sequence), 60)]))

    lines = '\n'.join(lines).encode('utf-8')
    return lines


def __download_all_alignments(df) -> bytes:
    blast_response = st.session_state.blast_response

    alignments = blast_response.alignments(indexes=df['id'])
    alignments = '\n\n\n\n'.join(alignments).encode('utf-8')
    return alignments


def __get_unique_keys(n=1) -> tuple:
    """Returns a tuple of n unique keys which are not already used by streamlit"""

    keys = set()
    while len(keys) < n:
        key = random.choices(ascii_letters, k=10)
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
