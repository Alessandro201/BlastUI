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


def read_genomes(uploaded_files) -> list[GenomeData]:
    genomes = []
    for uploaded_file in uploaded_files:
        # To convert to a string based IO
        genome_io = StringIO(uploaded_file.getvalue().decode("utf-8"))
        genome_str = genome_io.read()
        genome = GenomeData(name=uploaded_file.name, genome=genome_str)
        genomes.append(genome)

    return genomes


def save_uploads_to_disk(genomes: list[GenomeData], location: Path):
    for genome in genomes:
        with open(location / genome.genome_name, 'w') as f:
            f.write(genome.genome_str)


def set_sidebar():
    options = dict()
    return options


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='None')

    options = set_sidebar()

    st.title("Visualize gff")

    st.write('Upload the genomes you want to use. To clear the list, just refresh the page.')

    uploaded_files = st.file_uploader("Upload genomes", type=["gff", "gbf"], accept_multiple_files=True)
    if uploaded_files:
        with st.spinner('Reading files...'):
            st.session_state['genomes'] = read_genomes(uploaded_files)
        st.write(f'You have uploaded {len(st.session_state["genomes"])} genomes.')
    else:
        del st.session_state['genomes']


if __name__ == '__main__':
    main()
