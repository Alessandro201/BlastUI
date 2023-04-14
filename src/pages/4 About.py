from io import BytesIO
from pathlib import Path

import streamlit as st

from scripts import utils


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       # st.set_page_config() requires a BytesIO object instead of raw bytes
                       page_icon=BytesIO(Path(utils.resource_path('./icon.png')).read_bytes()))

    st.title('About')
    st.markdown("""
    This is version 1.0 of BlastUI. The source code is available on [GitHub](https://github.com/Alessandro201/BlastUI).
    
    Blast+ is a software package that can be used to perform local and remote sequence similarity searches.
    If you want to cite Blast+ you can use the following:
    >Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: 
    architecture and applications. BMC Bioinformatics. 2009 Dec 15;10:421. 
    doi: 10.1186/1471-2105-10-421. PMID: 20003500; PMCID: PMC2803857.
    """)


if __name__ == "__main__":
    main()
