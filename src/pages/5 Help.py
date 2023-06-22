import sys
from io import BytesIO
from pathlib import Path

import streamlit as st

# Needed to search for scripts in the parent folder when using PyInstaller
sys.path.append(str(Path(__file__).parent))
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

    st.header('Help page ❓')

    st.subheader('How to view results with the table')
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


if __name__ == "__main__":
    main()
