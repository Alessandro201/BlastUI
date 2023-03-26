import streamlit as st
import shutil
from string import whitespace
from io import StringIO
from pathlib import Path
from scripts.makeblastdb import *
from scripts.utils import *
from multiprocessing import cpu_count
from subprocess import CalledProcessError

from streamlit_extras.switch_page_button import switch_page
from streamlit_extras.stoggle import stoggle
from streamlit_extras.no_default_selectbox import selectbox as ndf_selectbox
from st_keyup import st_keyup


@st.cache_data(show_spinner=False)
def read_genomes(uploaded_files) -> list[GenomeData]:
    genomes = []
    for uploaded_file in uploaded_files:
        # To convert to a string based IO
        genome_io = StringIO(uploaded_file.getvalue().decode("utf-8"))
        genome_str = genome_io.read()
        genome = GenomeData(name=uploaded_file.name, genome=genome_str)
        genomes.append(genome)

    return genomes


def sidebar_options():
    st.sidebar.title('Options')
    st.sidebar.subheader('Blast database options')

    st.session_state['threads'] = st.sidebar.number_input('Threads to use: ',
                                                          min_value=1, max_value=cpu_count(),
                                                          value=round(cpu_count() / 2), step=1)


def check_db_name_validity(db_name: str):
    """
    Check if the name of the blast database is valid. The user could insert a path or a forbidden character
    which may lead to security issues. If the name is not valid it uses streamlit to inform the user.

    """
    blast_db_dir = Path.cwd() / 'BlastDatabases'
    test_path = Path(blast_db_dir / db_name).resolve()

    # symbols and whitespaces like tabs are not allowed, but spaces are
    unallowed_chars = ['\\', '/', ':', '*', '?', '"', '<', '>'] + [*whitespace.replace(' ', '')]
    if test_path.parent != Path(blast_db_dir).resolve() or \
            any([char in db_name for char in unallowed_chars]):
        st.error(f'Filename "{db_name!r}" is not valid. You cannot enter "\\/\:*?"<>|"')
        st.stop()

    if test_path.exists():
        st.warning(f'The database name "{db_name}" is already in use. Please choose another one.')


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=Path(resource_path('.'), 'icon.png').read_bytes())

    sidebar_options()
    st.title("Manage databases")

    # Check that BLAST is installed
    if 'blast_exec' not in st.session_state:
        print('Checking blast...')
        blast_exec = get_programs_path()
        st.session_state['blast_exec'] = blast_exec

    if st.session_state['blast_exec'] is None:
        st.error('Could not find BLAST. Please download it in the home section.')
        if st.button('Go to home'):
            switch_page('Home')
        st.stop()

    stoggle('❓ What is a blast database? ',
            """
    To blast a query against some genomes, you need to first create a blast database with them. 
    This process will take a few minutes, depending on the number of genomes and their size, 
    and it will only be done once. <br> <br>
    The database will be saved locally and no file will be shared online. As a matter of fact, 
    the application works entirely offline, apart from the first download of BLAST.
    """)

    st.write('')

    create_tab, manage_tab = st.tabs(['Create new database', 'Manage databases'])

    with create_tab:
        st.write('Upload the genomes you want to use. To clear the list refresh the page.')

        uploaded_files = st.file_uploader("Upload genomes", type=["fasta", "faa", 'fa'], accept_multiple_files=True)
        if uploaded_files:
            st.write(f'You have uploaded {len(uploaded_files)} genomes.')

        if uploaded_files:

            ##### OPTIONS #####
            with st.expander('DATABASE OPTIONS', expanded=True):

                new_db_name = st_keyup('Database name', value='my_database')

                if not new_db_name:
                    st.error('You must enter a database name')
                    st.stop()

                new_db_name = new_db_name.strip(whitespace)
                check_db_name_validity(new_db_name)
                st.session_state['new_db_name'] = new_db_name

                dbtype = st.radio('Database type:', ('Nucleotides', 'Proteins'))
                if dbtype == 'Nucleotides':
                    st.session_state['dbtype'] = 'nucl'
                else:
                    st.session_state['dbtype'] = 'prot'

                stoggle('❓ Why renaming headers?',
                        """
                        Blast requires that the headers of the fasta files are unique. If you have 
                        uploaded genomes with repeated headers, you can choose to rename them. 
                        Each contig will be renamed as follows: "[file_name]_NODE_[contig_number]" so 
                        be sure to not upload fasta files with the same name.
                        Ideally, the fasta files should have the same name as the genome they contain.
                        """)

                st.checkbox('Auto rename headers', key='rename_headers_checkbox', value=True)

                stoggle('❓ Why removing small contigs?',
                        """
                        You can choose to remove contigs under a specified length from your genomes before 
                        creating the blast database. This will remove many partial and broken contigs formed during 
                        the assembly which may introduce noise in the blast results. <br> <br>
                        If you have uploaded multifasta files of annotated proteins you should avoid 
                        removing small contigs, as you may remove actual proteins. 
                        Furthermore, you should have removed them before annotating the genome.
                        """)
                st.checkbox('Remove small contigs from the fasta files', key='remove_contigs_checkbox', value=True)

                st.number_input('Minimum contig length', key='min_length',
                                disabled=not st.session_state['remove_contigs_checkbox'],
                                value=1000, min_value=1, max_value=10000001, step=100)

            if st.button(label='Create database'):

                pbar = st.progress(0)
                with st.spinner('Reading files...'):
                    genomes = read_genomes(uploaded_files)
                blast_exec = st.session_state['blast_exec']
                makeblastdb = MakeBlastDB(genomes=genomes,
                                          db_name=st.session_state['new_db_name'],
                                          dbtype=st.session_state['dbtype'],
                                          threads=st.session_state['threads'],
                                          pbar=pbar,
                                          makeblastdb_exec=blast_exec['makeblastdb'],
                                          rename_headers=st.session_state['rename_headers_checkbox'])

                # Filter genomes
                if st.session_state['remove_contigs_checkbox']:
                    try:
                        makeblastdb.remove_small_contigs(st.session_state['min_length'])
                    except ValueError as e:
                        st.error(e)
                        st.stop()

                # Generating multifasta file to make blast db
                makeblastdb.generate_multifasta()

                # Creating blast database
                try:
                    with st.spinner('Creating blast database...'):
                        makeblastdb.run()
                except CalledProcessError as e:
                    st.error(e.stderr)
                    raise e

                st.success('Done!')

    with manage_tab:
        st.subheader('Here you can view the databases you have created:')

        Path(Path().cwd(), 'BlastDatabases').mkdir(exist_ok=True, parents=True)
        databases = list([path for path in Path(Path().cwd(), 'BlastDatabases').iterdir() if path.is_dir()])

        if databases:
            st.markdown(''.join(['🔹 ' + db.name + '<br>' for db in databases]), unsafe_allow_html=True)
        else:
            st.info('No databases found.')
            st.stop()

        ### RENAME DATABASE
        st.markdown("""
                    <br>

                    ##### Rename database
                    """, unsafe_allow_html=True)

        db = ndf_selectbox('Select database', [db.name for db in databases], key='choose_rename_db')

        if db:
            new_name = st_keyup('New name', value=f'{db}')
            check_db_name_validity(new_name)

        btn_disabled = False if db else True
        if st.button('Rename', disabled=btn_disabled):
            with st.empty():
                old_db_path = Path(Path().cwd(), 'BlastDatabases', db)
                new_db_path = Path(Path().cwd(), 'BlastDatabases', new_name.strip(whitespace))
                old_db_path.rename(new_db_path)

            st.session_state['database_renamed'] = True
            st.experimental_rerun()

        if 'database_renamed' in st.session_state:
            st.success('Database renamed!')
            del st.session_state['database_renamed']

        ### DELETE DATABASE
        st.markdown("""
                    <br>

                    ##### Delete database
                    """, unsafe_allow_html=True)

        db = ndf_selectbox('Select database', [db.name for db in databases], key='delete_rename_db')

        btn_disabled = False if db else True
        if st.button('Delete', disabled=btn_disabled):
            def delete_database():
                db_path = Path(Path().cwd(), 'BlastDatabases', db)
                shutil.rmtree(db_path)

                st.session_state['database_deleted'] = True

            st.warning(f"The database ***{db}*** will be deleted!")
            st.button('Confirm', on_click=delete_database)

        if 'database_deleted' in st.session_state:
            st.success('Database deleted!')
            del st.session_state['database_deleted']


if __name__ == "__main__":
    main()
