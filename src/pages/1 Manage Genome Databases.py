import streamlit as st
import shutil
from io import StringIO
from pathlib import Path
from scripts.makeblastdb import *
from scripts.utils import *
from multiprocessing import cpu_count
from subprocess import CalledProcessError


def save_uploads_to_disk(genomes: list[GenomeData], location: Path):
    for genome in genomes:
        with open(location / genome.genome_name, 'w') as f:
            f.write(genome.genome_str)


def read_genomes(uploaded_files) -> list[GenomeData]:
    genomes = []
    for uploaded_file in uploaded_files:
        # To convert to a string based IO
        genome_io = StringIO(uploaded_file.getvalue().decode("utf-8"))
        genome_str = genome_io.read()
        genome = GenomeData(name=uploaded_file.name, genome=genome_str)
        genomes.append(genome)

    return genomes


def set_sidebar():
    options = dict()

    st.sidebar.title('Options')

    st.sidebar.subheader('Blast database options')

    options['threads'] = st.sidebar.number_input('Threads to use: ',
                                                 min_value=1, max_value=cpu_count(),
                                                 value=round(cpu_count() / 2), step=1)

    return options


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='None')

    options = set_sidebar()

    st.title("Manage blast db")
    st.write('To blast a query against some genomes, you need to first create a blast database with them. '
             'This process will take a few minutes, depending on the number of genomes and their size, '
             'and it will only be done once but it will greatly speed up any blast you perform.')
    st.write('The database will be saved locally and no file will be shared online. As a matter of fact, '
             'the application works entirely offline, apart from the first download of blast executables.')

    create_tab, manage_tab = st.tabs(['Create blast database', 'Manage blast databases'])

    with create_tab:
        st.write('Upload the genomes you want to use. To clear the list, just refresh the page.')

        uploaded_files = st.file_uploader("Upload genomes", type=["fasta", "faa"], accept_multiple_files=True)
        if uploaded_files:
            with st.spinner('Reading files...'):
                st.session_state['genomes'] = read_genomes(uploaded_files)
            st.write(f'You have uploaded {len(st.session_state["genomes"])} genomes.')
        else:
            del st.session_state['genomes']

        expand = True if 'genomes' in st.session_state else False
        with st.expander('DATABASE OPTIONS', expanded=expand):

            if 'genomes' in st.session_state:
                st.text_input('Database name', key='new_db_name', value='my_database')

                dbtype = st.radio('Database type:', ('Nucleotides', 'Proteins'))
                if dbtype == 'Nucleotides':
                    st.session_state['dbtype'] = 'nucl'
                else:
                    st.session_state['dbtype'] = 'prot'

                st.write('You can choose to remove contigs under a specified length from your genomes '
                         'before creating the blast database. This will remove many partial contigs which may '
                         'introduce noise in the blast results. ')
                st.write('If you have uploaded multifasta files of annotated proteins you should avoid '
                         'removing small contigs, as you may remove actual proteins. '
                         'Furthermore, you should have removed them before annotating the genome.')

                st.checkbox('Remove small contigs', key='remove_contigs_checkbox', value=True)
                st.number_input('Minimum contig length', key='min_length',
                                disabled=not st.session_state['remove_contigs_checkbox'],
                                value=1000, min_value=1, max_value=10000001, step=100)

                st.write('Blast requires that the headers of the fasta files are unique. If you have '
                         'uploaded genomes with repeated headers, you can choose to rename them. '
                         'Each contig will be renamed as follows: "[genome_name]_NODE_[contig_number]" so '
                         'be sure to not upload genomes with the same name.')
                st.write('If you have uploaded multifasta files of annotated proteins you should avoid '
                         'renaming the headers as they contain the protein name or its function. '
                         'Blast requires the headers to be less than 50 character long, though.')
                st.checkbox('Rename headers', key='rename_headers_checkbox', value=True)

        if 'genomes' in st.session_state:
            if st.button(label='Create database'):

                pbar = st.progress(0)
                binaries_in = read_configs()['BLAST']['use_executables_in']
                makeblastdb_exec = get_program_path('makeblastdb', binaries_in=binaries_in)

                makeblastdb = MakeBlastDB(genomes=st.session_state['genomes'],
                                          db_name=st.session_state['new_db_name'],
                                          dbtype=st.session_state['dbtype'],
                                          threads=options['threads'],
                                          pbar=pbar,
                                          makeblastdb_exec=makeblastdb_exec,
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
        st.write('Here you can manage the blast databases you have created. '
                 'You can delete them or rename them.')

        st.subheader('Choose databases:')
        databases = list([path for path in Path(Path().cwd(), 'BlastDatabases').iterdir() if path.is_dir()])
        if databases:
            st.session_state['database'] = st.radio('Databases', [db.name for db in databases])
        else:
            st.info('No databases found.')

        st.markdown(f'##### Rename {st.session_state["database"]} database:')
        new_name = st.text_input('New name', value=f'{st.session_state["database"]}')
        if st.button('Rename'):
            with st.empty():
                old_db_path = Path(Path().cwd(), 'BlastDatabases', st.session_state['database'])
                new_db_path = Path(Path().cwd(), 'BlastDatabases', new_name)
                old_db_path.rename(new_db_path)

            st.success('Done!')

        st.markdown(f'##### Delete {st.session_state["database"]} database:')
        if st.button('Delete'):
            def delete_database():
                db_name = st.session_state['database']
                db_path = Path(Path().cwd(), 'BlastDatabases', db_name)
                shutil.rmtree(db_path)

                st.session_state['database_deleted'] = True

            st.warning(f"The database {st.session_state['database']} will be deleted!")
            st.button('Confirm', on_click=delete_database)

        if 'database_deleted' in st.session_state:
            st.success('Done!')
            del st.session_state['database_deleted']


if __name__ == "__main__":
    main()
