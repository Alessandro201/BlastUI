import streamlit as st
import shutil
from io import StringIO
from pathlib import Path
from scripts.blast_utilities import *


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
        genome = GenomeData(genome_name=uploaded_file.name, genome_str=genome_str)
        genomes.append(genome)

    return genomes


def clear_genomes_folder(location: Path):
    for file in location.iterdir():
        file.unlink()


def set_sidebar():
    options = dict()

    st.sidebar.title('Options')
    st.sidebar.write('You can choose to remove contigs under a specified length from your genomes '
                     'before creating the blast database. This will remove many partial contigs which only '
                     'introduce noise in the blast results. ')

    options['remove_small_contigs'] = False
    if st.sidebar.checkbox('Remove small contigs', key='remove_small_contigs', value=True):
        min_length = st.sidebar.number_input('Minimum contig length',
                                             value=1000, min_value=1, max_value=10000001, step=100)

        options['remove_small_contigs'] = True
        options['min_length'] = min_length

    st.sidebar.subheader('Blast database options')
    binaries = st.sidebar.radio('Use the blast binaries:', ('downloaded with BlastUI', 'in PATH'))

    if binaries == 'downloaded with BlastUI':
        options['binaries'] = 'BlastUI'
    else:
        options['binaries'] = 'PATH'

    dbtype = st.sidebar.radio('Database type:', ('Nucleotides', 'Proteins'))
    if dbtype == 'Nucleotides':
        options['dbtype'] = 'nucl'
    else:
        options['dbtype'] = 'prot'

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
             'this application works entirely offline.')

    create_tab, manage_tab = st.tabs(['Create blast database', 'Manage blast databases'])

    with create_tab:
        st.write('Upload the genomes you want to use. To clear the list, just refresh the page.')

        uploaded_files = st.file_uploader("Upload genomes", type=["fasta"], accept_multiple_files=True)
        if uploaded_files:
            with st.spinner('Reading files...'):
                st.session_state['genomes'] = read_genomes(uploaded_files)
            st.write(f'You have uploaded {len(st.session_state["genomes"])} genomes.')

        if 'genomes' in st.session_state:

            st.markdown('##### Insert the Database name')
            st.session_state['new_database'] = st.text_input('Database name', value='my_database')

            if st.button('Create blast database'):
                with st.empty():
                    # Filter genomes
                    genomes = st.session_state['genomes']
                    if st.session_state['remove_small_contigs']:
                        with st.spinner('Removing small contigs...'):
                            genomes: list[GenomeData] = remove_small_contigs(genomes, options['min_length'])

                    # Generating multifasta file to make blast db
                    with st.spinner('Generating multifasta file for makeblastdb...'):
                        multifasta_path = Path('tempData/multifasta.fasta')
                        make_multifasta_for_blastdb(genomes, multifasta_path)

                    # Creating blast database
                    with st.spinner('Creating blast database...'):
                        db_path = Path('BlastDatabases', st.session_state['new_database'], 'blastdb')
                        db_path.parent.mkdir(parents=True, exist_ok=True)
                        makeblastdb_exec = get_program_path('makeblastdb', binaries_in=options['binaries'])

                        make_blast_db(multifasta_path,
                                      database=db_path,
                                      title=st.session_state['new_database'],
                                      dbtype=options['dbtype'],
                                      makeblastdb_exec=makeblastdb_exec)

                        shutil.rmtree(multifasta_path.parent)


                st.success('Done!')

    with manage_tab:
        st.write('Here you can manage the blast databases you have created. '
                 'You can delete them or rename them.')

        st.subheader('Choose databases:')
        databases = list(Path(Path().cwd(), 'BlastDatabases').iterdir())
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

            st.success('Done!')

        if 'database_deleted' in st.session_state:
            st.success('Done!')
            del st.session_state['database_deleted']


if __name__ == "__main__":
    main()
