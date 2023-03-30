import sys
from pathlib import Path

# Needed to search for scripts in the parent folder
sys.path.append(str(Path(__file__).parent))

import shlex
import time
from io import StringIO, BytesIO
from multiprocessing import cpu_count
from datetime import datetime

from streamlit_option_menu import option_menu
from streamlit_extras.switch_page_button import switch_page

from scripts.blast_response import *
from scripts.utils import *


@st.cache_data(show_spinner=False)
def blast(query: str, blast_mode: str, db: str, threads: int = cpu_count() / 2, **kwargs) -> Path:
    """
    This function runs the blast command and returns the results as a dictionary.
    """

    additional_params = ''
    for key, value in kwargs.items():
        if value is None:
            continue

        if key == 'user_custom_commands':
            additional_params += f' {value}'
            continue

        additional_params += f' -{key} {value}'

    today = datetime.today()
    today_time = today.strftime("%Y%m%d_%H%M%S")

    # Write query to file to be used by blast. If it doesn't have a header add it
    query_file = f'./Analysis/{today_time}_query.fasta'
    query = query.strip()
    if query[0] != '>':
        query = '>Query_1\n' + query
    Path(query_file).parent.mkdir(parents=True, exist_ok=True)
    Path(query_file).write_text(query)

    out_file = Path(f'./Analysis/{today_time}_results.json')

    blast_exec = st.session_state['blast_exec']

    cmd = shlex.split(f'"{blast_exec[blast_mode]}" -query "{query_file}" -db "{db}" -outfmt 15 '
                      f'-out "{out_file}" -num_threads {threads} {additional_params}')

    run_command(cmd)
    return out_file


def choose_database(container=None):
    if not container:
        container = st

    Path('./BlastDatabases').mkdir(parents=True, exist_ok=True)
    dbs = [path.name for path in Path('./BlastDatabases').iterdir() if path.is_dir()]

    if dbs:

        previous_db = st.session_state.get('db', None)

        if previous_db:
            previous_db_index = dbs.index(previous_db.parent.name)
        else:
            previous_db_index = 0

        db = container.selectbox('Select Blast Database', dbs, index=previous_db_index)
        db = Path('./BlastDatabases', db, 'blastdb')
        return db

    return None


def read_query(uploaded_files: list[BytesIO]) -> str:
    """
    This function reads the query from the uploaded files and returns it as a string.

    :param uploaded_files: files uploaded with st.file_uploader to be read
    :return: the uploaded files read and joined as a single string
    """

    queries = list()
    for uploaded_file in uploaded_files:
        # To convert to a string based IO
        query_io = StringIO(uploaded_file.getvalue().decode("utf-8"))
        query = query_io.read().strip()
        queries.append(query)

    return '\n'.join(queries)


def set_advanced_options(container=None, blast_mode=None):
    if container is None:
        container = st

    if blast_mode is None:
        blast_mode = st.session_state.blast_mode

    options = dict()

    row1 = container.container()
    row1_col1, row1_col2 = row1.columns([1, 1])

    with row1_col1:
        options['evalue'] = st.select_slider('E-value: ', options=[10 ** i for i in range(-100, 4)], value=10)

    with row1_col2:
        options['max_target_seqs'] = st.number_input('Max sequences per query: ', min_value=1, value=500,
                                                     step=100)

    row2 = container.container()
    row2_col1, row2_col2 = row2.columns([1, 1])

    if blast_mode == 'blastn':
        with row2_col1:
            tasks = ['megablast', 'dc-megablast', 'blastn', 'blastn-short']
            options['task'] = st.selectbox('Task: ', options=tasks, index=tasks.index('megablast'))

            gp_megablast = ['Linear', 'Existence: 5 Extension: 2', 'Existence: 2 Extension: 2',
                            'Existence: 1 Extension: 2', 'Existence: 0 Extension: 2', 'Existence: 3 Extension: 1',
                            'Existence: 2 Extension: 1', 'Existence: 1 Extension: 1']

            gp_blastn = ['Existence: 4 Extension: 4', 'Existence: 2 Extension: 4',
                         'Existence: 0 Extension: 4', 'Existence: 3 Extension: 3', 'Existence: 6 Extension: 2',
                         'Existence: 5 Extension: 2', 'Existence: 4 Extension: 2', 'Existence: 2 Extension: 2']

            gap_penalty = {'megablast': gp_megablast,
                           'dc-megablast': gp_blastn,
                           'blastn': gp_blastn,
                           'blastn-short': gp_blastn}

            gp_defaults = {'megablast': 'Linear',
                           'dc-megablast': 'Existence: 5 Extension: 2',
                           'blastn': 'Existence: 5 Extension: 2',
                           'blastn-short': 'Existence: 5 Extension: 2'}

            penalty = st.selectbox('Gap penalty: ',
                                   options=gap_penalty[options['task']],
                                   index=gap_penalty[options['task']].index(
                                       gp_defaults[options['task']])
                                   )

            if penalty == 'Linear':
                gapopen, gapextend = 0, None
            else:
                # form: 'Existence: 5 Extension: 2' I only care about the numbers
                _, gapopen, _, gapextend = penalty.split()

            options['gapopen'] = gapopen
            options['gapextend'] = gapextend

        with row2_col2:
            word_sizes = {'megablast': [16, 20, 24, 28, 32, 48, 64, 128, 256],
                          'dc-megablast': [11, 12],
                          'blastn': [7, 11, 15],
                          'blastn-short': [7, 11, 15]}

            ws_defaults = {'megablast': 28, 'dc-megablast': 11, 'blastn': 11, 'blastn-short': 7}

            options['word_size'] = st.selectbox('Word size: ',
                                                options=word_sizes[options['task']],
                                                index=word_sizes[options['task']].index(ws_defaults[options['task']]))

            rp_options = ['Match: 1 Mismatch: -2', 'Match: 1 Mismatch: -3', 'Match: 1 Mismatch: -4',
                          'Match: 2 Mismatch: -3', 'Match: 4 Mismatch: -5', 'Match: 1 Mismatch: -1']

            rp_defaults = {'megablast': 'Match: 1 Mismatch: -2',
                           'dc-megablast': 'Match: 2 Mismatch: -3',
                           'blastn': 'Match: 2 Mismatch: -3',
                           'blastn-short': 'Match: 1 Mismatch: -3'}

            reward_penalty = st.selectbox('Reward/Penalty: ',
                                          options=rp_options,
                                          index=rp_options.index(rp_defaults[options['task']]))

            _, reward, _, penalty = reward_penalty.split()
            options['reward'] = reward
            options['penalty'] = penalty

    if blast_mode == 'blastp':
        with row2_col1:
            tasks = ['blastp', 'blastp-fast', 'blastp-short']
            options['task'] = st.selectbox('Task: ', options=tasks, index=tasks.index('blastp'))

            ws_defaults = {'blastp': 3, 'blastp-fast': 6, 'blastp-short': 2}
            ws_max_values = {'blastp': 7, 'blastp-fast': None, 'blastp-short': None}

            options['word_size'] = st.number_input('Word size: ', min_value=2, max_value=ws_max_values[options['task']],
                                                   value=ws_defaults[options['task']], step=1)

        with row2_col2:
            matrixes = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']

            matrix_defaults = {'blastp': 'BLOSUM62', 'blastp-fast': 'BLOSUM62', 'blastp-short': 'PAM30'}

            options['matrix'] = st.selectbox('Matrix: ', options=matrixes,
                                             index=matrixes.index(matrix_defaults[options['task']]))

            gap_penalty = ['Existence: 11 Extension: 2', 'Existence: 10 Extension: 2', 'Existence: 9 Extension: 2',
                           'Existence: 8 Extension: 2', 'Existence: 7 Extension: 2', 'Existence: 6 Extension: 2',
                           'Existence: 13 Extension: 1', 'Existence: 12 Extension: 1', 'Existence: 11 Extension: 1',
                           'Existence: 10 Extension: 1', 'Existence: 9 Extension: 1', ]

            gp_defaults = {'blastp': 'Existence: 11 Extension: 1',
                           'blastp-fast': 'Existence: 11 Extension: 1',
                           'blastp-short': 'Existence: 9 Extension: 1'}

            penalty = st.selectbox('Gap penalty: ',
                                   options=gap_penalty,
                                   index=gap_penalty.index(gp_defaults[options['task']]))

            if penalty == 'Linear':
                gapopen, gapextend = 0, None
            else:
                # form: 'Existence: 5 Extension: 2' I only care about the numbers
                _, gapopen, _, gapextend = penalty.split()

            options['gapopen'] = gapopen
            options['gapextend'] = gapextend

    if blast_mode == 'blastx':
        with row2_col1:
            tasks = ['blastx', 'blastx-fast']
            options['task'] = st.selectbox('Task: ', options=tasks, index=tasks.index('blastx'))

            ws_defaults = {'blastx': 3, 'blastx-fast': 6}
            ws_max_values = {'blastx': 7, 'blastx-fast': None}

            options['word_size'] = st.number_input('Word size: ', min_value=2, max_value=ws_max_values[options['task']],
                                                   value=ws_defaults[options['task']], step=1)

            genetic_code_options = {1: "Standard (1)",
                                    2: "Vertebrate Mitochondrial (2)",
                                    3: "Yeast Mitochondrial (3)",
                                    4: "Mold Mitochondrial; ... (4)",
                                    5: "Invertebrate Mitochondrial (5)",
                                    6: "Ciliate Nuclear; ... (6)",
                                    9: "Echinoderm Mitochondrial (9)",
                                    10: "Euplotid Nuclear (10)",
                                    11: "Bacteria and Archaea (11)",
                                    12: "Alternative Yeast Nuclear (12)",
                                    13: "Ascidian Mitochondrial (13)",
                                    14: "Flatworm Mitochondrial (14)",
                                    15: "Blepharisma Macronuclear (15)",
                                    16: "Chlorophycean Mitochondrial (16)",
                                    21: "Trematode Mitochondrial (21)",
                                    22: "Scenedesmus obliquus Mitochondrial (22)",
                                    23: "Thraustochytrium Mitochondrial (23)",
                                    24: "Pterobranchia Mitochondrial (24)",
                                    25: "Candidate Division SR1 and Gracilibacteria (25)",
                                    26: "Pachysolen tannophilus Nuclear (26)",
                                    27: "Karyorelict Nuclear (27)",
                                    28: "Condylostoma Nuclear (28)",
                                    29: "Mesodinium Nuclear (29)",
                                    30: "Peritrich Nuclear (30)",
                                    31: "Blastocrithidia Nuclear (31)"}

            options['query_gencode'] = st.selectbox('Query Genetic code: ',
                                                    options=genetic_code_options.keys(),
                                                    index=0,
                                                    format_func=lambda x: genetic_code_options[x])

        with row2_col2:
            matrixes = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']

            options['matrix'] = st.selectbox('Matrix: ', options=matrixes, index=matrixes.index('BLOSUM62'))

            gap_penalty = ['Existence: 11 Extension: 2', 'Existence: 10 Extension: 2', 'Existence: 9 Extension: 2',
                           'Existence: 8 Extension: 2', 'Existence: 7 Extension: 2', 'Existence: 6 Extension: 2',
                           'Existence: 13 Extension: 1', 'Existence: 12 Extension: 1', 'Existence: 11 Extension: 1',
                           'Existence: 10 Extension: 1', 'Existence: 9 Extension: 1', ]

            penalty = st.selectbox('Gap penalty: ', options=gap_penalty,
                                   index=gap_penalty.index('Existence: 11 Extension: 1'))

            if penalty == 'Linear':
                gapopen, gapextend = 0, None
            else:
                # form: 'Existence: 5 Extension: 2' I only care about the numbers
                _, gapopen, _, gapextend = penalty.split()

            options['gapopen'] = gapopen
            options['gapextend'] = gapextend

    if blast_mode == 'tblastn':
        with row2_col1:
            tasks = ['tblastn', 'tblastn-fast']
            options['task'] = st.selectbox('Task: ', options=tasks, index=tasks.index('tblastn'))

            ws_defaults = {'tblastn': 3, 'tblastn-fast': 6}
            ws_max_values = {'tblastn': 7, 'tblastn-fast': None}

            options['word_size'] = st.number_input('Word size: ', min_value=2, max_value=ws_max_values[options['task']],
                                                   value=ws_defaults[options['task']], step=1)

            genetic_code_options = {1: "Standard (1)",
                                    2: "Vertebrate Mitochondrial (2)",
                                    3: "Yeast Mitochondrial (3)",
                                    4: "Mold Mitochondrial; ... (4)",
                                    5: "Invertebrate Mitochondrial (5)",
                                    6: "Ciliate Nuclear; ... (6)",
                                    9: "Echinoderm Mitochondrial (9)",
                                    10: "Euplotid Nuclear (10)",
                                    11: "Bacteria and Archaea (11)",
                                    12: "Alternative Yeast Nuclear (12)",
                                    13: "Ascidian Mitochondrial (13)",
                                    14: "Flatworm Mitochondrial (14)",
                                    15: "Blepharisma Macronuclear (15)",
                                    16: "Chlorophycean Mitochondrial (16)",
                                    21: "Trematode Mitochondrial (21)",
                                    22: "Scenedesmus obliquus Mitochondrial (22)",
                                    23: "Thraustochytrium Mitochondrial (23)",
                                    24: "Pterobranchia Mitochondrial (24)",
                                    25: "Candidate Division SR1 and Gracilibacteria (25)",
                                    26: "Pachysolen tannophilus Nuclear (26)",
                                    27: "Karyorelict Nuclear (27)",
                                    28: "Condylostoma Nuclear (28)",
                                    29: "Mesodinium Nuclear (29)",
                                    30: "Peritrich Nuclear (30)",
                                    31: "Blastocrithidia Nuclear (31)"}

            options['db_gencode'] = st.selectbox('Database Genetic code: ',
                                                 options=genetic_code_options.keys(),
                                                 index=0,
                                                 format_func=lambda x: genetic_code_options[x])

        with row2_col2:
            matrixes = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']

            options['matrix'] = st.selectbox('Matrix: ', options=matrixes, index=matrixes.index('BLOSUM62'))

            if not options['matrix'] in ('BLOSUM45', 'BLOSUM50', 'PAM30', 'PAM70', 'PAM250'):
                gap_penalty = ['Existence: 11 Extension: 2', 'Existence: 10 Extension: 2', 'Existence: 9 Extension: 2',
                               'Existence: 8 Extension: 2', 'Existence: 7 Extension: 2', 'Existence: 6 Extension: 2',
                               'Existence: 13 Extension: 1', 'Existence: 12 Extension: 1', 'Existence: 11 Extension: 1',
                               'Existence: 10 Extension: 1', 'Existence: 9 Extension: 1', ]

                penalty = st.selectbox('Gap penalty: ', options=gap_penalty,
                                       index=gap_penalty.index('Existence: 11 Extension: 1'))

                if penalty == 'Linear':
                    gapopen, gapextend = 0, None
                else:
                    # form: 'Existence: 5 Extension: 2' I only care about the numbers
                    _, gapopen, _, gapextend = penalty.split()

                options['gapopen'] = gapopen
                options['gapextend'] = gapextend

    if blast_mode == 'tblastx':
        with row2_col1:
            options['word_size'] = st.number_input('Word size: ', min_value=2, value=3, step=1)

            genetic_code_options = {1: "Standard (1)",
                                    2: "Vertebrate Mitochondrial (2)",
                                    3: "Yeast Mitochondrial (3)",
                                    4: "Mold Mitochondrial; ... (4)",
                                    5: "Invertebrate Mitochondrial (5)",
                                    6: "Ciliate Nuclear; ... (6)",
                                    9: "Echinoderm Mitochondrial (9)",
                                    10: "Euplotid Nuclear (10)",
                                    11: "Bacteria and Archaea (11)",
                                    12: "Alternative Yeast Nuclear (12)",
                                    13: "Ascidian Mitochondrial (13)",
                                    14: "Flatworm Mitochondrial (14)",
                                    15: "Blepharisma Macronuclear (15)",
                                    16: "Chlorophycean Mitochondrial (16)",
                                    21: "Trematode Mitochondrial (21)",
                                    22: "Scenedesmus obliquus Mitochondrial (22)",
                                    23: "Thraustochytrium Mitochondrial (23)",
                                    24: "Pterobranchia Mitochondrial (24)",
                                    25: "Candidate Division SR1 and Gracilibacteria (25)",
                                    26: "Pachysolen tannophilus Nuclear (26)",
                                    27: "Karyorelict Nuclear (27)",
                                    28: "Condylostoma Nuclear (28)",
                                    29: "Mesodinium Nuclear (29)",
                                    30: "Peritrich Nuclear (30)",
                                    31: "Blastocrithidia Nuclear (31)"}

            options['query_gencode'] = st.selectbox('Query Genetic code: ',
                                                    options=genetic_code_options.keys(),
                                                    index=0,
                                                    format_func=lambda x: genetic_code_options[x])

        with row2_col2:
            matrixes = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']

            options['matrix'] = st.selectbox('Matrix: ', options=matrixes, index=matrixes.index('BLOSUM62'))

            options['db_gencode'] = st.selectbox('Database Genetic code: ',
                                                 options=genetic_code_options.keys(),
                                                 index=0,
                                                 format_func=lambda x: genetic_code_options[x])

    row3 = container.container()
    with row3:
        options['user_custom_commands'] = st.text_area('Insert additional commands:', height=50,
                                                       placeholder='Example: -strand both -ungapped ...')

    st.session_state.advanced_options = options


def sidebar_options():
    """
    This function creates the sidebar with the options and returns a dictionary of them.
    """

    st.sidebar.title("Options")

    st.session_state['threads'] = st.sidebar.number_input('Threads to use: ',
                                                          min_value=1, max_value=cpu_count(),
                                                          value=round(cpu_count() / 2), step=1)

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


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon=Path(resource_path('.'), 'icon.png').read_bytes())

    st.title('Blast queries against your local database!')
    sidebar_options()

    # Check that BLAST is installed
    if 'blast_exec' not in st.session_state:
        blast_exec = get_programs_path()
        st.session_state['blast_exec'] = blast_exec

    if st.session_state['blast_exec'] is None:
        st.error('Could not find BLAST. Please download it in the home section.')
        if st.button('Go to home'):
            switch_page('Home')
        st.stop()

    ###### BLAST MODE ######
    blast_modes = ["BLASTN", "BLASTP", "BLASTX", 'TBLASTN', 'TBLASTX']
    icons = ['list-task', 'list-task', "list-task", 'list-task', 'list-task']

    if 'blast_mode' in st.session_state:
        default = blast_modes.index(st.session_state.blast_mode.upper())
    else:
        default = 3
    st.session_state.blast_mode = option_menu('', options=blast_modes, icons=icons, menu_icon="gear",
                                              default_index=default, orientation="horizontal").lower()

    ###### QUERY ######
    query = st.text_area('Insert the queries: ', placeholder="Query...", height=200).strip()
    st.session_state.query = query

    uploaded_files = st.file_uploader("Alternatively upload queries in fasta format", type=["fasta", "faa"],
                                      accept_multiple_files=True)
    if uploaded_files:
        if st.session_state.query:
            st.warning('You have uploaded files and written a query. The query will be ignored.')
        st.session_state.query = read_query(uploaded_files)

    ###### BLAST OPTIONS ######
    db = choose_database()
    st.session_state.db = db

    if db is None:
        st.warning('No databases found. Please make one in the Manage Genome Databases section.')
        st.stop()

    exp = st.expander('Advanced options', expanded=False)
    set_advanced_options(exp)

    ###### BLAST ######
    if st.button('Blast query'):
        st.session_state.switch_to_result_page = False

        if 'db' not in st.session_state:
            st.warning('No databases found. Please make one in the Manage Genome Databases section.')
            st.stop()

        if len(st.session_state.get('query', '')) == 0:
            st.warning('Please enter a query!')
            st.stop()

        try:
            with st.spinner(f"Running {st.session_state.blast_mode}..."):
                st.markdown(f'Blast started at: {datetime.now()}')

                blast_output_file = blast(query=st.session_state['query'],
                                          blast_mode=st.session_state['blast_mode'],
                                          db=st.session_state['db'],
                                          threads=st.session_state['threads'],
                                          **st.session_state['advanced_options'])

                st.session_state['blast_output_file'] = blast_output_file

        except subprocess.CalledProcessError as e:
            stderr = e.stderr
            st.error(f'Error running blast: {stderr}')

            if 'BLAST Database error: No alias or index file found for nucleotide database' in stderr:
                st.info(f'It seems you were trying to do a ***{st.session_state["blast_mode"].upper()}*** '
                        f'which requires a nucleotide database, but ***{Path(st.session_state["db"]).parent.name}*** '
                        f'is a protein one.')
            elif 'BLAST Database error: No alias or index file found for protein database' in stderr:
                st.info(f'It seems you were trying to do a ***{st.session_state["blast_mode"].upper()}*** '
                        f'which requires a protein database, but ***{Path(st.session_state["db"]).parent.name}*** '
                        f'is a nucleotide one.')

            st.stop()

        with st.spinner('Parsing results...'):
            st.session_state['blast_response'] = load_analysis(blast_output_file)

        st.session_state['switch_to_result_page'] = True

    ##### SWITCH PAGE #####
    if st.session_state.get('switch_to_result_page', False):
        blast_response = st.session_state.blast_response

        if not blast_response.messages:
            st.session_state.switch_to_result_page = False
            switch_page('Results')

        for message in blast_response.messages:
            st.warning(message)

        if st.button('Go to results'):
            st.session_state.switch_to_result_page = False
            switch_page('Results')


if __name__ == "__main__":
    main()
