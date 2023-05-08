import shlex
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timedelta
from io import StringIO, BytesIO
from multiprocessing import cpu_count
from pathlib import Path

import psutil
import streamlit as st
from streamlit_extras.switch_page_button import switch_page
from streamlit_option_menu import option_menu

# Needed to search for scripts in the parent folder when using PyInstaller
sys.path.append(str(Path(__file__).parent))
from scripts.blast_parser import load_analysis, EmptyCSVError
from scripts import utils


def prepare_for_blast_command(query: str, blast_mode: str, db: str, threads: int = cpu_count() / 2, **kwargs):
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

    # Find a file name for query and result that doesn't exist
    today = datetime.today()
    query_file = f'./Analysis/{today:%Y%m%d_%H%M%S}_query.fasta'
    while Path(query_file).exists():
        today = (today + timedelta(seconds=1))
        query_file = f'./Analysis/{today:%Y%m%d_%H%M%S}_query.fasta'

    out_file = Path(f'./Analysis/{today:%Y%m%d_%H%M%S}_results.tsv')

    # Write query to file to be used by blast. If it doesn't have a header add it
    query = query.strip()
    if query[0] != '>':
        query = '>Query_1\n' + query

    # Write query to file
    Path(query_file).parent.mkdir(parents=True, exist_ok=True)
    Path(query_file).write_text(query)

    blast_exec = st.session_state['blast_exec']
    outfmt = "7 qaccver saccver nident pident qlen length qcovhsp gaps gapopen " \
             "mismatch positive ppos qstart qend sstart send qframe sframe score " \
             "evalue bitscore qseq sseq"

    cmd = shlex.split(f'"{blast_exec[blast_mode]}" -query "{query_file}" -db "{db}" -outfmt "{outfmt}" '
                      f'-out "{out_file}" -num_threads {threads} {additional_params}')

    return out_file, cmd


def kill_process_group(procid):
    """
    Terminate a process and all its children, if some fails it kills them.
    """

    try:
        parent = psutil.Process(procid)

        childrens = parent.children(recursive=True)
        for child in childrens:
            child.terminate()
        parent.terminate()

        gone, alive = psutil.wait_procs(childrens + [parent], timeout=3)
        for proc in alive:
            proc.kill()
    except psutil.NoSuchProcess as e:
        print(f'kill_process_group raised: NoSuchProcess{e}.')


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


def write_metadata(file: Path | str, metadata: dict):
    """
    This function writes the metadata to the output file of BLAST to preserve the params of the search.
    :param file: the file to write the metadata to
    :param metadata: the metadata to be written
    """

    file = Path(file)
    with open(file, 'a') as f:
        f.write("# [PARAMS]\n")
        f.write(f"# database:{st.session_state['db'].parent.name}\n")

        for key, value in metadata.items():
            f.write(f'# {key}:{value}\n')

        f.write("# [END PARAMS]\n")


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


def duplicated_headers(query: str) -> list | None:
    """
    This function checks if the query has duplicated headers.

    :param query: the query to be checked
    :return: The duplicated headers, None otherwise
    """

    headers = set()
    dup_headers = list()

    for line in query.splitlines():
        # Blank lines
        if len(line) == 0:
            continue

        line = line.strip()

        if line[0] == '>':
            # Remove the '>' from the header and remove leading and trailing whitespaces
            # Example: ">   header1" --> "header1"
            line = line[1:].strip()

            if line in headers:
                dup_headers.append(line)
            else:
                headers.add(line)

    return dup_headers if dup_headers else None


def duplicated_sequences(query: str) -> dict | None:
    """
    This function checks if the query has duplicated sequences.

    :param query: the query to be checked
    :return: The duplicated headers, None otherwise
    """

    seqs = defaultdict(list)

    queries = query.strip().split('>')

    for query in queries:
        # Blank lines
        if len(query) == 0:
            continue

        header, seq = query.split('\n', maxsplit=1)

        # Remove new lines to avoid missing identical sequences with different new lines in the middle
        seq = seq.replace('\n', '')

        seqs[seq].append(header)

    dup_seqs = {seq: headers for seq, headers in seqs.items() if len(headers) > 1}

    return dup_seqs if dup_seqs else None


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

            if options['matrix'] not in ('BLOSUM45', 'BLOSUM50', 'PAM30', 'PAM70', 'PAM250'):
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
        user_custom_commands = st.text_area('Insert additional commands or overwrite existing ones:', height=50,
                                            placeholder='Example: -strand both -ungapped ...')
        user_custom_commands = shlex.split(user_custom_commands)
        options['user_custom_commands'] = ''
        for index, item in enumerate(user_custom_commands):
            if item.startswith('-'):
                # if the next item is another option, then it is a flag
                if user_custom_commands[index + 1].startswith('-'):
                    options['user_custom_commands'] += item + ' '
                else:
                    options[item[1:]] = user_custom_commands[index + 1]

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
                       page_icon=BytesIO(utils.resource_path('./icon.png').read_bytes()))

    st.title('Blast queries against your local database!')
    sidebar_options()

    # Check that BLAST is installed
    if 'blast_exec' not in st.session_state:
        blast_exec = utils.get_programs_path()
        st.session_state['blast_exec'] = blast_exec

    if st.session_state['blast_exec'] is None:
        st.error('Could not find BLAST. Please download it in the home section.')
        if st.button('Go to home'):
            switch_page('Home')

        st.stop()

    ###### BLAST MODE ######
    blast_modes = ["BLASTN", "BLASTP", "BLASTX", 'TBLASTN', 'TBLASTX']
    icons = ['list-task', 'list-task', "list-task", 'list-task', 'list-task']
    DEFAULT_BLAST_MODE = 'BLASTN'
    default_index = blast_modes.index(DEFAULT_BLAST_MODE)

    st.session_state.blast_mode = option_menu('', options=blast_modes, icons=icons, menu_icon="gear",
                                              default_index=default_index, orientation="horizontal").lower()

    ###### QUERY ######
    query = st.text_area('Insert the queries: ', placeholder="Query...", height=200).strip()
    st.session_state.query = query

    uploaded_files = st.file_uploader("Alternatively upload queries in fasta format",
                                      type=["fasta", "faa"],
                                      accept_multiple_files=True)
    if uploaded_files:
        if st.session_state.query:
            st.warning('You have uploaded files and written a query. The query will be ignored.')
        st.session_state.query = read_query(uploaded_files)

    ###### BLAST OPTIONS ######
    db = choose_database()
    st.session_state.db = db

    if db is None:
        st.warning('No databases found. Please make one in the Manage Databases section.')
        st.stop()

    exp = st.expander('⚙️ Advanced options', expanded=False)
    set_advanced_options(exp)

    start_col, end_col, _ = st.columns([1, 1, 4])

    ###### RUN BLAST COMMAND IF PRESENT ######
    if 'command_to_run' in st.session_state:
        command = st.session_state['command_to_run']
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        st.session_state['process_pid'] = p.pid
        st.session_state['process'] = p
        st.session_state["blast_start_time"] = datetime.now()

        del st.session_state['command_to_run']

    ###### BLAST ######
    # Button disabled during blast process
    if start_col.button('Blast query',
                        disabled=bool(st.session_state.get('process_pid', False)),
                        use_container_width=True):
        st.session_state.switch_to_result_page = False

        if 'db' not in st.session_state:
            st.warning('No databases found. Please make one in the Manage Databases section.')
            st.stop()

        if len(st.session_state.get('query', '')) == 0:
            st.warning('Please enter a query!')
            st.stop()

        if headers := duplicated_headers(st.session_state['query']):
            if len(headers) > 1:
                headers = ' \n- '.join(headers)
            else:
                # Get a random element from the set, which is the only element
                headers = headers.pop()
            st.warning(f'The following headers are present more than one time: \n- {headers}')
            st.stop()

        if duplicates := duplicated_sequences(st.session_state['query']):
            for seq, headers in duplicates.items():
                headers = ' \n- '.join(headers)
                st.warning(f'The following sequences are identical: \n- {headers}')
            st.stop()

        st.markdown(f'Blast started at: {datetime.now():%d/%m/%Y %H:%M:%S}')

        blast_output_file, command = prepare_for_blast_command(query=st.session_state['query'],
                                                               blast_mode=st.session_state['blast_mode'],
                                                               db=st.session_state['db'],
                                                               threads=st.session_state['threads'],
                                                               **st.session_state['advanced_options'])

        # rerun to update button states and execute command
        st.session_state['command_to_run'] = command
        st.session_state['blast_output_file'] = blast_output_file
        st.experimental_rerun()

    # Button enabled during blast process
    if end_col.button('Stop process',
                      disabled=not bool(st.session_state.get('process_pid', False)),
                      use_container_width=True):
        process_pid = st.session_state.get('process_pid', None)
        proc = st.session_state.get('process', None)

        if not process_pid or not proc:
            st.experimental_rerun()

        if proc.poll():
            st.info('Process already stopped')

        else:
            proc = st.session_state['process']

            # A well-behaved application should finish communicating after it's killed to consume the output.
            print(f'INSIDE MAIN - Terminating process with pid: {process_pid}')

            kill_process_group(procid=process_pid)
            proc.communicate()

        del st.session_state['process_pid']
        del st.session_state['process']
        st.experimental_rerun()

    # If the process is running, show the spinner and wait for the process to finish
    process_pid = st.session_state.get('process_pid', None)
    if process_pid and psutil.pid_exists(process_pid):

        st.markdown(f"""
        Blast started at: {st.session_state["blast_start_time"]:%Y-%m-%d %H:%M:%S}\n
        """)

        with st.spinner(f"Running {st.session_state.blast_mode}..."):
            try:
                # communicate() will wait until the process terminates
                p = st.session_state['process']
                out, err = p.communicate()

                if p.returncode != 0:
                    raise subprocess.CalledProcessError(p.returncode, p.args, output=out, stderr=err)

                if err:
                    lines = err.splitlines()
                    if any('FASTA-Reader: Ignoring invalid residues at position(s):' in line for line in lines):
                        st.warning(f'The analysis finished but some residues are invalid. '
                                   f'Please check that you have selected the correct BLAST program: ')
                    else:
                        st.warning(f'The analysis finished but there were some errors: ')

                    st.write(f"Showing last {20 if len(lines) > 20 else len(lines)} lines of error output:")
                    st.code('\n'.join(lines[-20:]))

            except subprocess.CalledProcessError as e:
                stderr = e.stderr
                st.error(f'Error running blast: {stderr}')

                if 'BLAST Database error: No alias or index file found for nucleotide database' in stderr:
                    st.info(f'It seems you were trying to do a ***{st.session_state["blast_mode"].upper()}*** '
                            f'which requires a nucleotide database, but '
                            f'***{Path(st.session_state["db"]).parent.name}*** is a protein one.')
                elif 'BLAST Database error: No alias or index file found for protein database' in stderr:
                    st.info(f'It seems you were trying to do a ***{st.session_state["blast_mode"].upper()}*** '
                            f'which requires a protein database, but ***{Path(st.session_state["db"]).parent.name}*** '
                            f'is a nucleotide one.')
                elif "there's a line that doesn't look like plausible data, but it's not marked as defline" in stderr:
                    st.info(f"Error parsing blast results. It's likely that there is a wrong character "
                            f"in the query that BLAST does not know how to interpret. "
                            f"Please check the query and try again.")
                else:
                    raise e

                st.stop()

        del st.session_state['process_pid']
        del st.session_state['process']

        with st.spinner('Parsing results...'):

            blast_output_file = st.session_state['blast_output_file']
            write_metadata(blast_output_file, st.session_state['advanced_options'])

            try:
                st.session_state['blast_parser'] = load_analysis(blast_output_file)
            except EmptyCSVError:

                st.error(f"The analysis did not produce any result. No matches were found.")
                st.stop()

        st.session_state["blast_end_time"] = datetime.now()
        st.session_state['switch_to_result_page'] = True

        # By rerunning the app, we can re-enable the blast query button, disable the stop process button
        # and then switch to the results page
        st.experimental_rerun()

    ##### SWITCH PAGE #####
    if st.session_state.get('switch_to_result_page', False):

        if 'blast_start_time' in st.session_state and 'blast_end_time' in st.session_state:
            time_elapsed: timedelta = st.session_state["blast_end_time"] - st.session_state["blast_start_time"]
            st.markdown(f"""
            Blast started at: {st.session_state["blast_start_time"]:%Y-%m-%d %H:%M:%S}\n
            Blast finished at: {st.session_state["blast_end_time"]:%Y-%m-%d %H:%M:%S}\n
            Elapsed time: {utils.strfdelta(time_elapsed, "{H}h {M}m {S:02.0f}s")}
            """)

        blast_parser = st.session_state.blast_parser

        # Check whether any query did not produce any hits
        zero_hit_queries = list()
        for query in blast_parser.queries:
            query_title = query['query_title']
            hits = query['hits']
            if hits == 0:
                zero_hit_queries.append(query_title)

        st.markdown("""---""")
        if st.button('Go to results'):
            st.session_state.switch_to_result_page = False
            switch_page('Results')

        for query_title in zero_hit_queries:
            st.info(f'No hits found for query: \\\n*{query_title}*')


if __name__ == "__main__":
    main()
