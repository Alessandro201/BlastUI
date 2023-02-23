#!/usr/bin/env python3

from os.path import splitext
from multiprocessing import Pool, cpu_count
from pathlib import Path
import streamlit as st
import sys
import shutil
import shlex
from subprocess import Popen, PIPE, CalledProcessError


class GenomeData:
    def __init__(self, genome_name: str, genome_str: str):
        self.genome_name: str = genome_name
        self.genome_str: str = genome_str


class _WorkerData:
    def __init__(self, genome_name: str, genome_str: str, min_length: int):
        self.genome_name = genome_name
        self.genome_str = genome_str
        self.min_length = min_length


##### Remove small contigs #####
def _worker_remove_small_contigs(data: _WorkerData) -> GenomeData:
    """
    Remove all small contigs from a file and rewrites it.

    :param data: Dictionary containing the genome as string and the threshold length of the contigs.
    :return:
    """

    genome_name = data.genome_name
    genome_str = data.genome_str
    min_length = data.min_length

    # Don't keep the first empty sequence
    sequences = genome_str.split('>')[1:]

    filtered_sequences = list()
    for index, sequence in enumerate(sequences):
        header, seq = sequence.split('\n', maxsplit=1)
        seq = seq.replace('\n', '')
        if len(seq) > min_length:
            filtered_sequences.append(sequence)

    genome = GenomeData(genome_name=genome_name,
                        genome_str='>' + '>'.join(filtered_sequences))
    return genome


def remove_small_contigs(genomes: list[GenomeData], min_length) -> list[GenomeData]:
    """
    :param genomes: tuple containing the genome name and the genome as string
    :param min_length: threshold length of the contigs
    :return:
    """
    queue = list()

    for genome in genomes:
        # This is just a way to pass all the needed data to remove_small_contigs function which will be
        # called by the Pool of processes, as sharing data between processes is tricky.

        data = _WorkerData(
            genome_name=genome.genome_name,
            genome_str=genome.genome_str,
            min_length=min_length
        )

        queue.append(data)

    with Pool(cpu_count()) as pool:
        filtered_genomes = pool.map(_worker_remove_small_contigs, queue)

    return filtered_genomes


##### Generate multifasta for blastdb #####
def _rename_fasta_headers(genome: GenomeData) -> str:
    """
    Read the fasta file, renaming each scaffold with the name of the file and a counter.
    It's both for clarity (1) and a bit more security (2).

    1) When you use blast it will output the name of the contig in which a match was found, but the assembler I
    used, SPAdes, does not use any meaningful name for the scaffolds. This is a way to know straight away the genome
    which gave a match.
    2) It happened that multiple genomes had the same prefix for the contigs. If you use genomes from the same folder
    and with the same extension you are sure that you will not get any duplicate name.

    :param genome: GenomeData object containing the genome name and the genome as string
    :return:
    """

    count = 1
    lines = genome.genome_str.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if '>' == line[0]:
            name, _ = splitext(genome.genome_name)
            lines[index] = f'>{name}_NODE_{count}\n'
            count += 1

    return ''.join(lines)


def _generate_multifasta(genomes: list[GenomeData]):
    pool = Pool(cpu_count())
    results_list = pool.map(_rename_fasta_headers, genomes)
    return ''.join(results_list)


def make_multifasta_for_blastdb(genomes: list[GenomeData], output: Path):
    """
    Join all the genomes in a single file, renaming the headers to make them unique.
    This multifasta file will be used to create the blast database.

    :param genomes:
    :param output:
    :return:
    """

    multifasta = _generate_multifasta(genomes)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(multifasta)


##### Run makeblastdb #####
def run_command(command, error_description=''):
    """
    This function is a wrapper to run a command in the console.
    It prints the output as soon as it's given, prints the error and raises a subprocess.CalledProcessError exception
    to keep track of where the error was found and which command returned it.

    """

    if not error_description:
        error_description = f'ERROR FOUND'

    with Popen(command, stdout=PIPE, stderr=PIPE, bufsize=1,
                          universal_newlines=True) as popen:

        # Save the errors in a variable to print them outside the with statement in case the return code is not zero.
        # It needs to be this way because inside the with statement stderr works but no return code is given,
        # and outside stderr does not work but there is a return code
        stderr = ''.join(popen.stderr)

    if popen.returncode != 0:
        print(f'\n{error_description}: ')
        text = ''
        for line in stderr:
            print(line, end='')
            text += line

        st.error(text)
        print(f'')
        raise CalledProcessError(popen.returncode, popen.args)


def make_blast_db(multifasta, database, makeblastdb_exec='makeblastdb', dbtype='nucl', title='blastdb'):
    """
    Make blast database
    """

    command = f"\"{makeblastdb_exec}\" -in \"{multifasta}\" -input_type fasta -parse_seqids -dbtype {dbtype} " \
              f"-out \"{database}\" -title \"{title}\""

    print(command)
    command = shlex.split(command)
    print(command)
    run_command(command, error_description='ERROR FOUND WHEN MAKING A BLAST DB')


##### Other functions #####
def get_program_path(program_name, binaries_in='BlastUI'):
    """
    This function returns the path to the program passed as an argument.
    It checks if the program is in the PATH variable, and if not, it checks if it is in the Binaries folder.
    If it is not in the Binaries folder, it raises an exception.
    """

    if binaries_in == 'PATH':
        program_path = Path(shutil.which(program_name))
        if not program_path.exists():
            st.warning(f'{program_name} not found in PATH. Using the one in the Binaries folder')
            binaries_in = 'BlastUI'

    if binaries_in == 'BlastUI':
        if sys.platform == "linux" or sys.platform == "linux2":
            program_path = Path(Path().cwd(), 'Binaries/linux/bin', program_name)
        elif sys.platform == "win32":
            program_path = Path(Path().cwd(), 'Binaries/win/bin', program_name + '.exe')

    if not program_path.exists():
        st.error(f'{program_name} not found in PATH or in the Binaries folder')
        st.stop()

    return program_path
