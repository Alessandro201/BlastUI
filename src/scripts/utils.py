import shutil
import subprocess
import sys
from configparser import ConfigParser
from pathlib import Path
import streamlit as st
from io import BytesIO
import pandas as pd


class fragile(object):
    class Break(Exception):
        """Break out of the with statement"""

    def __init__(self, value):
        self.value = value

    def __enter__(self):
        return self.value.__enter__()

    def __exit__(self, etype, value, traceback):
        error = self.value.__exit__(etype, value, traceback)
        if etype == self.Break:
            return True
        return error


class GenomeData:
    def __init__(self, name: str, genome: str):
        self.name: str = name
        self.genome: str = genome


def read_configs(config_file='config.ini'):
    config = ConfigParser()
    config.read(config_file)
    return config


def write_configs(config, config_file='config.ini'):
    with open(config_file, 'w') as f:
        config.write(f)


def get_program_path(program_name, binaries_in='BlastUI'):
    """
    This function returns the path to the program passed as an argument.
    It checks if the program is in the PATH variable, and if not, it checks if it is in the Binaries folder.
    If it is not in the Binaries folder, it raises an exception.
    """

    if binaries_in == 'PATH':
        program_path = shutil.which(program_name)
        if program_path is None:
            st.warning(f'{program_name} not found in PATH. Using the one in the Binaries folder')
            binaries_in = 'BlastUI'

    if binaries_in == 'BlastUI':
        match sys.platform:
            case 'linux' | 'linux2' | 'darwin':
                program_path = Path('./Binaries/bin') / program_name
            case 'win32':
                program_path = Path('./Binaries/bin') / (program_name + '.exe')
            case _:
                raise OSError(f'Your platform ({sys.platform}) is not supported, there are no blast executables '
                              f'for your Operating System.')

    if not Path(program_path).exists():
        raise ValueError(f'{program_name} not found in the Binaries folder.')

    return program_path


def check_blast_executables_in_path():
    """
    Returns a list of paths to the blast executables found in $PATH.
    If a program is not found, it returns None.
    :return:
    """

    blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    for program in blast_programs:
        yield shutil.which(program)


def check_blast_executables_in_bin():
    """
    Returns a list of paths to the blast executables in the Binaries folder.
    If a program is not found, it returns None.
    :return:
    """

    blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    for program in blast_programs:
        program_path = Path('./Binaries/bin', program)
        if sys.platform == "win32":
            program_path = Path('./Binaries/bin', program + '.exe')

        yield program_path if program_path.exists() else None


def check_blast_version(program_path: Path) -> str:
    """
    Returns the version of the blast program passed as an argument.
    :param program_path: Path to the blast program
    :return: Version of the blast program
    """

    version = subprocess.run([program_path, '-version'], check=True, capture_output=True, text=True)
    return version.stdout.split()[1]


def run_command(command):
    subprocess.run(command, check=True, capture_output=True, text=True)


def generate_xlsx_table(df) -> bytes:
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine="xlsxwriter")

    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='Sheet1', index=False, engine='xlsxwriter')

    # Get the xlsxwriter workbook and worksheet objects.
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    # Get the dimensions of the dataframe.
    (max_row, max_col) = df.shape

    # Create a list of column headers, to use in add_table().
    column_settings = []
    for header in df.columns:
        column_settings.append({'header': header})

    # Add the table.
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings, 'style': 'Table Style Medium 11'})

    # Make the columns wider for clarity.
    worksheet.set_column(0, max_col - 1, 12)

    # Autofit columns
    worksheet.autofit()

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()
    output.seek(0)

    return output.getvalue()
