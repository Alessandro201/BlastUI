import os
import shutil
import subprocess
import sys
from io import BytesIO
from pathlib import Path

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


def get_programs_path(blast_programs: list[str] = None) -> dict[str, Path] | None:
    """
    Returns a dictionary with the paths to the blast executables.
    It gives priority to the programs in the Binaries folder, and if they are not found, it looks for them in $PATH.
    All programs must be found in the same directory, otherwise it returns None to avoid having different versions of
    blast.
    """

    if not blast_programs:
        blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx',
                          'makeblastdb', 'makembindex', 'tblastn_vdb']

    blast_exec = dict()
    platform = sys.platform

    # Find programs in Binaries folder
    for program in blast_programs:
        match platform:
            case 'linux' | 'linux2' | 'darwin':
                program_path = Path('./Binaries/bin', program)
            case "win32":
                program_path = Path('./Binaries/bin', program + '.exe')
            case _:
                raise OSError(f'Your platform ({platform}) is not supported, there are no blast executables '
                              f'for your Operating System.')

        blast_exec[program] = program_path if program_path.exists() else None

    if all(blast_exec.values()):
        return blast_exec

    # Find programs in $PATH
    for program in blast_programs:
        program_path = shutil.which(program)

        if program_path is None:
            blast_exec[program] = None
            continue

        blast_exec[program] = program_path if Path(program_path).exists() else None

    if all(blast_exec.values()):
        return blast_exec

    return None


def get_blast_installation_directory():
    executables_in = None

    if all(_check_blast_executables_in_bin()):
        executables_in = './Binaries/bin'

    elif all(_check_blast_executables_in_path()):
        executables_in = '$PATH'
    else:
        pass

    return executables_in


def _check_blast_executables_in_bin():
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


def _check_blast_executables_in_path():
    """
    Returns a list of paths to the blast executables found in $PATH.
    If a program is not found, it returns None.
    :return:
    """

    blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    for program in blast_programs:
        yield shutil.which(program)


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


def resource_path(relative_path='.') -> Path:
    """ Get absolute path to resources, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return Path(base_path, relative_path)
