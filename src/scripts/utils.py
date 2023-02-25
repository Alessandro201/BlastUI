import sys
from configparser import ConfigParser
from pathlib import Path
from subprocess import Popen, PIPE, CalledProcessError


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
                program_path = Path().cwd() / 'Binaries/bin' / program_name
            case 'win32':
                program_path = Path().cwd() / 'Binaries/bin' / (program_name + '.exe')
            case _:
                raise OSError(f'Your platform ({platform}) is not supported.')

    if not Path(program_path).exists():
        st.error(f'{program_name} not found in PATH or in the Binaries folder')
        st.stop()

    return program_path


def check_blast_executables():
    blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    for program in blast_programs:
        yield program, shutil.which(program)


def check_blast_executables_in_bin():
    blast_programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    for program in blast_programs:
        program_path = Path('./Binaries/bin', program)
        if sys.platform == "win32":
            program_path = Path('./Binaries/bin', program + '.exe')
        yield program, program_path if program_path.exists() else None


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
        for line in stderr:
            print(line, end='')
        print('')
        raise CalledProcessError(popen.returncode, popen.args)
