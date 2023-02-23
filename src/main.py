import sys
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT
import shlex

import webbrowser


def main():
    # Getting path to python executable (full path of deployed python on Windows)
    executable = sys.executable

    path_to_run = Path(Path(__file__).parent, "Home.py")

    # Running streamlit server in a subprocess and printing to console

    cmd = shlex.split(f'"{executable}" -m streamlit run "{path_to_run}"')

    with Popen(cmd, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True) as popen:
        print(f"Running streamlit server in a subprocess: {cmd}")

        # Force the opening (does not open automatically) of the browser tab after a brief delay to let
        # the streamlit server start.
        webbrowser.open("http://localhost:8501")

        # prints the output as soon as it's given
        for line in popen.stdout:
            print(line, end='', flush=True)
            sys.stdout.flush()

        # Save the errors in a variable to print them outside the with statement in case the return code is not zero.
        # It needs to be this way because inside the with statement stderr works but no return code is given,
        # and outside stderr does not work but there is a return code
        stderr = ''.join(popen.stderr)

    if popen.returncode != 0:
        error_description = f'ERROR IN RUNNING THE APPLICATION'

        print(f'\n{bcolors.FAIL}{error_description}: {bcolors.ENDC}')
        for line in stderr:
            print(bcolors.FAIL + line + bcolors.ENDC, end='')
        print('')
        raise subprocess.CalledProcessError(popen.returncode, popen.args)


if __name__ == "__main__":
    main()
