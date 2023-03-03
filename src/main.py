import shlex
import sys
import webbrowser
from pathlib import Path
from subprocess import Popen, PIPE, CalledProcessError


def main():
    # Getting path to python executable (full path of deployed python on Windows)
    executable = sys.executable

    path_to_run = Path(Path(__file__).parent, "Home.py")

    # Running streamlit server in a subprocess and printing to console

    cmd = shlex.split(f'"{executable}" -m streamlit run "{path_to_run}"')

    with Popen(cmd, stdout=sys.stdout, stderr=sys.stderr) as popen:
        print(f"Running streamlit server in a subprocess: {cmd}")

        # Force the opening (does not open automatically) of the browser tab after a brief delay to let
        # the streamlit server start.
        webbrowser.open("http://localhost:8501")


if __name__ == "__main__":
    main()
