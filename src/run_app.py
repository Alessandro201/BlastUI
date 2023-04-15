from multiprocessing import freeze_support
from typing import Any, Dict, List, Optional

import streamlit.web.bootstrap as bootstrap
from streamlit.web.cli import check_credentials


# This is a modified version of a function inside streamlit.web.cli to make streamlit run from within a file
# This way it works after having been packaged by PyInstaller
def run_streamlit(
        file,
        command_line: str,
        args: Optional[List[str]] = None,
        flag_options: Optional[Dict[str, Any]] = None,
) -> None:
    if args is None:
        args = []

    if flag_options is None:
        flag_options = {}

    check_credentials()

    bootstrap.run(file, command_line, args, flag_options)


if __name__ == '__main__':
    # The first time that streamlit is run it print a hello message which has an emoji.
    # If the terminal is not configured to support emojis, it will throw an error, thus we need to configure it
    # to use utf-8 encoding.
    import sys

    sys.stdout.reconfigure(encoding='utf-8')
    sys.stdin.reconfigure(encoding='utf-8')

    # Needed when you use multiprocessing/concurrent.futures modules and you package the app into an executable
    freeze_support()

    # This function was created inside our streamlit framework -> streamlit.web.cli
    # cli._main_run_cl(str(Path(cwd, 'Home.py')), 'streamlit run')
    run_streamlit('Home.py', 'streamlit run')
