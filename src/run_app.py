from streamlit.web import cli
from multiprocessing import freeze_support


if __name__ == '__main__':
    # Needed when you use multiprocessing/concurrent.futures modules and you package the app into an executable
    freeze_support()

    # This function was created inside our streamlit framework -> streamlit.web.cli
    # cli._main_run_cl(str(Path(cwd, 'Home.py')), 'streamlit run')
    cli._main_run_cl('Home.py', 'streamlit run')
