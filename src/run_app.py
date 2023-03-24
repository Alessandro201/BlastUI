from streamlit.web import cli
from multiprocessing import freeze_support
import sys
from pathlib import Path


def set_icon():
    if getattr(sys, 'frozen', False):
        # running as a bundled executable
        icon_path = Path(sys._MEIPASS, 'icon.ico')  # assuming the icon file is packaged in your PyInstaller bundle
    else:
        # running as a script
        icon_path = Path('icon.ico')

    if not icon_path.exists():
        return

    # set the icon
    if icon_path:

        if sys.platform == 'win32':
            try:
                import ctypes

                # set an AppUserModelID if necessary
                ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('BlastUI')
                ctypes.windll.user32.SetClassLongW(
                    ctypes.windll.kernel32.GetConsoleWindow(),
                    -14,
                    ctypes.windll.shell32.SHGetFileInfoW(
                        icon_path,
                        0,
                        ctypes.byref(ctypes.wintypes.SHFILEINFOW()),
                        ctypes.sizeof(ctypes.wintypes.SHFILEINFOW()),
                        ctypes.windll.shell32.SHGFI_ICON |
                        ctypes.windll.shell32.SHGFI_SMALLICON |
                        ctypes.windll.shell32.SHGFI_USEFILEATTRIBUTES)
                    .hIcon)
            except:
                pass


if __name__ == '__main__':
    set_icon()

    # Needed when you use multiprocessing/concurrent.futures modules and you package the app into an executable
    freeze_support()

    # This function was created inside our streamlit framework -> streamlit.web.cli
    cli._main_run_cl('Home.py', 'streamlit run')
