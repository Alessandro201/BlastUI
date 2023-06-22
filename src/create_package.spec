# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path
from sys import platform

import altair
import streamlit
from PyInstaller.building.api import COLLECT, EXE, PYZ
from PyInstaller.building.build_main import Analysis
from PyInstaller.utils.hooks import collect_data_files

APP_NAME = 'BlastUI'
VERSION = '1.0.0'

datas = []
datas += collect_data_files('st_aggrid')
datas += collect_data_files('streamlit_option_menu')
datas += collect_data_files('bokeh')
datas += collect_data_files('st_keyup')
datas += [
    (Path(altair.__path__[0], "vegalite/v4/schema/vega-lite-schema.json"), "./altair/vegalite/v4/schema/"),
    (Path(streamlit.__path__[0], "static"), "./streamlit/static"),
    (Path(streamlit.__path__[0], "runtime"), "./streamlit/runtime"),
]
datas += [
    ('./media', './media'),
    ("./icon.png", "."),
]

hidden_imports = [
    "streamlit_extras.switch_page_button",
    "streamlit_extras.stoggle",
    "st_aggrid",
    "streamlit_option_menu",
    "streamlit_extras.no_default_selectbox",
    "xlsxwriter",
    "pyarrow.vendored.version",
    "st_keyup"
]

# Add all python files in the pages and scripts folders to hidden imports so that PyInstaller can find any module
# that is imported in those files
scripts = ['./Home.py']
for folder in ('./pages', './scripts'):
    for file in Path(folder).glob('**/*.py'):
        scripts.append(str(file.relative_to('.')))

block_cipher = None

a = Analysis(
    ['run_app.py'] + scripts,
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=['./hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=scripts,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Remove the scripts from the Analysis object so that they are not included in the EXE
scripts_name = [Path(script).name for script in scripts]
for item in a.scripts:
    if item[0] in scripts_name:
        print(item)
        a.scripts.remove(item)

# The EXE contain the Python interpreter, the libraries and other files but the scripts need to be copied later
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [('icon.png', Path('icon.png').resolve(), 'DATA')],
          name=APP_NAME,
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True,
          icon='icon.png')

# Add the scripts to the COLLECT object so that they are copied to the dist folder
items_to_collect = list()
for file in scripts:
    items_to_collect.append((file, Path(file).resolve(), 'DATA'))
items_to_collect.append(('./.streamlit/config.toml', Path('./.streamlit/config.toml').resolve(), 'DATA'))

coll = COLLECT(exe,
               items_to_collect,
               strip=False,
               upx=True,
               upx_exclude=[],
               name=f"{APP_NAME}_{platform}_v{VERSION}",
               )

if platform == 'win32':
    # Remove the .exe file from the dist folder
    Path(f'dist/{APP_NAME}.exe').unlink()
elif platform == 'linux':
    # Remove the executable file from the dist folder
    Path(f'dist/{APP_NAME}').unlink()
