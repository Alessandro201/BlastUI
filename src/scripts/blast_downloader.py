import hashlib
import os
import re
import shutil
import stat
import sys
import tarfile
import urllib.request
from contextlib import closing
from pathlib import Path

import streamlit as st


def sizeof_fmt(num, suffix="B"):
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"


def md5(fname):
    """
    Calculate the md5 hash of a file and return it as hexadecimal.

    :param fname: file of which to calculate the md5 hash
    :return: hexadecimal hash
    """

    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


class BlastDownloader:
    def __init__(self,
                 pbar=None):

        self.URL = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'

        download_file, download_file_md5 = self.get_download_files_names()
        self.filename, self.filesize = download_file
        self.md5filename, self.md5filesize = download_file_md5

        # Delete the ./bin folder to remove any old file inside it and recreate it.
        self.download_folder = Path().cwd() / 'bin'
        shutil.rmtree(self.download_folder, ignore_errors=True)
        self.download_folder.mkdir(parents=True, exist_ok=True)

        self.download_path = self.download_folder / self.filename
        self.md5_download_path = self.download_folder / self.md5filename

        # Letting the user give a progress bar is useful because he can place it where it wants while
        # letting this class update it.
        if not pbar:
            pbar = st.progress(0)
        self.pbar = pbar

        self.downloaded_bytes = 0

        self.download()
        self.extract_bin()
        self.remove_unnecessary_executables()

        self.download_path.unlink()
        self.md5_download_path.unlink()

    def get_download_files_names(self):
        match platform := sys.platform:
            case 'linux' | 'linux2':
                os_name = 'linux'
            case 'win32':
                os_name = 'win64'
            case 'darwin':
                os_name = 'macosx'
            case _:
                raise OSError(f'Your platform ({platform}) is not supported.')

        with closing(urllib.request.urlopen(self.URL)) as r:
            page_text = r.read().decode('utf-8')
            matches = re.findall(r'<a href="(.+?)">.+?</a>', page_text)

        download_file = None
        download_file_md5 = None
        for file in matches:
            if file.endswith(f'{os_name}.tar.gz'):
                with closing(urllib.request.urlopen(self.URL + file)) as r:
                    size = r.length
                download_file = (file, size)

            elif file.endswith(f'{os_name}.tar.gz.md5'):
                with closing(urllib.request.urlopen(self.URL + file)) as r:
                    size = r.length
                download_file_md5 = (file, size)

        if not download_file and not download_file_md5:
            raise FileNotFoundError(f'Could not find BLAST for your platform ({platform}). '
                                    f'Please install it manually, either by downloading it from NCBI and placing '
                                    f'the executables inside the "BlastUI/src/bin" folder or by installing it '
                                    f'in the $PATH. If you can use BLAST from the terminal, the program should '
                                    f'be able to use it as well.')

        return download_file, download_file_md5

    def _update_progress(self, data):
        self.downloaded_bytes += len(data)
        percentage = round((self.downloaded_bytes / self.filesize) * 100)
        percentage = min(percentage, 100)
        self.pbar.progress(percentage, text=f'Downloading:    {self.filename}...      '
                                            f'{sizeof_fmt(self.downloaded_bytes)}/{sizeof_fmt(self.filesize)} '
                                            f'({percentage}%)')

    @staticmethod
    def _download_file(download_url, download_path, buf_size=1024 * 1024, callback=None, *args, **kwargs):
        with closing(urllib.request.urlopen(download_url)) as r:
            with open(download_path, 'wb') as f:
                while True:
                    buf = r.read(buf_size)
                    if not buf:
                        break

                    if callback:
                        callback(buf, *args, **kwargs)

                    f.write(buf)

    def download(self):
        self.pbar.progress(0, text=f'Downloading:    md5 hash...')
        self._download_file(self.URL + self.md5filename, self.md5_download_path)
        self._download_file(self.URL + self.filename, self.download_path, buf_size=1024 * 256,
                            callback=self._update_progress)

        if os.stat(self.download_path).st_size != self.filesize:
            raise ValueError(f'The downloaded file was corrupted (size does not match). Try again. '
                             f'downloaded file: {os.stat(self.download_path).st_size}, '
                             f'original file: {self.filesize}')

        if os.stat(self.md5_download_path).st_size != self.md5filesize:
            raise ValueError(f'The downloaded file was corrupted (size does not match). Try again. '
                             f'downloaded file: {os.stat(self.md5_download_path).st_size}, '
                             f'original file: {self.md5filesize}')

        self.check_hash()

    def check_hash(self):
        md5_of_downloaded_file = md5(self.download_path)
        md5_correct = Path(self.md5_download_path).read_text().split(' ')[0]
        if md5_of_downloaded_file != md5_correct:
            raise ValueError(f'The downloaded file was corrupted (hashes do not match). Try again.')

    def extract_bin(self):
        self.pbar.progress(100, text=f'Extracting {self.filename}...')

        # Extracting the tar file.
        with tarfile.open(self.download_path, "r:gz") as tar:
            tar.extractall(path=self.download_folder)

        # Moving the bin folder to the parent directory.
        bin_dir = self.download_folder / self.filename.split('-x64')[0] / 'bin'

        for file in bin_dir.iterdir():
            shutil.move(file, self.download_folder)
        shutil.rmtree(bin_dir.parent)

        self.pbar.progress(100, text=f'Done!')

    def remove_unnecessary_executables(self):
        for file in self.download_folder.iterdir():
            try:
                if file.suffix == '.dll':
                    continue
                if file.suffix == '.manifest' or file.suffix == '.pl' or file.suffix == '.py':
                    file.unlink()
                elif file.name.startswith(('rps', 'windowmasker', 'segmasker', 'psiblast', 'dustmasker',
                                           'blastdbcheck', 'blastdbcmd', 'blastdb_aliastool', 'blast_formatter',
                                           'blastn_vdb', 'blast_vdb_cmd', 'cleanup-blastdb-volumes', 'deltablast',
                                           'makeprofiledb', 'convert2blastmask', 'get_species_taxids.sh')):
                    file.unlink()

            except PermissionError:
                # Sometimes some folders get a stubborn read-only attribute which inhibits
                # os.rmdir from removing the directory. Changing the file permission should
                # do the trick

                # Change folder permissions to 0777:
                # stat.S_IRWXU Mask for file owner permissions.
                # stat.S_IRWXG Mask for group permissions.
                # stat.S_IRWXO Mask for permissions for others (not in group).
                os.chmod(file, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
                try:
                    file.unlink()
                except PermissionError:
                    pass
