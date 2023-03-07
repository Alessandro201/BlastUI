import ftplib
import hashlib
import shutil
import streamlit as st
import sys
import os
import tarfile
from ftplib import FTP
import streamlit as st
from pathlib import Path
import stat


class DownloadError(Exception):
    pass


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
    def __init__(self, host: str = None,
                 directory_url: str = None,
                 download_folder: Path = None,
                 force_download=False,
                 pbar=None):

        if not host:
            host = 'ftp.ncbi.nlm.nih.gov'
        if not directory_url:
            directory_url = 'blast/executables/blast+/LATEST'

        self.force_download = force_download

        try:
            self.ftp = FTP(host, encoding='latin-1')
            self.ftp.login()
            self.ftp.cwd(directory_url)
            self.filename, self.filesize, self.md5filename = self.get_download_file()
            self.filesize = int(self.filesize)
            self.downloaded_bytes = 0

            if not download_folder:
                download_folder = Path().cwd() / 'Binaries' / self.filename
            else:
                if Path(download_folder).is_file():
                    raise OSError(f'{download_folder} is not a file')

                download_folder = Path(download_folder) / self.filename

            self.download_path = Path(download_folder)
            self.download_path.parent.mkdir(parents=True, exist_ok=True)
            self.md5_download_path = self.download_path.parent / self.md5filename

            # Letting the user give a progress bar is useful because he can place it where it wants while
            # letting this class update it.
            if not pbar:
                pbar = st.progress(0)
            self.pbar = pbar.progress(0, text=f'Downloading {self.filename}...')

            if self.download_path.exists():
                if self.force_download:
                    self.download_path.unlink()
                else:
                    raise OSError(f'{self.filename} already downloaded')

            if self.md5_download_path.exists():
                if self.force_download:
                    self.md5_download_path.unlink()
                else:
                    raise OSError(f'{self.md5filename} already downloaded')

            self.file_handle = None
            self.download()

        finally:
            self.ftp.quit()

        self.extract_bin()
        self.remove_unnecessary_executables()

        self.download_path.unlink()
        self.md5_download_path.unlink()

    def get_download_file(self):
        for file in self.ftp.mlsd():
            file_name = file[0]
            file_type = file[1]['type']
            file_size = file[1]['size']

            platform = sys.platform
            if platform == "linux" or platform == "linux2":
                os_name = 'linux'
            elif platform == "win32":
                os_name = 'win64'
            elif platform == 'darwin':
                os_name = 'macosx'
            else:
                raise OSError(f'Your platform ({platform}) is not supported.')

            if not file_type == 'file':
                continue

            if file_name.endswith(f'{os_name}.tar.gz'):
                return file_name, file_size, file_name + '.md5'

    def _write_file(self, data):
        self.file_handle.write(data)

        self.downloaded_bytes += len(data)
        percentage = round((self.downloaded_bytes / self.filesize) * 100)
        percentage = min(percentage, 100)
        self.pbar.progress(percentage, text=f'Downloading:    {self.filename}...      '
                                            f'{sizeof_fmt(self.downloaded_bytes)}/{sizeof_fmt(self.filesize)} '
                                            f'({percentage}%)')

    def download(self):
        self.file_handle = open(self.download_path, 'wb')
        try:
            self.ftp.retrbinary("RETR " + self.filename, self._write_file, blocksize=8*1024)
            self.ftp.retrbinary("RETR " + self.md5filename, Path(self.md5_download_path).write_bytes)
        finally:
            self.file_handle.close()

        self.check_hash()

    def check_hash(self):
        md5_downloaded_blast = md5(self.download_path)
        md5_blast = Path(self.md5_download_path).read_text().split(' ')[0]
        if md5_downloaded_blast != md5_blast:
            raise DownloadError(f'The downloaded file was corrupted (hashes do not match). Try again.')

    def extract_bin(self):
        self.pbar.progress(100, text=f'Extracting {self.filename}...')

        # If the user wants to force the download, we need to remove the old blast folder and bin folder.
        blast_dir = self.download_path.parent / self.filename.split('-x64')[0]
        if blast_dir.exists() and self.force_download:
            shutil.rmtree(blast_dir)

        bin_dir = self.download_path.parent / 'bin'
        if bin_dir.exists() and self.force_download:
            shutil.rmtree(bin_dir)

        # Extracting the tar file.
        with tarfile.open(self.download_path, "r:gz") as tar:
            tar.extractall(path=self.download_path.parent)

        # Moving the bin folder to the parent directory.
        bin_dir = self.download_path.parent / self.filename.split('-x64')[0] / 'bin'
        shutil.move(bin_dir, self.download_path.parent)
        shutil.rmtree(bin_dir.parent)

        self.pbar.progress(100, text=f'Done!')

    def remove_unnecessary_executables(self):
        bin_dir = self.download_path.parent / 'bin'
        for file in bin_dir.iterdir():
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
                file.unlink()
