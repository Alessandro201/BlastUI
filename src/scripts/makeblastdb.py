from concurrent.futures import ProcessPoolExecutor, as_completed

from os.path import splitext

import shlex
import streamlit as st
from pathlib import Path

from scripts.utils import *


class MakeBlastDB:
    def __init__(self,
                 genomes: list[GenomeData],
                 db_name: Path,
                 makeblastdb_exec: Path = 'makeblastdb',
                 threads: int = 1,
                 pbar=None,
                 multifasta: Path = None,
                 dbtype: str = 'nucl',
                 rename_headers: bool = False):

        self.genomes = genomes
        self.threads = threads

        self.rename_headers = rename_headers
        self.pbar = pbar

        self.multifasta = Path(multifasta) if multifasta else Path('./BlastDatabases/multifasta.fasta')
        self.multifasta.parent.mkdir(parents=True, exist_ok=True)
        if self.multifasta.exists():
            self.multifasta.unlink()

        self.db = Path('./BlastDatabases', db_name, 'blastdb')

        self.exec = makeblastdb_exec
        self.dbtype = dbtype

    @staticmethod
    def _worker_rem_contigs(genome: GenomeData, min_length: int) -> GenomeData:
        """
        Remove all small contigs from a file and rewrites it.
        """

        # Split into contigs, Don't keep the first empty sequence
        contigs = genome.genome.split('>')[1:]

        indexes_to_remove = list()
        for index, contig in enumerate(contigs):

            header, seq = contig.split('\n', maxsplit=1)
            n_lines = seq.count('\n')

            if (len(seq) - n_lines) < min_length:
                indexes_to_remove.append(index)

        for index in reversed(indexes_to_remove):
            contigs.pop(index)

        if not contigs:
            raise ValueError(f'No contigs left after filtering for {genome.name}')

        return GenomeData(name=genome.name, genome='>' + '>'.join(contigs))

    def remove_small_contigs(self, min_length):
        """
        :param min_length: threshold length of the contigs
        :return:
        """

        filtered_genomes = list()

        self.pbar.progress(0, text=f"Removing small contigs... (0/{len(self.genomes)})")

        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = [executor.submit(self._worker_rem_contigs, genome, min_length) for genome in self.genomes]

            for i, future in enumerate(as_completed(futures)):
                if self.pbar:
                    percentage = round((i + 1) / len(futures) * 100)
                    self.pbar.progress(percentage, text=f"Removing small contigs... ({i + 1}/{len(futures)})")
                filtered_genomes.append(future.result())

        self.genomes = filtered_genomes

    @staticmethod
    def _rename_fasta_headers(genome: GenomeData) -> str:
        """
        Read the fasta file, renaming each scaffold with the name of the file and a counter.
        It's both for clarity (1) and a bit more security (2).

        1) When you use blast it will output the name of the contig in which a match was found, but the assembler I
        used, SPAdes, does not use any meaningful name for the scaffolds. This is a way to know straight away the genome
        which gave a match.
        2) It happened that multiple genomes had the same prefix for the contigs. If you use genomes from the same folder
        and with the same extension you are sure that you will not get any duplicate name.

        :param genome: GenomeData object containing the genome name and the genome as string
        :return:
        """

        name, _ = splitext(genome.name)

        # Split into contigs, Don't keep the first empty sequence
        contigs = genome.genome.split('>')[1:]

        for index, contig in enumerate(contigs):
            header, seq = contig.split('\n', maxsplit=1)
            header = f'>{name}_NODE_{index + 1}\n'
            contigs[index] = header + seq

        return ''.join(contigs)

    def generate_multifasta(self):

        if not self.rename_headers:
            self.multifasta.write_text(''.join([genome.genome for genome in self.genomes]))
            return

        self.pbar.progress(0, text=f"Joining the genomes... (0/{len(self.genomes)})")

        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = [executor.submit(self._rename_fasta_headers, genome) for genome in self.genomes]

            for i, future in enumerate(as_completed(futures)):
                if self.pbar:
                    percentage = round((i + 1) / len(futures) * 100)
                    self.pbar.progress(percentage, text=f"Joining the genomes... ({i + 1}/{len(futures)})")

                with open(self.multifasta, 'a') as f:
                    f.write(future.result())

    def run(self):
        """
        Make blast database
        """

        command = f'"{self.exec}" -in "{self.multifasta}" -input_type fasta -parse_seqids -dbtype {self.dbtype} ' \
                  f'-out "{self.db}" -title "{self.db.parent.name}" -hash_index -max_file_sz 2GB'

        command = shlex.split(command)
        try:
            run_command(command)
        finally:
            self.multifasta.unlink()
