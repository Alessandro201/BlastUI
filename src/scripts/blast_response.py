from pathlib import Path
from typing import Union, Optional, Iterable, Generator
import streamlit as st
import pandas as pd
import json


class BlastMatch:
    def __init__(self, **kwargs: dict):
        self.program: str = kwargs.get('program', '')

        self.id = kwargs.get('id', '')
        self.query_title: str = kwargs.get('query_title', '')
        self.query_len: int = kwargs.get('query_len', 0)
        self.strain: str = kwargs.get('strain', '')
        self.node: str = kwargs.get('node', '')

        self.qseq: str = kwargs.get('qseq', '')
        self.midline: str = kwargs.get('midline', '')
        self.sseq: str = kwargs.get('hseq', '')

        self.align_len: int = kwargs.get('align_len', 0)
        self.identity: int = kwargs.get('identity', 0)

        self.mismatch: int = self.align_len - self.identity
        self.q_gap = self.qseq.count('-')
        self.s_gap: int = kwargs.get('gaps', 0)

        self.perc_identity = round(self.identity / self.align_len * 100)
        self.perc_alignment = round((self.align_len + self.q_gap) / self.query_len * 100)
        self.perc_gaps = round(self.s_gap / self.align_len * 100)
        self.perc_mismatches = round(self.mismatch / self.align_len * 100)

        self.q_start: int = kwargs.get('query_from', 0)
        self.q_end: int = kwargs.get('query_to', 0)
        self.s_start: int = kwargs.get('hit_from', 0)
        self.s_end: int = kwargs.get('hit_to', 0)

        self.evalue: float = kwargs.get('evalue', 0.0)
        self.bit_score: float = kwargs.get('bit_score', 0.0)
        self.score: int = kwargs.get('score', 0)

        if self.program in ('tblastn', 'blastx', 'blastp', 'tblastx'):
            self.positive: int = kwargs.get('positive', 0)
            self.perc_positives = round(self.positive / self.align_len * 100)

        match self.program:
            case 'blastn':
                self.q_orient: str = kwargs.get('query_strand', '')
                self.s_orient: str = kwargs.get('hit_strand', '')
            case 'tblastn' | 'tblastx':
                self.s_frame = kwargs.get('hit_frame', 0)
            case 'blastx':
                self.q_frame = kwargs.get('query_frame', 0)
            case 'blastp' | _:
                pass

        self._test_args()

    def alignment(self) -> str:
        match self.program:
            case 'tblastn' | 'tblastx':
                alignment_text = f">{self.query_title} \n" \
                                 f"Strain = {self.strain}, Node = {self.node}\n" \
                                 f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                                 f"E-value = {self.evalue :.1g} \n" \
                                 f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                                 f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%) \n" \
                                 f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                                 f"\tMismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                                 f"Gaps = {self.s_gap}/{self.align_len} ({self.perc_gaps}%)\n" \
                                 f"\tFrame = {self.s_frame}\n\n"

            case 'blastn':
                alignment_text = f">{self.query_title} \n" \
                                 f"Strain = {self.strain}, Node = {self.node}\n" \
                                 f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                                 f"E-value = {self.evalue :.1g} \n" \
                                 f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                                 f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%) \n" \
                                 f"\tMismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                                 f"Gaps = {self.s_gap}/{self.align_len} ({self.perc_gaps}%)\n" \
                                 f"\tStrand = {self.q_orient}/{self.s_orient}\n\n"

            case 'blastp':
                alignment_text = f">{self.query_title} \n" \
                                 f"Strain = {self.strain}, Node = {self.node}\n" \
                                 f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                                 f"E-value = {self.evalue :.1g} \n" \
                                 f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                                 f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%) \n" \
                                 f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                                 f"\tMismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                                 f"Gaps = {self.s_gap}/{self.align_len} ({self.perc_gaps}%)\n\n"

            case 'blastx':
                alignment_text = f">{self.query_title} \n" \
                                 f"Strain = {self.strain}, Node = {self.node}\n" \
                                 f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                                 f"E-value = {self.evalue :.1g} \n" \
                                 f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                                 f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%) \n" \
                                 f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                                 f"\tMismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                                 f"Gaps = {self.s_gap}/{self.align_len} ({self.perc_gaps}%)\n" \
                                 f"\tQuery frame = {self.q_frame}\n\n"

            case _:
                alignment_text = f">{self.query_title} \n"

        # If the query or the sequence are translated, their positions in the alignment are multiplied by 3
        s_multiplier = 1
        q_multiplier = 1
        match self.program:
            case 'tblastn':
                s_multiplier = 3

            case 'blastx':
                q_multiplier = 3

            case 'tblastx':
                s_multiplier = 3
                q_multiplier = 3

            case 'blastp' | 'blastn' | _:
                s_multiplier = 1
                q_multiplier = 1

        # The number of padding spaces for the indexes in the alignment depends on the number of digits
        pad = max(len(str(self.s_start)),
                  len(str(self.q_start)),
                  len(str(self.s_end)),
                  len(str(self.q_end)))

        index = 0
        query_gap = 0
        seq_gap = 0
        prev_query_gap = 0
        prev_seq_gap = 0

        while index < len(self.sseq):
            query_gap += self.qseq[index:index + 60].count('-')
            seq_gap += self.sseq[index:index + 60].count('-')

            # Query
            if self.q_start < self.q_end:
                # Forward strand
                q_start = self.q_start + (index - prev_query_gap) * q_multiplier
                q_end = self.q_start + (index + 60 - query_gap) * q_multiplier - 1
                q_end = min(q_end, self.q_end)

            else:
                # Reverse strand
                q_start = self.q_end - (index - prev_query_gap) * q_multiplier
                q_end = self.q_end - (index + 60 - query_gap) * q_multiplier + 1
                q_end = max(q_end, self.q_start)

            # Subject
            if self.s_start < self.s_end:
                # Forward strand
                s_start = self.s_start + (index - prev_seq_gap) * s_multiplier
                s_end = self.s_start + (index + 60 - seq_gap) * s_multiplier - 1
                s_end = min(s_end, self.s_end)

            else:
                # Reverse strand
                s_start = self.s_start - (index - prev_seq_gap) * s_multiplier
                s_end = self.s_start - (index + 60 - seq_gap) * s_multiplier + 1
                s_end = max(s_end, self.s_end)

            prev_query_gap = query_gap
            prev_seq_gap = seq_gap

            alignment_text += f"Query  {q_start :{pad}}  {self.qseq[index:index + 60]}  {q_end :{pad}}\n" \
                              f"       {' ' * pad}  {self.midline[index:index + 60]} \n" \
                              f"Sbjct  {s_start :{pad}}  {self.sseq[index:index + 60]}  {s_end :{pad}}\n\n"
            index += 60

        return alignment_text

    def _test_args(self):
        assert 0 <= self.perc_identity <= 100, f"Identity value out of range [0:100]: {self.perc_identity}"
        assert self.align_len >= 1, f"length value out of range [>= 1]: {self.align_len}"
        assert self.mismatch >= 0, f"mis value out of range [>= 0]: {self.mismatch}"
        assert self.s_gap >= 0, f"gap value out of range [>= 0]: {self.s_gap}"
        assert self.q_start >= 0, f"q_start value out of range [>= 1]: {self.q_start}"
        assert self.q_end >= 0, f"q_end value out of range [>= 1]: {self.q_end}"
        assert self.s_start >= 0, f"s_start value out of range [>= 1]: {self.s_start}"
        assert self.s_end >= 0, f"s_end value out of range [>= 1]: {self.s_end}"
        assert self.evalue >= 0, f"evalue value out of range [>= 0]: {self.evalue}"
        assert self.bit_score >= 0, f"bit_score value out of range [>= 0]: {self.bit_score}"


class BlastResponse:
    headers = {
        'tblastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'positive', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                    's_end', 'evalue', 'bit_score', 'sseq', 'id'],
        'blastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                   's_end', 'evalue', 'bit_score', 'sseq', 'id'],
        'blastp': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'positive', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                   's_end', 'evalue', 'bit_score', 'sseq', 'id'],
        'blastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                   's_end', 'evalue', 'bit_score', 'sseq', 'id'],
        'tblastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'mismatch', 'gap_opens', 'q_start', 'q_end', 's_start',
                    's_end', 'evalue', 'bit_score', 'sseq', 'id']
    }

    def __init__(self, json_file: Union[dict, Path]):
        self.json: dict = self._read(json_file)
        self.reports = [report['report'] for report in self.json['BlastOutput2']]
        self.program = self.reports[0]['program']
        self.matches: dict = dict()
        self.parse_json()
        self.df = self.generate_df()
        self.reindexed_df = self.reindex_df()

    def _read(self, json_file: Union[dict, Path]):
        if isinstance(json_file, Path):
            with open(json_file, 'r') as f:
                return json.load(f)
        elif isinstance(json_file, dict):
            return json_file
        else:
            raise TypeError(f"Expected a dict or Path object, got {type(json_file)} instead")

    def parse_json(self):
        index = 0
        for report in self.reports:
            # If you don't insert a header in the query sequence, there won't be a query_title
            try:
                query_title: str = report['results']['search']['query_title']
            except KeyError:
                query_title: str = report['results']['search']['query_id']

            query_len: int = report['results']['search']['query_len']

            if "message" in report['results']['search']:
                st.warning(f"BLAST returned the following message: {report['results']['search']['message']}")
                st.stop()

            for hits in report['results']['search']['hits']:
                for hit in hits['hsps']:
                    hit['id'] = index
                    hit['query_title'] = query_title
                    hit['query_len'] = query_len
                    strain_node = hits['description'][0]['accession']
                    hit['strain'], hit['node'] = strain_node.split('_NODE_')
                    hit['program'] = self.program

                    self.matches[index] = BlastMatch(**hit)
                    index += 1

    def alignments(self, indexes: Optional[list[int]] = None) -> Generator[str, None, None]:
        if indexes is None:
            indexes = range(len(self.matches))

        for index in indexes:
            yield self.matches[index].alignment()

    def generate_df(self) -> pd.DataFrame:
        df = pd.DataFrame().from_records([match.__dict__ for match in self.matches.values()])
        df.sort_values(by=['perc_identity', 'perc_alignment'], ascending=[False, False], inplace=True)

        return df

    def reindex_df(self, columns: Optional[Iterable[str]] = None) -> pd.DataFrame:
        """
        Return a copy of the original dataframe with the columns reordered to match the columns given,
        or the default columns for the program used if no columns are given.
        """

        if not columns:
            columns = self.headers[self.program]

        df = self.df
        columns_to_drop = set(self.df.columns) - set(columns)
        df = df.drop(columns=columns_to_drop)
        df = df.reindex(columns=columns)
        # df.insert(0, 'row', pd.array(range(1, len(df))))

        return df

    def filtered_df(self, identity: float = 0, query_cov: float = 0) -> pd.DataFrame:
        """
        Filter the results of the blast analysis to keep only the results with query_cov and identity higher
        than the ones provided.
        """

        if not 0 <= identity <= 100:
            raise ValueError(f"Identity must be between 0 and 100, not {identity}")

        if not 0 <= query_cov <= 100:
            raise ValueError(f"Query coverage must be between 0 and 100, not {query_cov}")

        df = self.reindexed_df
        return df[(df['perc_identity'] >= identity) & (df['perc_alignment'] >= query_cov)]
