import re
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd
import regex_spm
import streamlit as st
from bokeh.models import (ColumnDataSource, HoverTool)
from bokeh.models import Legend, Rect
from bokeh.plotting import figure

from scripts.substitution_matrix import get_matrix, get_substitution_score


class EmptyCSVError(Exception):
    pass


class Alignment:
    def __init__(self, row, program, subtitution_matrix='BLOSUM62'):
        self.program = program.lower()

        if self.program not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
            raise ValueError(f"Invalid program: {self.program}")

        self.query_title: str = row['query_title']
        self.query_len: int = row['query_len']
        self.strain: str = row['strain']
        self.node: int = int(row['node'])

        self.align_len: int = row['align_len']
        self.identity: int = row['identity']

        self.gaps: int = row['gaps']
        self.gap_opens: int = row['gap_open']
        self.mismatch: int = row['mismatch']

        self.perc_identity: float = row['perc_identity']
        self.perc_alignment: float = row['perc_alignment']
        self.perc_gaps: float = np.round(self.gaps / self.align_len * 100)
        self.perc_mismatch = np.round(self.mismatch / (self.align_len - self.gaps) * 100)

        self.q_start: int = row['query_start']
        self.q_end: int = row['query_end']
        self.s_start: int = row['seq_start']
        self.s_end: int = row['seq_end']

        self.evalue: float = row['evalue']
        self.bit_score: float = row['bit_score']
        self.score: int = row['score']

        if self.program in ('tblastn', 'blastx', 'blastp', 'tblastx'):
            self.matrix = get_matrix(subtitution_matrix)
            self.positive: int = row['positive']
            self.perc_positives: float = round(self.positive / self.align_len * 100)

        self.qseq: str = row['qseq']
        self.sseq: str = row['sseq']
        self.midline: str = self._get_midline()

        # Find orientation of the strands
        self.q_frame = row['query_frame']
        self.s_frame = row['seq_frame']
        self.q_orient: str = "forward" if int(self.q_frame) >= 0 else "reverse"
        self.s_orient: str = "forward" if int(self.s_frame) >= 0 else "reverse"

        # If the query or the sequence are translated, their positions in the alignment are multiplied by 3
        if self.program == 'tblastn':
            self.q_multiplier = 1
            self.s_multiplier = 3
        elif self.program == 'blastx':
            self.q_multiplier = 3
            self.s_multiplier = 1
        elif self.program == 'tblastx':
            self.s_multiplier = 3
            self.q_multiplier = 3
        elif self.program in ('blastp', 'blastn'):
            self.s_multiplier = 1
            self.q_multiplier = 1

        # The number of padding spaces for the indexes in the alignment depends on the number of digits
        self.pad = max(len(str(self.s_start)),
                       len(str(self.q_start)),
                       len(str(self.s_end)),
                       len(str(self.q_end)))

    def align(self):
        match self.program:
            case 'blastn':
                return self._blastn_header() + self._alignment()

            case 'blastp':
                return self._blastp_header() + self._alignment()

            case 'blastx':
                return self._blastx_header() + self._alignment()

            case 'tblastn':
                return self._tblastn_header() + self._alignment()

            case 'tblastx':
                return self._tblastx_header() + self._alignment()

    def _tblastn_header(self):
        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatch}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tFrame = {self.s_frame}\n\n"

        return alignment_text

    def _blastn_header(self):
        q_strand = "Plus" if self.q_orient == 'forward' else "Minus"
        s_strand = "Plus" if self.s_orient == 'forward' else "Minus"

        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatch}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tStrand = {q_strand}/{s_strand}\n\n"

        return alignment_text

    def _blastx_header(self):
        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatch}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tQuery frame = {self.q_frame}\n\n"

        return alignment_text

    def _blastp_header(self):
        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatch}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n\n"

        return alignment_text

    def _tblastx_header(self):
        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"\tPositives = {self.positive}/{self.align_len} ({self.perc_positives}%), " \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatch}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tFrame = {self.q_frame}/{self.s_frame}\n\n"

        return alignment_text

    def _alignment(self, alignment_text: str = ''):
        # position in the sequences relative to the start of the alignment
        index = 0

        # Gaps in the current portion of the alignment
        query_gaps = 0
        seq_gaps = 0

        # Gaps from the start of the alignment
        prev_query_gaps = 0
        prev_seq_gaps = 0

        while index < len(self.sseq):
            query_gaps += self.qseq[index:index + 60].count('-')
            seq_gaps += self.sseq[index:index + 60].count('-')

            # Query
            if self.q_orient == 'forward':
                # Forward strand or positive frame
                q_start = self.q_start + (index - prev_query_gaps) * self.q_multiplier
                q_end = self.q_start + (index + 60 - query_gaps) * self.q_multiplier - 1
                q_end = min(q_end, self.q_end)

            else:
                # Reverse strand or negative frame
                q_start = self.q_start - (index - prev_query_gaps) * self.q_multiplier
                q_end = self.q_start - (index + 60 - query_gaps) * self.q_multiplier + 1
                q_end = max(q_end, self.q_end)

            # Subject
            if self.s_orient == 'forward':
                # Forward strand or positive frame
                s_start = self.s_start + (index - prev_seq_gaps) * self.s_multiplier
                s_end = self.s_start + (index + 60 - seq_gaps) * self.s_multiplier - 1
                s_end = min(s_end, self.s_end)

            else:
                # Reverse strand or negative frame
                s_start = self.s_start - (index - prev_seq_gaps) * self.s_multiplier
                s_end = self.s_start - (index + 60 - seq_gaps) * self.s_multiplier + 1
                s_end = max(s_end, self.s_end)

            prev_query_gaps = query_gaps
            prev_seq_gaps = seq_gaps

            alignment_text += f"Query  {q_start :{self.pad}}  {self.qseq[index:index + 60]}  {q_end :{self.pad}}\n" \
                              f"       {' ' * self.pad}  {self.midline[index:index + 60]} \n" \
                              f"Sbjct  {s_start :{self.pad}}  {self.sseq[index:index + 60]}  {s_end :{self.pad}}\n\n"
            index += 60

        return alignment_text

    def _get_midline(self):
        """
        Get midline between query and alignment to know which residues are the same and which are different.
        """

        midline = []

        if self.program == 'blastn':
            for q, s in zip(self.qseq.upper(), self.sseq.upper()):
                if q == s:
                    midline.append('|')
                else:
                    # Gap in either sequences or mismatch
                    midline.append(' ')
        else:
            for q, s in zip(self.qseq.upper(), self.sseq.upper()):
                if q == s:
                    midline.append(q)
                elif q == '-' or s == '-':
                    # Gap in either sequences
                    midline.append(' ')
                else:
                    score = get_substitution_score(q, s, self.matrix)
                    midline.append('+' if score > 0 else ' ')

        return ''.join(midline)


class BlastParser:
    headers = {
        'tblastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'positive', 'mismatch', 'gap_open', 'query_start', 'query_end', 'seq_start',
                    'seq_end', 'evalue', 'bit_score', 'query_frame', 'seq_frame', 'id'],
        'blastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_open', 'query_start', 'query_end', 'seq_start',
                   'seq_end', 'evalue', 'bit_score', 'query_frame', 'seq_frame', 'id'],
        'blastp': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'positive', 'mismatch', 'gap_open', 'query_start', 'query_end', 'seq_start',
                   'seq_end', 'evalue', 'bit_score', 'query_frame', 'seq_frame', 'id'],
        'blastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_open', 'query_start', 'query_end', 'seq_start',
                   'seq_end', 'evalue', 'bit_score', 'query_frame', 'seq_frame', 'id'],
        'tblastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'mismatch', 'gap_open', 'query_start', 'query_end', 'seq_start',
                    'seq_end', 'evalue', 'bit_score', 'query_frame', 'seq_frame', 'id']
    }
    pd_columns_dtypes = {
        'query_title': 'string',
        'strain': 'string',
        'identity': 'UInt64',
        'perc_identity': 'Float64',
        'query_len': 'UInt64',
        'align_len': 'UInt64',
        'perc_alignment': 'Float64',
        'gaps': 'UInt64',
        'gap_open': 'UInt64',
        'mismatch': 'UInt64',
        'positive': 'UInt64',
        'perc_positive': 'Float64',
        'query_start': 'UInt64',
        'query_end': 'UInt64',
        'seq_start': 'UInt64',
        'seq_end': 'UInt64',
        'query_frame': 'Int8',
        'seq_frame': 'Int8',
        'score': 'UInt64',
        'evalue': 'Float64',
        'bit_score': 'Float64',
        'qseq': 'string',
        'sseq': 'string',
        'qseqid': 'string',
    }

    def __init__(self, file: Path | str, params: dict = None):
        """
        :param file: Path to the blast result file
        :param params: Additional parameters of the blast analysis. They overwrite the ones found in the file, if any
        """

        if file.suffix != '.tsv':
            raise ValueError(f'File {file} is not a tsv file')

        self.file: Path = Path(file)
        self.query_file: Path = Path(str(self.file.with_suffix('.fasta')).replace('_result', '_query'))

        self.metadata = self._read_metadata()
        self.program: str = self.metadata['program']
        self.queries: list = self.metadata['queries']

        if params:
            self.metadata['params'].update(params)

        self.whole_df: pd.DataFrame = self._parse_csv()

    @property
    def df(self):
        """
        Subset of the main columns depending on the blast program used (blastn, blastp, etc.)
        """

        if self.whole_df.empty:
            df = pd.DataFrame()
            return df

        columns = self.headers[self.program]
        return self.whole_df[columns]

    def _read_metadata(self):
        blast_program = 'blastn'
        version = ''
        database = ''
        # fields = ''
        params = dict()
        query_titles = list()
        hits_found = list()

        with open(self.file, 'r', encoding='utf8') as f:
            params_block = False
            for line in f:
                if line.startswith('# [PARAMS]'):
                    params_block = True

                elif line.startswith('# [END PARAMS]'):
                    params_block = False

                elif params_block and line[0] == '#':
                    re_matches = re.search(r'# ([\w ]+): *([\w ]*)$', line)
                    if re_matches is None:
                        continue
                    key = re_matches.group(1).strip()
                    value = re_matches.group(2).strip()
                    params[key] = value

                elif line[0] == '#':

                    match regex_spm.search_in(line):
                        case r'# (BLASTN|BLASTP|BLASTX|TBLASTN|TBLASTX) (\d+\.\d+.\d+\++)\n' as m:
                            blast_program, version = m[1], m[2]
                        case r'# Query: (.+)\n' as m:
                            query_titles.append(m[1])
                        case r'# Database: (.+)\n' as m:
                            database = m[1]
                        # case r'# Fields: (.+)\n' as m:
                        #     fields = m[1]
                        case r'# (\d+) hits found' as m:
                            hits_found.append(m[1])
                        case '_':
                            pass

        if len(query_titles) != len(hits_found):
            raise ValueError(f'Number of queries ({len(query_titles)}) is not equal to number of '
                             f'hits found ({len(hits_found)}). There is an error in the blast result file.')

        queries = list()
        for query_title, hits in zip(query_titles, hits_found):
            queries.append({
                'query_title': query_title,
                'hits': int(hits),
            })

        metadata = {
            'program': blast_program.lower(),
            'version': version,
            'blast_db': database,
            'queries': queries,
            'params': params,
        }

        return metadata

    def _parse_csv(self) -> (pd.DataFrame, dict):
        whole_df = pd.read_csv(self.file, sep='\t', engine='c', encoding='utf8',
                               skip_blank_lines=True, comment='#', header=None, index_col=False,
                               dtype=self.pd_columns_dtypes, names=list(self.pd_columns_dtypes.keys()))

        if whole_df.empty:
            raise EmptyCSVError(f"CSV is empty: {self.file}")

        if '_NODE_' in whole_df['strain'].iloc[0]:
            strain_node_df = whole_df['strain'].str.split('_NODE_', expand=True, regex=False)
            whole_df['strain'] = strain_node_df[0]
            whole_df.insert(2, 'node', value=strain_node_df[1])

        query_titles = list()
        for query in self.metadata['queries']:
            query_titles.extend(repeat(query['query_title'], times=query['hits']))

        if query_titles:
            whole_df['query_title'] = pd.Series(query_titles)

        whole_df['id'] = whole_df.index.copy()

        return whole_df

    def filtered_df(self, identity: float = 60.0, query_cov: float = 50.0) -> pd.DataFrame:
        """
        Filter the results of the blast analysis to keep only the results with query_cov and identity higher
        than the ones provided.
        """

        if not 0 <= identity <= 100:
            raise ValueError(f"Identity must be between 0 and 100, not {identity}")

        if not 0 <= query_cov <= 100:
            raise ValueError(f"Query coverage must be between 0 and 100, not {query_cov}")

        return self.df[(self.df['perc_identity'] >= identity) & (self.df['perc_alignment'] >= query_cov)]

    def alignments(self, indexes=None) -> pd.Series:
        if indexes is None:
            indexes = self.df['id']

        substitution_matrix = self.metadata['params'].get('matrix', 'BLOSUM62')

        return self.whole_df.iloc[indexes].apply(
            lambda row: Alignment(row, self.program, substitution_matrix).align(), axis=1)

    def plot_alignments_bokeh(self, indexes: pd.Series, height: int = None, max_hits: int = 100, sort_by=None):
        """
        Plot the alignments against the query

        :param indexes:
        :param height:
        :param max_hits:
        :param sort_by:
        :return:
        """

        if sort_by is None:
            sort_by = ['evalue', 'perc_alignment']

        def get_color(perc_identity: int):

            if perc_identity >= 90:
                return '#054A29'  # dark green
            elif 70 <= perc_identity < 90:
                return '#37A5BE'  # cyan
            elif 50 <= perc_identity < 70:
                return '#F8E430'  # yellow
            elif 20 <= perc_identity < 50:
                return '#FF9C1A'  # orange
            else:
                return '#FF0A0A'  # red

        # All the indexes have the same query, thus the same query len
        query_title = self.df.iloc[indexes.iloc[0]].query_title
        query_len = self.df.iloc[indexes.iloc[0]].query_len

        rows = self.df.iloc[indexes]
        rows.insert(0, 'table_row_index', value=indexes.index + 1)

        rows = rows.sort_values(by=sort_by, inplace=False, ascending=[True, False])
        rows = rows[:max_hits]

        # Width is +1 because query_start starts at 1 and not at 0, meaning that if query_start = 1 and query_end = 10
        # the width is 10 but 10-1=9. For the same reason the x coordinate is +0.5 because it is the middle of the
        # rectangle
        rows.insert(0, 'graph_index', value=range(1, len(rows) + 1))
        rows.set_index('graph_index', inplace=True)
        rows.insert(0, 'x', value=(rows['query_end'] - rows['query_start']) / 2 + 0.5)
        rows.insert(0, 'y', value=len(rows) - rows.index)
        rows.insert(0, 'width', value=rows['query_end'] - rows['query_start'] + 1)
        rows.insert(0, 'height', value=0.8)
        rows.insert(0, 'fill_color', value=rows['perc_identity'].apply(get_color))

        hits_source = ColumnDataSource(rows)

        query_source = pd.DataFrame.from_dict({
            "x": [query_len / 2],
            "y": [len(rows) + 1],
            "width": [query_len],
            "height": [1],
            "fill_color": ["#58C7C7"],
            "query_title": [query_title],
            "query_len": [query_len],
        })

        height = len(rows) * 10 if height is None else height
        height = max(height, 200)

        plot = figure(width=1000,
                      height=height,
                      y_range=(-2, len(rows) + 5),
                      x_range=(-(query_len / 100 * 5), query_len + (query_len / 100 * 10)),
                      tools="save")

        # Query rectangle
        query_rect = plot.rect(x="x",
                               y="y",
                               width="width",
                               height='height',
                               fill_color="fill_color",
                               fill_alpha=1,
                               source=query_source,
                               line_color=None)

        # Subject rectangles
        sbjct_rect = plot.rect(x="x",
                               y="y",
                               width="width",
                               height='height',
                               fill_color="fill_color",
                               fill_alpha=1,
                               source=hits_source,
                               line_color=None)

        tooltips_query = [
            ("Query", "@query_title"),
            ("Query length", "@query_len"),
        ]

        tooltips_rect = [
            ("Strain", "@strain"),
            ("Node", "@node"),
            ("Row Index", "@table_row_index"),
            ("Perc Identity", "@perc_identity"),
            ("Perc Alignment", "@perc_alignment"),
            ("Evalue", "@evalue"),
            # ("Seq start", "@seq_start"),
            # ("Seq end", "@seq_end"),
            # ("Query start", "@query_start"),
            # ("Query end", "@query_end"),
            ("Alignment length", "@width"),
        ]

        query_hover_tool = HoverTool(renderers=[query_rect], tooltips=tooltips_query, point_policy="follow_mouse")
        plot.add_tools(query_hover_tool)

        sbjct_hover_tool = HoverTool(renderers=[sbjct_rect], tooltips=tooltips_rect, point_policy="follow_mouse")
        plot.add_tools(sbjct_hover_tool)

        # Legend
        rect0 = Rect(width=1, height=1, line_color="#054A29", fill_color="#054A29", fill_alpha=1)
        rect1 = Rect(width=1, height=1, line_color="#37A5BE", fill_color="#37A5BE", fill_alpha=1)
        rect2 = Rect(width=1, height=1, line_color="#F8E430", fill_color="#F8E430", fill_alpha=1)
        rect3 = Rect(width=1, height=1, line_color="#FF9C1A", fill_color="#FF9C1A", fill_alpha=1)
        rect4 = Rect(width=1, height=1, line_color="#FF0A0A", fill_color="#FF0A0A", fill_alpha=1)
        rect0 = plot.add_glyph(rect0)
        rect1 = plot.add_glyph(rect1)
        rect2 = plot.add_glyph(rect2)
        rect3 = plot.add_glyph(rect3)
        rect4 = plot.add_glyph(rect4)

        legend = Legend(title="Percentage identity",
                        orientation='vertical',
                        location='top',
                        items=[
                            (">= 90", [rect0]),
                            ("70 - 90", [rect1]),
                            ("50 - 70", [rect2]),
                            ("20 - 50", [rect3]),
                            ("<20", [rect4]),
                        ])

        plot.add_layout(legend, "right")
        # display legend in top left corner (default is top right corner)

        plot.xaxis.axis_label = "Nucleotides/Aminoacids"
        plot.xgrid.grid_line_color = "black"

        # change things only on the y-grid
        plot.xgrid.grid_line_alpha = 0.2
        plot.xgrid.grid_line_dash = [6, 4]

        plot.yaxis.visible = False
        plot.ygrid.visible = False
        plot.outline_line_color = None
        plot.toolbar.logo = None

        return plot


@st.cache_data(show_spinner=False)
def load_analysis(file: Path, params=None) -> BlastParser:
    """
    Load analysis from file
    :param file: analysis to load
    :param params: Additional parameters of the blast analysis. They overwrite the ones found in the file, if any
    :return:
    """

    return BlastParser(file=file, params=params)
