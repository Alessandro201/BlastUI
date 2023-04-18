import json as json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import streamlit as st
from bokeh.models import (ColumnDataSource, HoverTool)
from bokeh.models import Legend, Rect
from bokeh.plotting import figure


class Alignment:
    def __init__(self, row, program):
        self.program = program

        if self.program not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
            raise ValueError(f"Invalid program: {self.program}")

        self.query_title: str = row['query_title']
        self.query_len: int = row['query_len']
        self.strain: str = row['strain']
        self.node: int = int(row['node'])

        self.qseq: str = row['qseq']
        self.midline: str = row['midline']
        self.sseq: str = row['hseq']

        self.align_len: int = row['align_len']
        self.identity: int = row['identity']

        self.gaps: int = row['gaps']
        self.q_gap: int = row['q_gaps']
        self.s_gap: int = row['s_gaps']
        self.gap_opens: int = row['gap_opens']
        self.mismatch: int = row['mismatch']

        self.perc_identity: float = row['perc_identity']
        self.perc_alignment: float = row['perc_alignment']
        self.perc_gaps: float = row['perc_gaps']
        self.perc_mismatches = row['perc_mismatches']

        self.q_from: int = row['query_from']
        self.q_to: int = row['query_to']
        self.s_from: int = row['hit_from']
        self.s_to: int = row['hit_to']

        self.evalue: float = row['evalue']
        self.bit_score: float = row['bit_score']
        self.score: int = row['score']

        if self.program in ('tblastn', 'blastx', 'blastp', 'tblastx'):
            self.positive: int = row['positive']
            self.perc_positives: float = round(self.positive / self.align_len * 100)

        match self.program:
            case 'blastn':
                self.q_orient: str = row['query_strand']
                self.s_orient: str = row['hit_strand']
            case 'blastp':
                self.q_orient: str = "forward"
                self.s_orient: str = "forward"
            case 'blastx':
                self.q_frame = row['query_frame']
                self.q_orient: str = "forward" if int(self.q_frame) > 0 else "reverse"
                self.s_orient: str = "forward"
            case 'tblastn':
                self.s_frame = row['hit_frame']
                self.q_orient: str = "forward"
                self.s_orient: str = "forward" if int(self.s_frame) > 0 else "reverse"
            case 'tblastx':
                self.q_frame = row['query_frame']
                self.s_frame = row['hit_frame']
                self.q_orient: str = "forward" if int(self.q_frame) > 0 else "reverse"
                self.s_orient: str = "forward" if int(self.s_frame) > 0 else "reverse"

        # If the query or the sequence are translated, their positions in the alignment are multiplied by 3
        match self.program:
            case 'tblastn':
                self.q_multiplier = 1
                self.s_multiplier = 3

            case 'blastx':
                self.q_multiplier = 3
                self.s_multiplier = 1

            case 'tblastx':
                self.s_multiplier = 3
                self.q_multiplier = 3

            case 'blastp' | 'blastn' | _:
                self.s_multiplier = 1
                self.q_multiplier = 1

        # The number of padding spaces for the indexes in the alignment depends on the number of digits
        self.pad = max(len(str(self.s_from)),
                       len(str(self.q_from)),
                       len(str(self.s_to)),
                       len(str(self.q_to)))

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
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tFrame = {self.s_frame}\n\n"

        return alignment_text

    def _blastn_header(self):
        alignment_text = f">{self.query_title} \n" \
                         f"Strain = {self.strain}, Node = {self.node}\n" \
                         f"\tScore = {round(self.bit_score)} bits ({self.score}), " \
                         f"E-value = {self.evalue :.3g} \n" \
                         f"\tIdentities = {self.identity}/{self.align_len} ({self.perc_identity}%), " \
                         f"Query coverage = {self.align_len}/{self.query_len} ({self.perc_alignment}%), " \
                         f"Gap opens = {self.gap_opens}\n" \
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
                         f"Gaps = {self.gaps}/{self.align_len} ({self.perc_gaps}%)\n" \
                         f"\tStrand = {self.q_orient}/{self.s_orient}\n\n"

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
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
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
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
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
                         f"Mismatches = {self.mismatch}/{self.align_len} ({self.perc_mismatches}%), " \
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
                q_start = self.q_from + (index - prev_query_gaps) * self.q_multiplier
                q_end = self.q_from + (index + 60 - query_gaps) * self.q_multiplier - 1
                q_end = min(q_end, self.q_to)

            else:
                # Reverse strand or negative frame
                q_start = self.q_to - (index - prev_query_gaps) * self.q_multiplier
                q_end = self.q_to - (index + 60 - query_gaps) * self.q_multiplier + 1
                q_end = max(q_end, self.q_from)

            # Subject
            if self.s_orient == 'forward':
                # Forward strand or positive frame
                s_start = self.s_from + (index - prev_seq_gaps) * self.s_multiplier
                s_end = self.s_from + (index + 60 - seq_gaps) * self.s_multiplier - 1
                s_end = min(s_end, self.s_to)

            else:
                # Reverse strand or negative frame
                s_start = self.s_to - (index - prev_seq_gaps) * self.s_multiplier
                s_end = self.s_to - (index + 60 - seq_gaps) * self.s_multiplier + 1
                s_end = max(s_end, self.s_from)

            prev_query_gaps = query_gaps
            prev_seq_gaps = seq_gaps

            alignment_text += f"Query  {q_start :{self.pad}}  {self.qseq[index:index + 60]}  {q_end :{self.pad}}\n" \
                              f"       {' ' * self.pad}  {self.midline[index:index + 60]} \n" \
                              f"Sbjct  {s_start :{self.pad}}  {self.sseq[index:index + 60]}  {s_end :{self.pad}}\n\n"
            index += 60

        return alignment_text


class BlastResponse:
    headers = {
        'tblastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'positive', 'mismatch', 'gap_opens', 'query_from', 'query_to', 'hit_from',
                    'hit_to', 'evalue', 'bit_score', 'query_id', 'id'],
        'blastn': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_opens', 'query_from', 'query_to', 'hit_from',
                   'hit_to', 'evalue', 'bit_score', 'query_id', 'id'],
        'blastp': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'positive', 'mismatch', 'gap_opens', 'query_from', 'query_to', 'hit_from',
                   'hit_to', 'evalue', 'bit_score', 'query_id', 'id'],
        'blastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                   'identity', 'mismatch', 'gap_opens', 'query_from', 'query_to', 'hit_from',
                   'hit_to', 'evalue', 'bit_score', 'query_id', 'id'],
        'tblastx': ['query_title', 'strain', 'node', 'perc_identity', 'perc_alignment', 'query_len', 'align_len',
                    'identity', 'mismatch', 'gap_opens', 'query_from', 'query_to', 'hit_from',
                    'hit_to', 'evalue', 'bit_score', 'query_id', 'id']
    }

    def __init__(self, json_file: Path):
        self.json_file = json_file
        self.whole_df: pd.DataFrame = pd.DataFrame()
        self.metadata: dict = {}
        self.messages: list[str] = []
        self.queries: list[dict] = []
        self.program: str = ''

        json_data = self._read_json(json_file)
        self.whole_df, metadata = self._parse_json(json_data)
        del json_data

        self.metadata = metadata
        self.program: str = metadata['program']
        self.queries: list[dict] = metadata['queries']
        self.messages: list[str] = metadata['blast_messages']

        # Keep only the columns we need depending on the program (blastn, blastp, etc.)
        columns = self.headers[self.program]
        if self.whole_df.empty:
            self.df: pd.DataFrame = pd.DataFrame()
            return

        self.df: pd.DataFrame = self.whole_df[columns]

    @staticmethod
    def _read_json(json_file: Path):
        if isinstance(json_file, Path):
            with open(json_file, 'r') as f:
                return json.load(f)
        else:
            raise TypeError(f"Expected a dict or Path object, got {type(json_file)} instead")

    @staticmethod
    def _parse_json(json: dict) -> (pd.DataFrame, dict):
        index = 0
        metadata = {
            'queries': list(),
            'blast_messages': list(),
            'params': json['BlastOutput2'][0]['report']['params'],
            'program': json['BlastOutput2'][0]['report']['program'],
            'version': json['BlastOutput2'][0]['report']['version'],
            'reference': json['BlastOutput2'][0]['report']['reference'],
            'search_target': json['BlastOutput2'][0]['report']['search_target'],
        }

        whole_df = pd.DataFrame()

        for report in json['BlastOutput2']:
            report = report['report']

            # If you don't insert a header in the query sequence, there won't be a query_title
            query_id: str = report['results']['search']['query_id']
            query_len: int = report['results']['search']['query_len']
            try:
                query_title: str = report['results']['search']['query_title']
            except KeyError:
                query_title = query_id

            metadata['queries'].append({'query_id': query_id,
                                        'query_title': query_title,
                                        'query_len': query_len})

            if "message" in report['results']['search']:
                metadata['blast_messages'].append(f"Query ***\"{query_title}\"***: \n"
                                                  f"{report['results']['search']['message']}")
                continue

            hits = list()
            for strain_hits in report['results']['search']['hits']:
                strain_node = strain_hits['description'][0]['accession']
                strain, node = strain_node.split('_NODE_')
                for hit in strain_hits['hsps']:
                    hit['id'] = index
                    hit['strain'] = strain
                    hit['node'] = node
                    index += 1

                hits.extend(strain_hits['hsps'])

            df = pd.DataFrame().from_records(hits)
            df['query_title'] = query_title
            df['query_id'] = query_id
            df['query_len']: int = query_len

            whole_df = pd.concat((whole_df, df), ignore_index=True)

            whole_df['q_gaps'] = whole_df['qseq'].str.count('-')
            whole_df['s_gaps'] = whole_df['gaps'] - whole_df['q_gaps']
            whole_df['mismatch'] = whole_df['align_len'] - whole_df['identity'] - whole_df['gaps']

            count_gap_opens = lambda x: len([match for match in re.findall('-+', x) if match])
            qseq_gapopens = whole_df['qseq'].apply(count_gap_opens)
            sseq_gapopens = whole_df['hseq'].apply(count_gap_opens)
            whole_df['gap_opens'] = qseq_gapopens + sseq_gapopens

            whole_df['perc_identity'] = np.round(whole_df['identity'] / whole_df['align_len'] * 100, 1)
            whole_df['perc_alignment'] = np.round((whole_df['align_len'] - whole_df['q_gaps']) /
                                                  whole_df['query_len'] * 100, 1)
            whole_df['perc_gaps'] = np.round(whole_df['gaps'] / whole_df['align_len'] * 100)
            whole_df['perc_mismatches'] = np.round(whole_df['mismatch'] /
                                                   (whole_df['align_len'] - whole_df['gaps']) * 100)

        if not whole_df.empty:
            whole_df['index'] = whole_df['id'].copy()
            whole_df.set_index('index', inplace=True)

        return whole_df, metadata

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

        return self.whole_df.iloc[indexes].apply(lambda row: Alignment(row, self.program).align(), axis=1)

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
        for query in self.metadata['queries']:
            if query['query_title'] == query_title:
                query_len = query['query_len']
                break

        rows = self.df.iloc[indexes]
        rows.insert(0, 'table_row_index', value=indexes.index + 1)

        rows = rows.sort_values(by=sort_by, inplace=False, ascending=[True, False])
        rows = rows[:max_hits]

        # Width is +1 because query_from starts at 1 and not at 0, meaning that if query_from = 1 and query_to = 10
        # the width is 10 but 10-1=9. For the same reason the x coordinate is +0.5 because it is the middle of the
        # rectangle
        rows.insert(0, 'graph_index', value=range(1, len(rows) + 1))
        rows.set_index('graph_index', inplace=True)
        rows.insert(0, 'x', value=(rows['query_to'] - rows['query_from']) / 2 + 0.5)
        rows.insert(0, 'y', value=len(rows) - rows.index)
        rows.insert(0, 'width', value=rows['query_to'] - rows['query_from'] + 1)
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
            ("Seq start", "@hit_from"),
            ("Seq end", "@hit_to"),
            ("Query start", "@query_from"),
            ("Query end", "@query_to"),
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
def load_analysis(json_path: Path) -> BlastResponse:
    """
    Placing the loading in a separate function allows streamlit (st) to cache the data
    :param json_path:
    :return:
    """

    return BlastResponse(json_file=json_path)
