import json as json
import pandas as pd
import numpy as np
import streamlit as st
import re
from pathlib import Path
from typing import Iterable

from bokeh.models import (ColumnDataSource, HoverTool)
from bokeh.models import Legend, Rect
from bokeh.plotting import figure


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

        return self.whole_df.iloc[indexes].apply(self._alignment, axis=1)

    def _alignment(self, row) -> str:
        query_title: str = row['query_title']
        query_len: int = row['query_len']
        strain: str = row['strain']
        node: int = int(row['node'])

        qseq: str = row['qseq']
        midline: str = row['midline']
        sseq: str = row['hseq']

        align_len: int = row['align_len']
        identity: int = row['identity']

        gaps: int = row['gaps']
        q_gap: int = qseq.count('-')
        s_gap: int = gaps - q_gap
        regex_matches = re.findall('-*', qseq + '|' + sseq)
        gap_opens: int = len([match for match in regex_matches if match])
        mismatch: int = align_len - identity - gaps

        perc_identity: float = row['perc_identity']
        perc_alignment: float = row['perc_alignment']
        perc_gaps: float = row['perc_gaps']
        perc_mismatches = row['perc_mismatches']

        q_start: int = row['query_from']
        q_end: int = row['query_to']
        s_start: int = row['hit_from']
        s_end: int = row['hit_to']

        evalue: float = row['evalue']
        bit_score: float = row['bit_score']
        score: int = row['score']

        if self.program in ('tblastn', 'blastx', 'blastp', 'tblastx'):
            positive: int = row['positive']
            perc_positives: float = round(positive / align_len * 100)

        match self.program:
            case 'blastn':
                q_orient: str = row['query_strand']
                s_orient: str = row['hit_strand']
            case 'tblastn' | 'tblastx':
                s_frame = row['hit_frame']
            case 'blastx':
                q_frame = row['query_frame']
            case 'blastp' | _:
                pass

        match self.program:
            case 'tblastn' | 'tblastx':
                alignment_text = f">{query_title} \n" \
                                 f"Strain = {strain}, Node = {node}\n" \
                                 f"\tScore = {round(bit_score)} bits ({score}), " \
                                 f"E-value = {evalue :.3g} \n" \
                                 f"\tIdentities = {identity}/{align_len} ({perc_identity}%), " \
                                 f"Query coverage = {align_len}/{query_len} ({perc_alignment}%), " \
                                 f"Gap opens = {gap_opens}\n" \
                                 f"\tPositives = {positive}/{align_len} ({perc_positives}%), " \
                                 f"Mismatches = {mismatch}/{align_len} ({perc_mismatches}%), " \
                                 f"Gaps = {gaps}/{align_len} ({perc_gaps}%)\n" \
                                 f"\tFrame = {s_frame}\n\n"

            case 'blastn':
                alignment_text = f">{query_title} \n" \
                                 f"Strain = {strain}, Node = {node}\n" \
                                 f"\tScore = {round(bit_score)} bits ({score}), " \
                                 f"E-value = {evalue :.3g} \n" \
                                 f"\tIdentities = {identity}/{align_len} ({perc_identity}%), " \
                                 f"Query coverage = {align_len}/{query_len} ({perc_alignment}%), " \
                                 f"Gap opens = {gap_opens}\n" \
                                 f"Mismatches = {mismatch}/{align_len} ({perc_mismatches}%), " \
                                 f"Gaps = {gaps}/{align_len} ({perc_gaps}%)\n" \
                                 f"\tStrand = {q_orient}/{s_orient}\n\n"

            case 'blastp':
                alignment_text = f">{query_title} \n" \
                                 f"Strain = {strain}, Node = {node}\n" \
                                 f"\tScore = {round(bit_score)} bits ({score}), " \
                                 f"E-value = {evalue :.3g} \n" \
                                 f"\tIdentities = {identity}/{align_len} ({perc_identity}%), " \
                                 f"Query coverage = {align_len}/{query_len} ({perc_alignment}%), " \
                                 f"Gap opens = {gap_opens}\n" \
                                 f"\tPositives = {positive}/{align_len} ({perc_positives}%), " \
                                 f"Mismatches = {mismatch}/{align_len} ({perc_mismatches}%), " \
                                 f"Gaps = {gaps}/{align_len} ({perc_gaps}%)\n\n"

            case 'blastx':
                alignment_text = f">{query_title} \n" \
                                 f"Strain = {strain}, Node = {node}\n" \
                                 f"\tScore = {round(bit_score)} bits ({score}), " \
                                 f"E-value = {evalue :.3g} \n" \
                                 f"\tIdentities = {identity}/{align_len} ({perc_identity}%), " \
                                 f"Query coverage = {align_len}/{query_len} ({perc_alignment}%), " \
                                 f"Gap opens = {gap_opens}\n" \
                                 f"\tPositives = {positive}/{align_len} ({perc_positives}%), " \
                                 f"Mismatches = {mismatch}/{align_len} ({perc_mismatches}%), " \
                                 f"Gaps = {gaps}/{align_len} ({perc_gaps}%)\n" \
                                 f"\tQuery frame = {q_frame}\n\n"

            case _:
                alignment_text = f">{query_title} \n"

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
        pad = max(len(str(s_start)),
                  len(str(q_start)),
                  len(str(s_end)),
                  len(str(q_end)))

        index = 0
        query_gap = 0
        seq_gap = 0
        prev_query_gap = 0
        prev_seq_gap = 0

        while index < len(sseq):
            query_gap += qseq[index:index + 60].count('-')
            seq_gap += sseq[index:index + 60].count('-')

            # Query
            if q_start < q_end:
                # Forward strand
                q_start = q_start + (index - prev_query_gap) * q_multiplier
                q_end = q_start + (index + 60 - query_gap) * q_multiplier - 1
                q_end = min(q_end, q_end)

            else:
                # Reverse strand
                q_start = q_end - (index - prev_query_gap) * q_multiplier
                q_end = q_end - (index + 60 - query_gap) * q_multiplier + 1
                q_end = max(q_end, q_start)

            # Subject
            if s_start < s_end:
                # Forward strand
                s_start = s_start + (index - prev_seq_gap) * s_multiplier
                s_end = s_start + (index + 60 - seq_gap) * s_multiplier - 1
                s_end = min(s_end, s_end)

            else:
                # Reverse strand
                s_start = s_start - (index - prev_seq_gap) * s_multiplier
                s_end = s_start - (index + 60 - seq_gap) * s_multiplier + 1
                s_end = max(s_end, s_end)

            prev_query_gap = query_gap
            prev_seq_gap = seq_gap

            alignment_text += f"Query  {q_start :{pad}}  {qseq[index:index + 60]}  {q_end :{pad}}\n" \
                              f"       {' ' * pad}  {midline[index:index + 60]} \n" \
                              f"Sbjct  {s_start :{pad}}  {sseq[index:index + 60]}  {s_end :{pad}}\n\n"
            index += 60

        return alignment_text

    def plot_alignments_bokeh(self, indexes: Iterable[int], height: int = None, max_hits: int = 100, sort_by=None):
        """
        Plot the alignments against the query

        :param indexes:
        :param height:
        :param max_hits:
        :param sort_by:
        :return:
        """

        # I don't care about the index, I just want the values
        if isinstance(indexes, pd.Series):
            indexes.reset_index(drop=True, inplace=True)

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
        query_title = self.df.iloc[indexes[0]].query_title
        for query in self.metadata['queries']:
            if query['query_title'] == query_title:
                query_len = query['query_len']
                break

        rows = self.df.iloc[indexes]
        rows = rows.sort_values(by=sort_by, inplace=False, ascending=[True, False])
        rows = rows[:max_hits]

        rows.insert(0, 'index', value=range(1, len(rows) + 1))
        rows.set_index('index', inplace=True)
        rows.insert(0, 'x', value=(rows['query_to'] + rows['query_from']) // 2)
        rows.insert(0, 'y', value=-rows.index)
        rows.insert(0, 'width', value=(rows['query_to'] - rows['query_from']))
        rows.insert(0, 'height', value=0.8)
        rows.insert(0, 'color', value=rows['perc_identity'].apply(get_color))

        source = ColumnDataSource(rows)

        height = 100 + len(rows) * 10 if height is None else height

        plot = figure(width=1000, height=height, tools="save")

        # Query rectangle
        plot.rect(x=query_len // 2,
                  y=1,
                  width=query_len,
                  height=1,
                  fill_color="#58c7c7",
                  line_color=None)

        # Subject rectangles
        sbjct_rect = plot.rect(x="x",
                               y="y",
                               width="width",
                               height='height',
                               fill_color="color",
                               fill_alpha=1,
                               source=source,
                               line_color=None)

        tooltips = [
            ("Strain", "@strain"),
            ("Node", "@node"),
            ("Perc Identity", "@perc_identity"),
            ("Perc Alignment", "@perc_alignment"),
            ("Evalue", "@evalue"),
            ("Seq start", "@hit_from"),
            ("Seq end", "@hit_to"),
            ("Query start", "@query_from"),
            ("Query end", "@query_to"),
        ]

        sbjct_hover_tool = HoverTool(renderers=[sbjct_rect], tooltips=tooltips, point_policy="follow_mouse")
        plot.add_tools(sbjct_hover_tool)

        rect0 = Rect(x="x", y="color", width=1, height=1, line_color="#054A29", fill_color="#054A29", fill_alpha=1)
        rect1 = Rect(x="x", y="color", width=1, height=1, line_color="#37A5BE", fill_color="#37A5BE", fill_alpha=1)
        rect2 = Rect(x="x", y="color", width=1, height=1, line_color="#F8E430", fill_color="#F8E430", fill_alpha=1)
        rect3 = Rect(x="x", y="color", width=1, height=1, line_color="#FF9C1A", fill_color="#FF9C1A", fill_alpha=1)
        rect4 = Rect(x="x", y="color", width=1, height=1, line_color="#FF0A0A", fill_color="#FF0A0A", fill_alpha=1)
        rect0 = plot.add_glyph(source, rect0)
        rect1 = plot.add_glyph(source, rect1)
        rect2 = plot.add_glyph(source, rect2)
        rect3 = plot.add_glyph(source, rect3)
        rect4 = plot.add_glyph(source, rect4)

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
