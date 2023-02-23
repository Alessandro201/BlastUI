import json
from pathlib import Path
from typing import Union

import pandas as pd

HEADER_TBLASTN = ['Query', 'Strain', 'Node', '%_identity', '%_alignment', 'Query_length', 'Alignment_length',
                  'Identities', 'Positives', 'Mismatches', 'Gap_opens', 'Query_start', 'Query_end', 'Match_start',
                  'Match_end', 'Evalue', 'Bit_score', 'Sequence']

HEADER_BLASTN = ['Query', 'Strain', 'Node', '%_identity', '%_alignment', 'Query_length', 'Alignment_length',
                 'Identities', 'Mismatches', 'Gap_opens', 'Query_start', 'Query_end', 'Match_start',
                 'Match_end', 'Evalue', 'Bit_score', 'Sequence']


def json_to_df(json_file: Union[str, Path]) -> pd.DataFrame:
    with open(json_file, 'r') as f:
        data = json.load(f)

    data = get_matches(data)
    df = pd.DataFrame.from_dict(data)
    return data


def read_json(json_file: Union[str, Path]) -> dict:
    with open(json_file, 'r') as f:
        data = json.load(f)

    data = add_extra_info_to_matches(data)
    return data


def add_extra_info_to_matches(data: dict) -> dict:
    index = 0

    for results_of_query in data['BlastOutput2']:
        query_title: str = results_of_query['report']['results']['search']['query_title']
        query_len: int = results_of_query['report']['results']['search']['query_len']

        hits: list = results_of_query['report']['results']['search']['hits']

        for hit in hits:
            for match in hit['hsps']:
                match['index'] = index
                match['query'] = query_title
                match['query_len'] = query_len
                strain_node = hit['description'][0]['accession']
                match['strain'], match['node'] = strain_node.split('_NODE_')

                index += 1

    return data

def get_matches(data: dict) -> dict:
    match_id = 0

    # all_matches = []

    for results_of_query in data['BlastOutput2']:
        program = results_of_query['report']['program']
        query_title: str = results_of_query['report']['results']['search']['query_title']
        query_len: int = results_of_query['report']['results']['search']['query_len']

        hits: list = results_of_query['report']['results']['search']['hits']

        for hit in hits:
            for match in hit['hsps']:
                match['index'] = match_id
                match['query'] = query_title
                match['query_len'] = query_len
                strain_node = hit['description'][0]['accession']
                match['strain'], match['node'] = strain_node.split('_NODE_')

                alignment_length = match['align_len']
                identities = match['identity']
                match['mismatches'] = alignment_length - identities

                gap_opens_query = match['qseq'].count('-')

                match['perc_identity'] = round(identities / alignment_length * 100)
                match['perc_query_cov'] = round((alignment_length + gap_opens_query) / query_len * 100)

                if program == 'blastn':
                    query_strand = match['query_strand']
                    hit_strand = match['query_strand']
                    match['query_orient'] = True if query_strand == 'Plus' else False
                    match['match_orient'] = True if hit_strand == 'Plus' else False

                elif program == 'tblastn':
                    match_frame = match['hit_frame']
                    match['query_orient'] = 'Plus'
                    match['match_orient'] = 'Plus' if match_frame > 0 else 'Minus'

                # all_matches.append(match)
                match_id += 1

    return data


def filter_results(df: pd.DataFrame, identity: float, query_cov: float) -> pd.DataFrame:
    """
    Filter the results of the blast analysis to keep only the results with query_cov and identity higher
    than the ones provided.
    """

    if identity > 100 or identity < 0:
        raise ValueError(f"Identity must be between 0 and 100, not {identity}")

    if query_cov > 100 or query_cov < 0:
        raise ValueError(f"Query coverage must be between 0 and 100, not {query_cov}")

    new_df = df[(df['perc_identity'] >= identity) & (df['perc_alignment'] >= query_cov)]

    return new_df


def generate_tblastn_df(data: dict) -> pd.DataFrame:
    table = []

    for results_of_query in data['BlastOutput2']:

        query = results_of_query['report']['results']['search']['query_title']
        query_len = results_of_query['report']['results']['search']['query_len']
        hits: list = results_of_query['report']['results']['search']['hits']

        for hit in hits:
            strain_node = hit['description'][0]['id']
            strain, node = strain_node.split('_NODE_')

            for match in hit['hsps']:
                index = match['index']

                alignment_length = match['align_len']

                identities = match['identity']
                positives = match['positive']
                mismatches = alignment_length - identities
                gap_opens = match['gaps']

                perc_identity = round(identities / alignment_length * 100)
                perc_query_cov = round(alignment_length / query_len * 100)

                query_start = match['query_from']
                query_end = match['query_to']
                match_start = match['hit_from']
                match_end = match['hit_to']

                bit_score = match['bit_score']
                evalue = match['evalue']

                sequence = match['hseq']

                table.append([query,
                              strain,
                              node,
                              perc_identity,
                              perc_query_cov,
                              query_len,
                              alignment_length,
                              identities,
                              positives,
                              mismatches,
                              gap_opens,
                              query_start,
                              query_end,
                              match_start,
                              match_end,
                              evalue,
                              bit_score,
                              sequence,
                              index, ])

    header = HEADER_TBLASTN + ['index']
    df = pd.DataFrame(table, columns=header)
    return df


def generate_blastn_df(data: dict) -> pd.DataFrame:
    table = []

    for results_of_query in data['BlastOutput2']:

        query = results_of_query['report']['results']['search']['query_title']
        query_len = results_of_query['report']['results']['search']['query_len']
        hits: list = results_of_query['report']['results']['search']['hits']

        for hit in hits:
            strain_node = hit['description'][0]['id']
            strain, node = strain_node.split('_NODE_')

            for match in hit['hsps']:
                index = match['index']

                alignment_length = match['align_len']

                identities = match['identity']
                mismatches = alignment_length - identities
                gap_opens = match['gaps']

                perc_identity = round(identities / alignment_length * 100)
                perc_query_cov = round(alignment_length / query_len * 100)

                query_start = match['query_from']
                query_end = match['query_to']
                match_start = match['hit_from']
                match_end = match['hit_to']

                bit_score = match['bit_score']
                evalue = match['evalue']

                sequence = match['hseq']

                table.append([query,
                              strain,
                              node,
                              perc_identity,
                              perc_query_cov,
                              query_len,
                              alignment_length,
                              identities,
                              mismatches,
                              gap_opens,
                              query_start,
                              query_end,
                              match_start,
                              match_end,
                              evalue,
                              bit_score,
                              sequence,
                              index])

    header = HEADER_BLASTN + ['index']
    df = pd.DataFrame(table, columns=header)
    return df


def write_table(df: pd.DataFrame, output_file: Union[str, Path]):
    df.to_csv(output_file, sep='\t', header=True, index=False)


def generate_tblastn_alignment(match: dict):
    index = match['index']
    strain = match['strain']
    node = match['node']
    query_title = match['query']
    query_len = match['query_len']
    alignment_length = match['align_len']

    query_start = match['query_from']
    query_end = match['query_to']
    match_start = match['hit_from']
    match_end = match['hit_to']

    identities = match['identity']
    positives = match['positive']
    mismatches = alignment_length - identities
    gap_opens = match['gaps']

    bit_score = match['bit_score']
    score = match['score']
    evalue = match['evalue']

    query = match['qseq']
    midline = match['midline']
    sequence = match['hseq']

    gap_opens_query = query.count('-')

    perc_identity = round(identities / alignment_length * 100)
    perc_query_cov = round((alignment_length + gap_opens_query) / query_len * 100)
    perc_gaps = round(gap_opens / alignment_length * 100)
    perc_mismatches = round(mismatches / alignment_length * 100)
    perc_positives = round(positives / alignment_length * 100)

    match_frame = match['hit_frame']

    # Match statistics
    alignment_text = f">{query_title} \n" \
                     f"Strain = {strain}, Node = {node}\n" \
                     f"\tScore = {round(bit_score)} bits ({score}), E-value = {evalue :.1g} \n" \
                     f"\tIdentities = {identities}/{alignment_length} ({perc_identity}%), " \
                     f"Query coverage = {alignment_length + gap_opens_query}/{query_len} ({perc_query_cov}%) \n" \
                     f"\tPositives = {positives}/{alignment_length} ({perc_positives}%), " \
                     f"Mismatches = {mismatches}/{alignment_length} ({perc_mismatches}%), " \
                     f"Gaps = {gap_opens}/{alignment_length} ({perc_gaps}%)\n" \
                     f"\tFrame = {match_frame}\n\n"

    # The query is a protein sequence but the match is a DNA sequence which blast translated
    # Hence the indexes are multiplied by 3 to get the correct position in the DNA sequence
    dna_to_prot = 3
    prot_to_dna = 1

    query_forward_strand = True
    match_forward_strand = True if match_frame > 0 else False

    # The number of padding spaces for the indexes in the alignment depends on the number of digits
    pad = max(len(str(match_start)), len(str(query_start)))

    index = 0
    query_gap = 0
    match_gap = 0
    prev_query_gap = 0
    prev_match_gap = 0
    while index < len(sequence):

        query_gap += query[index:index + 60].count('-')
        match_gap += sequence[index:index + 60].count('-')

        # Query
        if query_forward_strand:
            # Forward strand
            q_start = query_start + (index - prev_query_gap) * prot_to_dna
            q_end = query_start + (60 + index - query_gap) * prot_to_dna - 1
            q_end = min(q_end, query_end)

        else:
            # Reverse strand
            q_start = query_start - (index - prev_query_gap) * prot_to_dna
            q_end = query_start - (60 + index - query_gap) * prot_to_dna + 1
            q_end = max(q_end, query_start)

        # Subject
        if match_forward_strand:
            # Forward strand
            s_start = match_start + (index - prev_match_gap) * dna_to_prot
            s_end = match_start + (60 + index - match_gap) * dna_to_prot - 1
            s_end = min(s_end, match_end)

        else:
            # Reverse strand
            s_start = match_end - (index - prev_match_gap) * dna_to_prot
            s_end = match_end - (60 + index - match_gap) * dna_to_prot + 1
            s_end = max(s_end, match_start)

        prev_query_gap = query_gap
        prev_match_gap = match_gap

        alignment_text += f"Query  {q_start :{pad}}  {query[index:index + 60]}  {q_end :{pad}}\n" \
                          f"       {' ' * pad}  {midline[index:index + 60]} \n" \
                          f"Sbjct  {s_start :{pad}}  {sequence[index:index + 60]}  {s_end :{pad}}\n\n"
        index += 60

    return alignment_text


def generate_blastn_alignment(match: dict):
    index = match['index']
    strain = match['strain']
    node = match['node']
    query_title = match['query']
    query_len = match['query_len']
    alignment_length = match['align_len']

    query_start = match['query_from']
    query_end = match['query_to']
    match_start = match['hit_from']
    match_end = match['hit_to']

    identities = match['identity']
    mismatches = alignment_length - identities
    gap_opens = match['gaps']

    perc_identity = round(identities / alignment_length * 100)
    perc_query_cov = round(alignment_length / query_len * 100)
    perc_gaps = round(gap_opens / alignment_length * 100)
    perc_mismatches = round(mismatches / alignment_length * 100)

    bit_score = match['bit_score']
    score = match['score']
    evalue = match['evalue']

    query = match['qseq']
    midline = match['midline']
    sequence = match['hseq']

    query_strand = match['query_strand']
    hit_strand = match['hit_strand']

    alignment_text = f">{query_title} \n" \
                     f"Strain = {strain}, Node = {node}\n" \
                     f"\tScore = {round(bit_score)} bits ({score}), E-value = {evalue :.1g} \n" \
                     f"\tIdentities = {identities}/{alignment_length} ({perc_identity}%), " \
                     f"Query coverage = {alignment_length}/{query_len} ({perc_query_cov}%) \n" \
                     f"\tMismatches = {mismatches}/{alignment_length} ({perc_mismatches}%), " \
                     f"Gaps = {gap_opens}/{alignment_length} ({perc_gaps}%)\n" \
                     f"\tStrand = {query_strand}/{hit_strand}\n\n"

    # The query is a protein sequence but the match is a DNA sequence which blast translated
    # Hence the indexes are multiplied by 3 to get the correct position in the DNA sequence
    dna_to_prot = 1
    prot_to_dna = 1

    query_forward_strand = True if query_strand == 'Plus' else False
    match_forward_strand = True if hit_strand == 'Plus' else False

    # The number of padding spaces for the indexes in the alignment depends on the number of digits
    pad = max(len(str(match_start)), len(str(query_start)))

    index = 0
    query_gap = 0
    match_gap = 0
    prev_query_gap = 0
    prev_match_gap = 0
    while index < len(sequence):

        query_gap += query[index:index + 60].count('-')
        match_gap += sequence[index:index + 60].count('-')

        # Query
        if query_forward_strand:
            # Forward strand
            q_start = query_start + (index - prev_query_gap) * prot_to_dna
            q_end = query_start + (60 + index - query_gap) * prot_to_dna - 1
            q_end = min(q_end, query_end)

        else:
            # Reverse strand
            q_start = query_start - (index - prev_query_gap) * prot_to_dna
            q_end = query_start - (60 + index - query_gap) * prot_to_dna + 1
            q_end = max(q_end, query_end)

        # Subject
        if match_forward_strand:
            # Forward strand
            s_start = match_start + (index - prev_match_gap) * dna_to_prot
            s_end = match_start + (60 + index - match_gap) * dna_to_prot - 1
            s_end = min(s_end, match_end)

        else:
            # Reverse strand
            s_start = match_start - (index - prev_match_gap) * dna_to_prot
            s_end = match_start - (60 + index - match_gap) * dna_to_prot + 1
            s_end = max(s_end, match_end)

        prev_query_gap = query_gap
        prev_match_gap = match_gap

        alignment_text += f"Query  {q_start :{pad}}  {query[index:index + 60]}  {q_end :{pad}}\n" \
                          f"       {' ' * pad}  {midline[index:index + 60]} \n" \
                          f"Sbjct  {s_start :{pad}}  {sequence[index:index + 60]}  {s_end :{pad}}\n\n"
        index += 60

    return alignment_text


def tblastn_alignments(data: dict):
    text = ''
    for results_of_query in data['BlastOutput2']:

        query_title = results_of_query['report']['results']['search']['query_title']
        hits: list = results_of_query['report']['results']['search']['hits']

        text += f"Query: {query_title}\n\n"

        for hit in hits:
            strain_node = hit['description'][0]['id']

            text += f">{strain_node}\n\n"

            for match in hit['hsps']:
                text += generate_tblastn_alignment(match) + '\n'


def print_blastn_alignments(data: dict):
    text = ''
    for results_of_query in data['BlastOutput2']:

        query_title = results_of_query['report']['results']['search']['query_title']
        hits: list = results_of_query['report']['results']['search']['hits']

        text += f"Query: {query_title}\n\n"

        for hit in hits:
            strain_node = hit['description'][0]['id']

            text += f">{strain_node}\n\n"

            for match in hit['hsps']:
                text += generate_blastn_alignment(match) + '\n'


def get_matches_by_id(data: dict, match_id: str, query: str):
    for results_of_query in data['BlastOutput2']:
        hits: list = results_of_query['report']['results']['search']['hits']
        query_title: str = results_of_query['report']['results']['search']['query_title']

        if query_title != query:
            continue

        for hit in hits:
            strain_node_id = hit['description'][0]['id']
            if strain_node_id == match_id:
                for match in hit['hsps']:
                    yield match


def get_matches_by_index(data, match_indexes):
    # match_index will keep track of the indexes
    match_indexes = tuple(match_indexes)

    # matches will keep the matches found in the order given by match_indexes
    matches = list(match_indexes)
    for results_of_query in data['BlastOutput2']:
        hits: list = results_of_query['report']['results']['search']['hits']

        for hit in hits:
            for match in hit['hsps']:
                if match['index'] in match_indexes:
                    list_i = match_indexes.index(match['index'])
                    matches[list_i] = match

    return matches


if __name__ == '__main__':

    # data = read_json('../Analysis/tblastn_results.json')
    # data = add_query_length_to_matches(data)
    # data = filter_results(data, identity=90, query_cov=60)
    # table = generate_tblastn_tsv(data)
    # write_table(table, HEADER_TBLASTN, '../Analysis/tblastn_results.tsv')
    # print_tblastn_alignments(data)

    # query_title = 'treP PTS system trehalose-specific EIIBC component'
    # subject_id = 'Sta_69_NODE_30'
    #

    data = read_json('../Analysis/blastn_results.json')
    data = add_query_length_to_matches(data)
    # # data = filter_results(data, identity=90, query_cov=60)
    # table = generate_blastn_tsv(data)
    # write_table(table, HEADER_BLASTN, '../Analysis/blastn_results.tsv')
    # print_blastn_alignments(data)

    query_title = 'TestQuery'
    subject_id = 'Sta_89_NODE_9'

    matches = get_matches_by_id(data, match_id=subject_id, query=query_title)
    print(f"Query: {query_title}")
    print(f"> {subject_id}\n")
    for match in matches:
        print(generate_blastn_alignment(match))
