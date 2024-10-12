import aux_mod.dna_rna_aux_mod as dnarna
import aux_mod.filter_fastq_aux_mod as filter
import os
from typing import Union, TextIO

instruments = {
    "transcribe": dnarna.transcribe,
    "reverse": dnarna.reverse,
    "complement": dnarna.complement,
    "reverse_complement": dnarna.reverse_complement,
}


def run_dna_rna_tools(*args: str) -> Union[str, list]:
    """Function run_dna_rna_tools

    Take DNA sequence and return transcribed, reverse,
    complement or reverse complement sequence

    *Args: str

    Returns: str | list
    """
    result = []
    data = args[:-1]
    tool = instruments[args[-1]]
    for sequence in data:
        if dnarna.is_dna(sequence
                         ) is False and dnarna.is_rna(sequence) is False:
            result.append("This is not DNA or RNA sequence!")
        else:
            result.append(tool(sequence))
    if len(result) == 1:
        return result[0]
    else:
        return result


def filter_fastq(
    input_fastq: TextIO,
    output_fastq: Union[str, os.PathLike],
    gc_bounds: Union[tuple, float] = (0, 100),
    length_bounds: Union[tuple, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """Function filter_fastq

    Take a dictionary with FASTQ sequences and
    discard the sequences which do not fit the
    boundaries of GC content, length or quality threshold

    Args: seqs: dict, gc_bounds: float | tuple=(0, 100),
                 length_bounds: int | tuple=(0, 2 ** 32),
                 quality_threshold: int=0

    Returns: dict
    """
    filtered_seqs = dict()
    if type(gc_bounds) is not tuple:
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) is not tuple:
        length_bounds = (0, length_bounds)
    seqs = filter.fastq_to_dict(input_fastq)
    for name in seqs:
        sequence, quality = seqs[name]
        if (
            len(sequence) >= length_bounds[0]
            and len(sequence) <= length_bounds[1]
            and filter.gc_content(sequence) >= gc_bounds[0]
            and filter.gc_content(sequence) <= gc_bounds[1]
            and filter.seq_quality(quality) >= quality_threshold
        ):
            filtered_seqs[name] = (sequence, quality)
    filter.write_in_file(filtered_seqs, output_fastq)
    return
