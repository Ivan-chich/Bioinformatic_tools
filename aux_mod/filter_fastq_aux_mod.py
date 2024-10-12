import os
from typing import Union, TextIO


def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def seq_quality(quality_line: str) -> float:
    qscore = 0
    for symbol in quality_line:
        qscore += ord(symbol) - 33
    return round(qscore / len(quality_line), ndigits=1)


def fastq_to_dict(file: TextIO) -> dict:
    lst = []
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            lst.append(line.strip())
    dct = dict()
    for i in range(0, len(lst), 4):
        name = lst[i]
        sequence = lst[i + 1]
        quality = lst[i + 3]
        dct[name] = (sequence, quality)
    return dct


def write_in_file(dct: dict,
                  output_file_name: Union[str, os.PathLike]) -> None:
    filename = f"./filtered/{output_file_name}"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "a", encoding="utf-8") as output:
        print(dct, file=output)
