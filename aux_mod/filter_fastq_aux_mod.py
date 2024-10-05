def gc_content(seq: str) -> float:

    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, ndigits=1)


def seq_quality(quality_line: str) -> float:

    qscore = 0
    for symbol in quality_line:
        qscore += ord(symbol) - 33
    return round(qscore / len(quality_line), ndigits=1)
