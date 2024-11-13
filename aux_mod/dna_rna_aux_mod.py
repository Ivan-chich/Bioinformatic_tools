dct_dna = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}
dct_rna = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
}


def is_dna(seq: str) -> str:
    return set(seq.upper()).issubset({"A", "T", "G", "C"})


def is_rna(seq: str) -> str:
    return set(seq.upper()).issubset({"A", "U", "G", "C"})


def transcribe(seq: str) -> str:
    for i in (("T", "U"), ("t", "u")):
        seq = seq.replace(*i)
    return seq


def reverse(seq: str) -> str:
    return seq[::-1]


def complement(seq: str) -> str:
    bases = list(seq)
    if is_dna(seq) is True:
        bases = [dct_dna[letter] for letter in bases]
    if is_rna(seq) is True:
        bases = [dct_rna[letter] for letter in bases]
    return "".join(bases)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]
