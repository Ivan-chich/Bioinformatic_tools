from typing import TextIO, Union


def convert_multiline_fasta_to_oneline(
    input_fasta: Union[str | TextIO], output_fasta: Union[str | TextIO]
) -> None:
    """
    Function convert_multiline_fasta_to_oneline

    Accepts 2 arguments as input: input_fasta and output_fasta.
    Reads the fasta file submitted as input_fasta, in which the
    sequence (DNA/RNA/protein/ ... ) can be split into several
    lines, after which it saves it to a new fasta file in which
    each sequence fits into one line

    Arguments: input_fasta: 'input_file_name',
    output_fasta: 'output_file_name'

    Returns None, creates output .fasta file
    """
    with open(input_fasta, "r", encoding="utf-8") as file:
        lst = []
        for line in file:
            lst.append(line.strip())
    with open(f"{output_fasta}.fasta", "a", encoding="utf-8") as output:
        s = ""
        print(lst[0], file=output)
        for i in range(1, len(lst)):

            if lst[i].startswith(">"):
                print(s, file=output)
                print(lst[i], file=output)
                s = ""
            else:
                s += lst[i]
        print(s, file=output)
    return


def parse_blast_output(input_file: TextIO, output_file: TextIO) -> None:
    """
    Function parse_blast_output

    Accepts 2 arguments as input: input_file and output_file.
    The function reads the specified file, selects the first
    row from the Description column for each QUERY (Sequences
    producing significant alignments). Saves the list of obtained
    proteins to a new file in one column sorted alphabetically.

    Arguments: input_file: 'input_file_name',
    output_file: 'output_file_name'

    Returns None, creates output .txt file
    """
    lst = []
    with open(input_file, "r", encoding="utf-8") as file:
        for line in file:
            if "Description" in line:
                s = file.readline().strip().split("...")[0]
                if "  " in s:
                    lst.append(s.split("  ")[0])
                else:
                    lst.append(s)

    with open(f"{output_file}.txt", "a", encoding="utf-8") as file:
        for line in sorted([item.capitalize() for item in lst]):
            print(line, file=file)
    return
