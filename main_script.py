from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC


class SequenceError(Exception):
    """Ошибка, возникающая при вводе неправильной последовательности."""

    pass


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        """Возвращает длину последовательности."""
        pass

    @abstractmethod
    def __getitem__(self, index):
        """Получает элемент по индексу или срез последовательности."""
        pass

    @abstractmethod
    def __str__(self):
        """Возвращает строковое представление последовательности."""
        pass

    @abstractmethod
    def __is_valid__(self):
        """Проверяет последовательность на недопустимые символы."""
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, input_item):
        self.sequence = input_item.upper()  # Приводим к верхнему регистру
        if type(input_item) not in ["__main__.DNASequence", "__main__.RNASequence"]:
            raise NotImplementedError(
                "Определите последовательность через класс DNASequence или RNASequence"
            )

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return f"Nucleic Acid Sequence: {self.sequence}"

    def __is_valid__(sequence, alphabet="AUCGT"):
        return all(base in alphabet for base in sequence)

    def complement(self):
        """Возвращает комплементарную последовательность."""
        complement_map = str.maketrans("ACGTU", "TGCAA")
        return self.sequence.translate(complement_map)

    def reverse(self):
        """Возвращает обратную последовательность."""
        return self.sequence[::-1]

    def reverse_complement(self):
        """Возвращает обратную комплементарную последовательность."""
        return self.reverse().translate(str.maketrans("ACGTU", "TGCAA"))


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        self.sequence = sequence.upper()  # Приводим к верхнему регистру
        if not NucleicAcidSequence.__is_valid__(sequence, alphabet="ACGT"):
            raise SequenceError("В последовательности ДНК есть недопустимые символы!")

    def transcribe(self):
        """Возвращает транскрибированную РНК-последовательность."""
        return self.sequence.replace("T", "U")


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        self.sequence = sequence.upper()  # Приводим к верхнему регистру
        if not NucleicAcidSequence.__is_valid__(sequence, alphabet="ACGU"):
            raise SequenceError("В последовательности РНК есть недопустимые символы!")

    def reverse_transcribe(self):
        """Возвращает обратно транскрибированную ДНК-последовательность."""
        return self.sequence.replace("U", "T")


class AminoAcidSequence(BiologicalSequence):
    # Молекулярные массы аминокислот (в г/моль)
    amino_acid_masses = {
        "A": 89.09,  # Alanine
        "C": 121.15,  # Cysteine
        "D": 133.10,  # Aspartic acid
        "E": 147.13,  # Glutamic acid
        "F": 165.19,  # Phenylalanine
        "G": 75.07,  # Glycine
        "H": 155.16,  # Histidine
        "I": 131.17,  # Isoleucine
        "K": 146.19,  # Lysine
        "L": 131.17,  # Leucine
        "M": 149.21,  # Methionine
        "N": 132.12,  # Asparagine
        "P": 115.13,  # Proline
        "Q": 146.15,  # Glutamine
        "R": 174.20,  # Arginine
        "S": 105.09,  # Serine
        "T": 119.12,  # Threonine
        "V": 117.15,  # Valine
        "W": 204.23,  # Tryptophan
        "Y": 181.19,  # Tyrosine
    }

    def __init__(self, sequence):
        self.sequence = sequence.upper()  # Приводим к верхнему регистру
        if not all(
            residue in AminoAcidSequence.amino_acid_masses.keys()
            for residue in self.sequence
        ):
            raise SequenceError("В последовательности белка есть недопустимые символы!")

    def __is_valid__(self):
        return all(
            residue in AminoAcidSequence.amino_acid_masses.keys()
            for residue in self.sequence
        )

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return f"Amino Acid Sequence: {self.sequence}"

    def molecular_weight(self):
        """Возвращает молекулярную массу аминокислотной последовательности."""
        return round(sum(self.amino_acid_masses.get(aa, 0) for aa in self.sequence))


def filter_fastq(
    input_file, output_file, min_length=0, min_quality=0, min_gc_content=0
):
    """
    Фильтрует записи из FASTQ-файла по длине, качеству и GC-составу.

    :param input_file: Путь к входному FASTQ-файлу.
    :param output_file: Путь к выходному FASTQ-файлу.
    :param min_length: Минимальная длина последовательности.
    :param min_quality: Минимальное среднее качество.
    :param min_gc_content: Минимальный процент GC-содержания (в диапазоне от 0 до 100).
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            # Проверка длины
            if len(record.seq) < min_length:
                continue

            # Проверка качества
            quality_scores = record.letter_annotations["phred_quality"]
            avg_quality = sum(quality_scores) / len(quality_scores)
            if avg_quality < min_quality:
                continue

            # Проверка GC-содержания
            gc_content = GC(record.seq)
            if gc_content < min_gc_content:
                continue

            # Если все условия выполнены, записываем в выходной файл
            SeqIO.write(record, out_handle, "fastq")

