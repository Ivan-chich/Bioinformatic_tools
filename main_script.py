import argparse
import logging
import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC


# Настройка логирования
logging.basicConfig(
    filename='filter_fastq.log',  # Имя файла для логов
    level=logging.DEBUG,  # Уровень логирования
    format='%(asctime)s - %(levelname)s - %(message)s',  # Формат сообщений
    filemode='w'  # Перезаписываем файл лога при каждом запуске
)


# Настройка логирования для ошибок
error_logging = logging.getLogger('error_logger')
error_handler = logging.FileHandler('filter_fastq_errors.log', mode='w')  # Перезаписываем файл ошибокerror_handler.setLevel(logging.ERROR)
error_handler.setLevel(logging.ERROR)
error_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
error_handler.setFormatter(error_formatter)
error_logging.addHandler(error_handler)


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
        if not input_item:  # Проверка на пустую строку
            raise SequenceError("Последовательность не может быть пустой!")
        
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
        if not sequence:  # Проверка на пустую строку
            raise SequenceError("Последовательность не может быть пустой!")
        
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

    # обработка ошибки, когда не указаны имена для файлов ввода и вывода
    if not input_file or not output_file:
       error_logging.error("Не переданы аргументы input_file или output_file.")
       raise ValueError("Не переданы аргументы input_file или output_file.")
    
    # обработка ошибки, когда файл ввода не существует
    if not os.path.exists(input_file):  
       error_logging.error(f"Файл {input_file} не существует.")
       raise FileNotFoundError(f"Файл {input_file} не существует.")
    
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            # Проверка длины
            if len(record.seq) < min_length:
                logging.warning(f"Запись {record.id} пропущена из-за недостаточной длины.")
                continue

            # Проверка качества
            quality_scores = record.letter_annotations["phred_quality"]
            avg_quality = sum(quality_scores) / len(quality_scores)
            if avg_quality < min_quality:
                logging.warning(f"Запись {record.id} пропущена из-за низкого качества.")
                continue

            # Проверка GC-содержания
            gc_content = GC(record.seq)
            if gc_content < min_gc_content:
                logging.warning(f"Запись {record.id} пропущена из-за низкого GC-содержания.")
                continue

            # Если все условия выполнены, записываем в выходной файл
            SeqIO.write(record, out_handle, "fastq")

if __name__ == "__main__":
    
    # Создание парсера аргументов командной строки.
    parser = argparse.ArgumentParser(description="Фильтрация FASTQ-файлов.")
    parser.add_argument("input_file", help="Путь к входному FASTQ-файлу.")
    parser.add_argument("output_file", help="Путь к выходному FASTQ-файлу.")
    parser.add_argument("--min_length", type=int,
                        default=0,
                        help="Минимальная длина последовательности (по умолчанию: %(default)s).")
    parser.add_argument("--min_quality", type=float,
                        default=0,
                        help="Минимальное среднее качество (по умолчанию: %(default)s).")
    parser.add_argument("--min_gc_content", type=float,
                        default=0,
                        help="Минимальный процент GC-содержания (по умолчанию: %(default)s).")

    try:
       args = parser.parse_args()
    except SystemExit as e:
       error_logging.error(f"Ошибка при разборе аргументов командной строки: {e}")
       raise

    # Вызов функции фильтрации с аргументами из командной строки.
    try:
        filter_fastq(args.input_file,args.output_file,args.min_length,args.min_quality,args.min_gc_content)
    except Exception as e:
        logging.error(f"Произошла ошибка: {e}")