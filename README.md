# Bioinformatic tools

The file `main_script.py` contains classes and functions for work with
nucleic acid and amino acid sequences and filtering FASTQ files

## Classes:

### Abstract class BiologicalSequence

This class includes abstract methods:

- `__len__`: Возвращает длину последовательности.
- `__getitem__`: Получает элемент по индексу или срез последовательности.
- `__str__`: Возвращает строковое представление последовательности.
- `__is_valid__`: Проверяет последовательность на недопустимые символы.

### Class NucleicAcidSequence

This class is inherent to BiologicalSequence. Includes the following methods:

- `__init__`
- `__len__`
- `__getitem__`
- `__str__`
- `__is_valid__`
- `complement`: Возвращает комплементарную последовательность
- `reverse`: Возвращает обратную последовательность
- `reverse_complement`: Возвращает обратную комплементарную последовательность

### Class DNASequence

This class is inherent to NucleicAcidSequence. Accepts nucleotide sequence
as argument. Includes the following methods:

- `__init__`: Осуществляет проверку входящей последовательности на недопустимые символы. Если находит, выдаёт ошибку.
- `transcribe`: Возвращает транскрибированную РНК-последовательность

### Class RNASequence

This class is inherent to NucleicAcidSequence. Accepts nucleotide sequence
as argument. Includes the following methods:

- `__init__`: Осуществляет проверку входящей последовательности на недопустимые символы. Если находит, выдаёт ошибку.
- `reverse_transcribe`: Возвращает транскрибированную ДНК-последовательность

### Class AminoAcidSequence

This class is inherent to BiologicalSequence. Accepts amino acid sequence
as argument. Includes the following methods:

- `__init__`: Осуществляет проверку входящей последовательности на недопустимые символы. Если находит, выдаёт ошибку.
- `__len__`
- `__getitem__`
- `__str__`
- `__is_valid__`
- `molecular_weight`: Возвращает молекулярную массу аминокислотной последовательности

This class has one attribute:

- `amino_acid_masses`: словарь, содержащий однобуквенный код аминокислот в качестве ключей и соответствующие молекулярные массы в качестве значений.

## Functions:

### filter_fastq

Фильтрует записи из FASTQ-файла по длине, качеству и GC-составу. Записывает фильтрованные записи в новый файл.

- :param `input_file`: Путь к входному FASTQ-файлу.
- :param `output_file`: Путь к выходному FASTQ-файлу.
- :param `min_length`: Минимальная длина последовательности.
- :param `min_quality`: Минимальное среднее качество.
- :param `min_gc_content`: Минимальный процент GC-содержания (в диапазоне от 0 до 100).
