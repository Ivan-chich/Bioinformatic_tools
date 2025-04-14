import unittest
import os

# Предположим, что ваши функции и классы находятся в файле your_script.py
from main_script import filter_fastq, DNASequence, SequenceError

class TestBiologicalSequences(unittest.TestCase):

    def test_dna_sequence_valid(self):
        """Тест на создание валидной ДНК последовательности."""
        seq = DNASequence("ATCG")
        self.assertEqual(str(seq), "Nucleic Acid Sequence: ATCG")

    def test_dna_sequence_invalid(self):
        """Тест на создание невалидной ДНК последовательности."""
        with self.assertRaises(SequenceError):
            DNASequence("ATCGX")

    def test_dna_sequence_empty(self):
        """Тест на создание пустой ДНК последовательности."""
        with self.assertRaises(SequenceError):
            DNASequence("")  # Ожидаем ошибку при создании пустой последовательности

    def test_dna_sequence_with_spaces(self):
        """Тест на создание ДНК последовательности с пробелами."""
        with self.assertRaises(SequenceError):
            DNASequence("A T C G")  # Ожидаем ошибку из-за пробелов

class TestFilterFastq(unittest.TestCase):

    def setUp(self):
        """Создание временного FASTQ файла для тестирования."""
        self.input_file = "test_input.fastq"
        self.output_file = "test_output.fastq"
        
        with open(self.input_file, "w") as f:
            f.write("@SEQ_ID\n")
            f.write("ACGT\n")
            f.write("+\n")
            f.write("IIII\n")  # Качество: все высокое

    def tearDown(self):
        """Удаление временных файлов после тестов."""
        if os.path.exists(self.input_file):
            os.remove(self.input_file)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def test_filter_fastq_success(self):
        """Тест на успешную фильтрацию FASTQ файла."""
        filter_fastq(self.input_file, self.output_file, min_length=1)
        
        # Проверяем, что выходной файл создан и содержит данные
        self.assertTrue(os.path.exists(self.output_file))
        
        with open(self.output_file) as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 4)  # Должно быть 4 строки (1 запись)

    def test_filter_fastq_min_length_error(self):
        """Тест на фильтрацию с ошибкой по минимальной длине."""
        filter_fastq(self.input_file, self.output_file, min_length=5)
        
        # Проверяем, что выходной файл не создан (поскольку запись была отфильтрована)
        self.assertFalse(os.path.exists(f'../{self.output_file}'))

    def test_filter_fastq_logging_error(self):
        """Тест на логирование ошибки при отсутствии входного файла."""
        with self.assertRaises(FileNotFoundError):
            filter_fastq("non_existent.fastq", self.output_file)

    def test_filter_fastq_no_records(self):
        """Тест на фильтрацию при отсутствии записей в входном файле."""
        empty_input_file = "empty_input.fastq"
        
        with open(empty_input_file, "w") as f:
            pass  # Создаем пустой файл
        
        filter_fastq(empty_input_file, self.output_file, min_length=1)
        
        # Проверяем, что выходной файл не создан (поскольку входной файл пуст)
        self.assertFalse(os.path.exists(f'../{self.output_file}'))
        
        os.remove(empty_input_file)  # Удаляем временный пустой файл

if __name__ == "__main__":
    unittest.main()