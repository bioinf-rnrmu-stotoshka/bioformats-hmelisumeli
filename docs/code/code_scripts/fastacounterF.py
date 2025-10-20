class FastaAnalyzer:
    """
Этот модуль предоставляет класс для анализа FASTA файлов.
Содержит функциональность для чтения последовательностей,
подсчета статистики и генерации данных для последующего
вывода в HTML формате.

Основные компоненты:
- Класс FastaAnalyzer: главный класс анализатора
- Метод fasta_sequence_generator: генератор последовательностей
- Метод fasta_counter: подсчет статистики

Входные данные: FASTA файл с биологическими последовательностями
Выходные данные: статистика по количеству и длине последовательностей
"""

class FastaAnalyzer:
    """
    КЛАСС FASTAANALYZER
    ===================
    Основной класс для анализа FASTA файлов.
    
    Атрибуты:
    - fasta_file: путь к входному FASTA файлу
    
    Методы:
    - fasta_sequence_generator: извлекает последовательности
    - fasta_counter: вычисляет статистику
    
    Возвращает данные для HTML отчета:
    - Среднюю длину последовательностей
    - Общее количество последовательностей
    """
    
    def __init__(self, fasta_file):
        """
        КОНСТРУКТОР КЛАССА
        ==================
        Инициализирует анализатор FASTA файлов.
        
        Параметры:
        - fasta_file: строка, путь к FASTA файлу для анализа
        
        Сохраняет путь к файлу для последующей обработки
        в других методах класса.
        """
        self.fasta_file = fasta_file
    
    def fasta_sequence_generator(self):
        """
        ГЕНЕРАТОР ПОСЛЕДОВАТЕЛЬНОСТЕЙ
        ==============================
        Построчно читает FASTA файл и извлекает последовательности.
        
        Логика работы:
        1. Читает файл построчно
        2. Определяет заголовки (строки начинающиеся с '>')
        3. Собирает последовательности до следующего заголовка
        4. Возвращает кортежи (заголовок, последовательность)
        
        Выходные данные:
        - header: строка, название последовательности (без '>')
        - sequence: строка, полная биологическая последовательность
        
        Используется для последующего анализа и HTML визуализации.
        """
        current_header = None
        current_sequence = []
        
        with open(self.fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:  # пропускаем пустые строки
                    continue
                    
                if line.startswith('>'):
                    if current_header is not None:
                        full_sequence = ''.join(current_sequence)
                        yield current_header, full_sequence
                        current_sequence = []  # сбрасываем для новой последовательности
                    
                    current_header = line[1:]  # убираем '>'
                else:
                    current_sequence.append(line)
            
            #обрабатываем последнюю последовательность отдельно, потому что нет тригера >
            if current_header is not None and current_sequence:
                full_sequence = ''.join(current_sequence)
                yield current_header, full_sequence
    
    def fasta_counter(self):
        """
        СТАТИСТИЧЕСКИЙ АНАЛИЗАТОР
        ==========================
        Вычисляет базовую статистику по FASTA файлу.
        
        Вычисления:
        - total_length: общая длина всех последовательностей
        - sequence_count: количество последовательностей (n)
        - average_length: средняя длина последовательности
        
        Возвращаемые данные для HTML отчета:
        - average_length: float, средняя длина (для отображения в таблице)
        - sequence_count: int, количество (для отображения в статистике)
        
        Использует генератор для эффективной обработки больших файлов.
        """
        total_length = 0 
        sequence_count = 0  #n
        
        #исп генератор вместо прямого чтения
        for header, sequence in self.fasta_sequence_generator():
            total_length += len(sequence)
            sequence_count += 1
        
        #вычисляем среднюю длину если есть последовательности
        if sequence_count > 0:
            average_length = total_length / sequence_count
        else:
            average_length = 0
        
        return average_length, sequence_count


"""
HTML TEMPLATE INTEGRATION READY
===============================
Данный модуль готов для интеграции с HTML шаблонизаторами:

Возможные варианты использования:
1. Jinja2 templates: {{ average_length }}, {{ sequence_count }}
2. Django templates: {{ context_data }}
3. Прямое форматирование в HTML строки

Пример HTML вывода:
- Количество последовательностей: <span class="count">{sequence_count}</span>
- Средняя длина: <span class="avg-length">{average_length:.2f}</span>

Для использования с HTML добавьте импорт:
from jinja2 import Template  # для шаблонизации
"""

#использование
#if __name__ == "__main__":
    #analyzer = FastaAnalyzer("example.fasta")
    
    #статистика
    #average_length, sequence_count = analyzer.fasta_counter()
    
    #print(f"Количество последовательностей: {sequence_count}")
    #print(f"Средняя длина последовательности: {average_length:.2f}")