"""
Модуль для анализа SAM файлов - формата хранения выравниваний последовательностей.
Предоставляет функциональность для чтения, парсинга и статистического анализа
данных выравнивания в SAM формате.

Основные компоненты:
- Класс SAMfile: основной класс для работы с SAM файлами
- Методы для анализа заголовков и выравниваний
- Функции для статистики по хромосомам и регионам

Входные данные: SAM файлы (.sam)
Выходные данные: статистика, таблицы с данными выравниваний
"""

import pandas as pd

class SAMfile:
    """
    КЛАСС SAMFILE - АНАЛИЗ SAM ФАЙЛОВ
    =================================
    Основной класс для работы с SAM файлами формата.
    
    Логика работы класса:
    - Чтение и разделение SAM файла на заголовки и выравнивания
    - Статистический анализ распределения выравниваний
    - Поиск выравниваний в специфичных геномных регионах
    
    Атрибуты: не требуют явной инициализации, используются внутренние переменные методов
    
    Методы:
    - parse_sam_file: разделение файла на заголовки и выравнивания
    - get_header_groups: группировка заголовков по типам
    - print_headers: вывод информации о заголовках
    - count_alignments: подсчет количества выравниваний
    - get_chromosome_stats: статистика по хромосомам
    - calculate_alignment_end: вычисление конечной позиции выравнивания
    - find_alignments_in_region: поиск выравниваний в заданном регионе
    """
    
    def parse_sam_file(self, filename):
        """
        ПАРСИНГ SAM ФАЙЛА
        =================
        Читаем файл и разделяем на заголовки и выравнивания.
        
        Логика работы:
        1. Открывает файл в текстовом режиме
        2. Итерируется по всем строкам файла
        3. Разделяет строки на заголовки (начинающиеся с '@') и выравнивания
        4. Игнорирует пустые строки
        
        Входные данные:
        - filename: строка, путь к SAM файлу
        
        Выходные данные:
        - headers: список строк заголовков
        - alignments: список строк выравниваний
        """
        headers = []
        alignments = []
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("@"):
                    headers.append(line)
                elif line:
                    alignments.append(line)
        return headers, alignments

    def get_header_groups(self, headers):
        """
        ГРУППИРОВКА ЗАГОЛОВКОВ ПО ТИПАМ
        ===============================
        Группируем заголовки по типам (@HD, @SQ, @RG, @PG).
        
        Логика работы:
        1. Для каждого заголовка определяет тип по первому полю (до первой табуляции)
        2. Создает словарь с ключами - типами заголовков
        3. Добавляет заголовки в соответствующие группы
        
        Входные данные:
        - headers: список строк заголовков
        
        Выходные данные:
        - groups: словарь {тип_заголовка: [список_заголовков]}
        """
        groups = {}
        for header in headers:
            header_type = header.split("\t")[0]
            if header_type not in groups:
                groups[header_type] = []
            groups[header_type].append(header)
        return groups

    def print_headers(self, headers):
        """
        ВЫВОД ИНФОРМАЦИИ О ЗАГОЛОВКАХ
        =============================
        Выводит структурированную информацию о заголовках SAM файла.
        
        Логика работы:
        1. Группирует заголовки по типам
        2. Для каждого типа выводит количество записей
        3. Отображает содержимое каждого заголовка
        
        Входные данные:
        - headers: список строк заголовков
        
        Выходные данные:
        - Текстовый вывод в консоль с форматированием
        """
        groups = self.get_header_groups(headers)
        print("Заголовки sam-файла:\n")
        for header_type, lines in groups.items():
            print(f"{header_type} - найдено {len(lines)} записей:")
            for line in lines:
                print(f"  {line}")
            print()

    def count_alignments(self, alignments):
        """
        ПОДСЧЕТ КОЛИЧЕСТВА ВЫРАВНИВАНИЙ
        ===============================
        Считаем общее количество выравниваний в SAM файле.
        
        Логика работы:
        - Простой подсчет элементов в списке выравниваний
        
        Входные данные:
        - alignments: список строк выравниваний
        
        Выходные данные:
        - int: количество выравниваний
        """
        return len(alignments)

    def get_chromosome_stats(self, alignments):
        """
        СТАТИСТИКА ПО ХРОМОСОМАМ
        ========================
        Анализирует распределение выравниваний по хромосомам.
        
        Логика работы:
        1. Извлекает информацию о хромосоме и позиции для каждого выравнивания
        2. Вычисляет минимальную и максимальную позицию для каждой хромосомы
        3. Создает DataFrame с подсчетом выравниваний по хромосомам
        
        Входные данные:
        - alignments: список строк выравниваний
        
        Выходные данные:
        - DataFrame с колонками: Хромосома, Количество, Мин_позиция, Макс_позиция
        """
        chromosomes = []
        positions = {}
        for alignment in alignments:
            fields = alignment.split("\t")
            chrom = fields[2]
            pos = int(fields[3])
            chromosomes.append(chrom)
            if chrom not in positions:
                positions[chrom] = {'min': pos, 'max': pos}
            else:
                positions[chrom]['min'] = min(positions[chrom]['min'], pos)
                positions[chrom]['max'] = max(positions[chrom]['max'], pos)
        df = pd.DataFrame({"Хромосома": chromosomes})
        stats = df["Хромосома"].value_counts().reset_index()
        stats.columns = ["Хромосома", "Количество"]
        # Добавляем столбцы мин/макс позиций
        stats["Мин_позиция"] = stats["Хромосома"].map(lambda c: positions[c]['min'])
        stats["Макс_позиция"] = stats["Хромосома"].map(lambda c: positions[c]['max'])
        return stats

    def calculate_alignment_end(self, position, cigar):
        """
        ВЫЧИСЛЕНИЕ КОНЕЧНОЙ ПОЗИЦИИ ВЫРАВНИВАНИЯ
        ========================================
        Вычисляет конечную позицию выравнивания на основе CIGAR строки.
        
        Логика работы:
        1. Анализирует CIGAR строку посимвольно
        2. Суммирует длины операций, которые consume референс (M, D, N, =, X)
        3. Игнорирует операции, которые не двигают позицию на референсе (I, S, H, P)
        
        Входные данные:
        - position: начальная позиция выравнивания (1-based)
        - cigar: CIGAR строка (например, "50M", "10M1I40M")
        
        Выходные данные:
        - int: конечная позиция выравнивания (1-based)
        """
        if cigar == "*":
            return position
        length = 0
        current_number = ""
        for char in cigar:
            if char.isdigit():
                current_number += char
            else:
                num = int(current_number)
                if char in "MDN=X":
                    length += num
                current_number = ""
        return position + length - 1

    def find_alignments_in_region(self, alignments, chromosome, start, end):
        """
        ПОИСК ВЫРАВНИВАНИЙ В ГЕНОМНОМ РЕГИОНЕ
        =====================================
        Находим выравнивания в заданном геномном регионе.
        
        Логика работы:
        1. Для каждого выравнивания проверяет принадлежность к целевой хромосоме
        2. Вычисляет конечную позицию выравнивания через CIGAR парсинг
        3. Проверяет пересечение с заданным регионом [start, end]
        4. Собирает информацию о подходящих выравниваниях
        
        Входные данные:
        - alignments: список строк выравниваний
        - chromosome: целевая хромосома
        - start: начальная позиция региона
        - end: конечная позиция региона
        
        Выходные данные:
        - DataFrame с колонками: Название, Хромосома, Начало, Конец, CIGAR
        """
        results = []
        for alignment in alignments:
            fields = alignment.split("\t")
            read_name = fields[0]
            chrom = fields[2]
            position = int(fields[3])
            cigar = fields[5]
            if chrom != chromosome:
                continue
            align_end = self.calculate_alignment_end(position, cigar)
            if not (align_end < start or position > end):
                results.append({
                    "Название": read_name,
                    "Хромосома": chrom,
                    "Начало": position,
                    "Конец": align_end,
                    "CIGAR": cigar,
                })
        return pd.DataFrame(results)
