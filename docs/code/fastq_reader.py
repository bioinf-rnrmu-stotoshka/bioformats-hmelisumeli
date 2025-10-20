from jinja2 import Template, Environment, FileSystemLoader
import json
"""
Модуль для комплексного анализа FASTQ файлов с генерацией графиков и статистики.

Основные компоненты:
- Record: хранение одной записи FASTQ
- FastqReader: чтение и парсинг FASTQ файлов
- FastqAnalyzer: анализ данных и визуализация

Входные данные: FASTQ файлы (.fastq, .fq, .gz)
Выходные данные: статистика + графики в PNG/HTML форматах

Логика работы модуля:
1. Чтение FASTQ файла построчно с обработкой gzip сжатия
2. Парсинг записей на компоненты: идентификатор, последовательность, качество
3. Статистический анализ: подсчет, распределение длин, качество
4. Визуализация результатов в multiple форматах
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from typing import Iterator, List, Dict
import gzip
import os
from collections import defaultdict

class Record:
    """
    КЛАСС RECORD - ОДНА ЗАПИСЬ FASTQ
    ================================
    Хранит информацию об одной последовательности из FASTQ файла.
    
    Логика работы:
    - Принимает три строковых параметра при инициализации
    - Сохраняет их как атрибуты объекта для последующего анализа
    
    Атрибуты:
    - seq_id: идентификатор последовательности (без '@')
    - sequence: биологическая последовательность (ATCGN)
    - quality: строка качества в символьном формате
    
    Выходные данные:
    - Объект с структурированными данными одной прочтения
    """
    def __init__(self, seq_id: str, sequence: str, quality: str):
        self.seq_id = seq_id
        self.sequence = sequence
        self.quality = quality

class FastqReader:
    """
    КЛАСС FASTQ READER - ЧТЕНИЕ FASTQ ФАЙЛОВ
    ========================================
    Обеспечивает чтение FASTQ файлов через контекстный менеджер.
    
    Логика работы:
    - Автоматическое определение типа файла (gzip/plain text)
    - Построчное чтение с обработкой структуры FASTQ формата
    - Использование генератора для эффективной работы с большими файлами
    
    Особенности:
    - Автоматически определяет gzip сжатие
    - Использует генераторы для больших файлов
    - Контекстный менеджер для безопасной работы
    
    Выходные данные:
    - Генератор объектов Record для последовательной обработки
    """
    def __init__(self, filename: str):
        self.filename = filename
        self.file_handle = None
        self._is_gzipped = filename.endswith('.gz')

    def __enter__(self):
        """
        ОТКРЫТИЕ ФАЙЛА В КОНТЕКСТНОМ МЕНЕДЖЕРЕ
        =====================================
        Логика работы:
        - Проверяет расширение файла для определения типа
        - Открывает соответствующим методом (gzip/open)
        - Возвращает self для использования в with-блоке
        """
        if self._is_gzipped:
            self.file_handle = gzip.open(self.filename, 'rt')
        else:
            self.file_handle = open(self.filename, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        ЗАКРЫТИЕ ФАЙЛА ПРИ ВЫХОДЕ ИЗ КОНТЕКСТА
        =====================================
        Логика работы:
        - Автоматически вызывается при выходе из with-блока
        - Гарантирует закрытие файлового дескриптора
        - Обрабатывает исключения без потери ресурсов
        """
        self.close()

    def read(self) -> Iterator[Record]:
        """
        ГЕНЕРАТОР ЗАПИСЕЙ FASTQ
        =======================
        Построчно читает файл и возвращает Record объекты.
        
        Логика работы:
        1. Читает заголовок (начинается с '@')
        2. Читает последовательность
        3. Пропускает '+' строку
        4. Читает строку качества
        5. Создает объект Record и возвращает через yield
        
        Алгоритм:
        - Чтение продолжается до пустой строки (конец файла)
        - Обработка каждой четверки строк как одной записи
        - Автоматическое определение конца файла
        
        Выходные данные:
        - Iterator[Record]: поток объектов для последовательной обработки
        """
        if self.file_handle is None:
            raise RuntimeError("Use 'with FastqReader(...)'")
        while True:
            header = self.file_handle.readline().strip()
            if not header:
                break
            sequence = self.file_handle.readline().strip()
            self.file_handle.readline()  # плюс-строка
            quality = self.file_handle.readline().strip()
            seq_id = header[1:] if header.startswith('@') else header
            yield Record(seq_id, sequence, quality)

    def close(self):
        """
        РУЧНОЕ ЗАКРЫТИЕ ФАЙЛА
        =====================
        Логика работы:
        - Проверяет наличие открытого файлового дескриптора
        - Закрывает файл и сбрасывает указатель на None
        - Предотвращает утечки ресурсов
        """
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None

    @staticmethod
    def quality_to_scores(quality_str: str, phred_offset: int = 33) -> List[int]:
        """
        ПРЕОБРАЗОВАНИЕ КАЧЕСТВА В Phred SCORES
        ======================================
        Конвертирует символьное представление качества в числовые значения.
        
        Логика работы:
        - Для каждого символа в строке качества вычисляет числовое значение
        - Использует формулу Phred: score = ASCII(code) - offset
        - По умолчанию используется Phred+33 кодировка
        
        Формула: Phred score = ord(char) - phred_offset
        
        Выходные данные:
        - List[int]: список числовых значений качества для каждой позиции
        - Качество 20 = 1% ошибка, 30 = 0.1% ошибка, 40 = 0.01% ошибка
        """
        return [ord(char) - phred_offset for char in quality_str]

class FastqAnalyzer:
    """
    КЛАСС FASTQ ANALYZER - АНАЛИЗ И ВИЗУАЛИЗАЦИЯ
    ============================================
    Выполняет полный анализ FASTQ файлов и генерирует графики.
    
    Логика работы:
    - Использует FastqReader для чтения данных
    - Применяет статистические методы для анализа
    - Генерирует multiple форматы визуализации
    - Сохраняет результаты в указанную директорию
    
    Основные функции:
    - Подсчет статистики (количество, средняя длина)
    - Анализ качества по позициям
    - Анализ нуклеотидного состава
    - Генерация графиков в разных форматах
    
    Выходные данные:
    - Статистические показатели (числа, словари)
    - Графические файлы (PNG, HTML)
    - Текстовые отчеты о выполнении
    """
    def __init__(self, filename: str, graphs_dir: str):
        self.reader = FastqReader(filename)
        self.graphs_dir = graphs_dir

    def get_sequence_count(self) -> int:
        """
        ПОДСЧЕТ КОЛИЧЕСТВА ПОСЛЕДОВАТЕЛЬНОСТЕЙ
        ======================================
        Возвращает общее количество записей в FASTQ файле.
        
        Логика работы:
        - Итерируется по всем записям через генератор
        - Инкрементирует счетчик для каждой записи
        - Возвращает итоговое количество
        
        Выходные данные:
        - int: общее количество последовательностей в файле
        """
        count = 0
        with self.reader:
            for _ in self.reader.read():
                count += 1
        return count

    def get_average_sequence_length(self) -> float:
        """
        ВЫЧИСЛЕНИЕ СРЕДНЕЙ ДЛИНЫ ПОСЛЕДОВАТЕЛЬНОСТИ
        ===========================================
        Рассчитывает среднюю длину всех последовательностей.
        
        Логика работы:
        - Суммирует длины всех последовательностей
        - Подсчитывает общее количество
        - Вычисляет среднее арифметическое
        
        Выходные данные:
        - float: средняя длина последовательности в базовых парах
        - 0.0 если файл пустой
        """
        total_length = 0
        count = 0
        with self.reader:
            for record in self.reader.read():
                total_length += len(record.sequence)
                count += 1
        return total_length / count if count > 0 else 0.0

    def get_sequence_length_distribution(self) -> Dict[int, int]:
        """
        РАСПРЕДЕЛЕНИЕ ДЛИН ПОСЛЕДОВАТЕЛЬНОСТЕЙ
        ======================================
        Создает словарь с распределением длин последовательностей.
        
        Логика работы:
        - Для каждой последовательности определяет длину
        - Увеличивает счетчик для соответствующей длины
        - Использует defaultdict для автоматической инициализации
        
        Выходные данные:
        - Dict[int, int]: словарь где ключ - длина, значение - количество
        - Пример: {75: 1500, 76: 4500, 100: 200}
        """
        length_dist = defaultdict(int)
        with self.reader:
            for record in self.reader.read():
                length_dist[len(record.sequence)] += 1
        return dict(length_dist)

    def calculate_per_base_quality(self):
        """
        РАСЧЕТ КАЧЕСТВА ПО ПОЗИЦИЯМ
        ===========================
        Анализирует качество последовательностей для каждой позиции.
        
        Логика работы:
        1. Для каждой записи преобразует строку качества в числовые scores
        2. Группирует scores по позициям в последовательности
        3. Вычисляет статистику для каждой позиции
        
        Выходные данные:
        - mean_qualities: список средних значений качества по позициям
        - all_qualities: список списков со всеми значениями для каждой позиции
        - Позиции нумеруются с 0 до максимальной длины
        """
        position_qualities = defaultdict(list)
        with self.reader:
            for record in self.reader.read():
                qualities = FastqReader.quality_to_scores(record.quality)
                for pos, qual in enumerate(qualities):
                    position_qualities[pos].append(qual)
        max_pos = max(position_qualities.keys()) if position_qualities else 0
        mean_qualities, all_qualities = [], []
        for pos in range(max_pos + 1):
            quals = position_qualities[pos]
            mean_qualities.append(np.mean(quals) if quals else 0)
            all_qualities.append(quals)
        return mean_qualities, all_qualities

    def calculate_per_base_content(self) -> Dict[str, List[float]]:
        """
        АНАЛИЗ НУКЛЕОТИДНОГО СОСТАВА ПО ПОЗИЦИЯМ
        ========================================
        Рассчитывает процентное содержание A,C,G,T для каждой позиции.
        
        Логика работы:
        1. Для каждой позиции подсчитывает количество каждого нуклеотида
        2. Вычисляет процентное соотношение относительно общего количества
        3. Игнорирует нестандартные нуклеотиды (кроме N)
        
        Выходные данные:
        - Dict[str, List[float]]: словарь с ключами 'A','C','G','T'
        - Каждый ключ содержит список процентов для каждой позиции
        - Пример: {'A': [25.5, 26.1, ...], 'C': [24.8, 24.5, ...]}
        """
        position_counts = defaultdict(lambda: {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'total': 0})
        with self.reader:
            for record in self.reader.read():
                for pos, base in enumerate(record.sequence.upper()):
                    if base in 'ACGTN':
                        position_counts[pos][base] += 1
                        position_counts[pos]['total'] += 1
        max_pos = max(position_counts.keys()) if position_counts else 0
        result = {'A': [], 'C': [], 'G': [], 'T': []}
        for pos in range(max_pos + 1):
            counts = position_counts[pos]
            total = counts['total']
            for base in ['A', 'C', 'G', 'T']:
                result[base].append((counts[base] / total) * 100 if total > 0 else 0)
        return result

    def plot_per_base_quality_matplotlib(self):
        """
        ГРАФИК КАЧЕСТВА ПО ПОЗИЦИЯМ (MATPLOTLIB)
        ========================================
        Создает детализированный boxplot-подобный график качества.
        
        Логика работы:
        - Вычисляет перцентили (25%, медиана, 75%, 10%, 90%) для каждой позиции
        - Рисует прямоугольники для межквартильного размаха
        - Добавляет линии для медианы и процентилей
        - Размечает цветовые зоны качества (зеленый/оранжевый/красный)
        
        Выходные данные:
        - PNG файл с детальным графиком качества
        - Визуализация распределения качества по позициям
        """
        mean_qualities, all_qualities = self.calculate_per_base_quality()
        if not mean_qualities:
            print("Нет данных для графика")
            return
        positions = list(range(1, len(mean_qualities) + 1))
        percentiles = []
        for quals in all_qualities:
            if quals:
                percentiles.append({
                    'q25': np.percentile(quals, 25),
                    'median': np.median(quals),
                    'q75': np.percentile(quals, 75),
                    'q10': np.percentile(quals, 10),
                    'q90': np.percentile(quals, 90),
                    'mean': np.mean(quals)
                })
            else:
                percentiles.append({k: 0 for k in ['q25', 'median', 'q75', 'q10', 'q90', 'mean']})
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.axhspan(28, 42, facecolor='green', alpha=0.1)
        ax.axhspan(20, 28, facecolor='orange', alpha=0.1)
        ax.axhspan(0, 20, facecolor='red', alpha=0.1)
        box_width = 0.8
        for pos, perc in zip(positions, percentiles):
            ax.add_patch(plt.Rectangle((pos - box_width/2, perc['q25']),
                                      box_width, perc['q75'] - perc['q25'],
                                      facecolor='yellow', edgecolor='black', linewidth=0.5))
            ax.plot([pos - box_width/2, pos + box_width/2],
                   [perc['median'], perc['median']], 'r-', linewidth=1)
            ax.plot([pos, pos], [perc['q75'], perc['q90']], 'k-', linewidth=0.5)
            ax.plot([pos, pos], [perc['q25'], perc['q10']], 'k-', linewidth=0.5)
        ax.plot(positions, [p['mean'] for p in percentiles], 'b-', linewidth=1, label='Среднее')
        ax.set_xlabel('Позиция в прочтении (bp)', fontsize=12)
        ax.set_ylabel('Phred качество', fontsize=12)
        ax.set_title('Per Base Sequence Quality', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 42)
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        output_path = os.path.join(self.graphs_dir, 'per_base_quality_matplotlib.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"График сохранен: {output_path}")

    def plot_per_base_quality_seaborn(self):
        """
        ГРАФИК КАЧЕСТВА ПО ПОЗИЦИЯМ (SEABORN)
        =====================================
        Создает чистый линейный график среднего качества.
        
        Логика работы:
        - Использует стиль whitegrid от seaborn для чистого вида
        - Рисует линию среднего качества по позициям
        - Добавляет цветовые зоны для интерпретации качества
        
        Выходные данные:
        - PNG файл с упрощенным графиком качества
        - Чистая визуализация тренда качества
        """
        mean_qualities, _ = self.calculate_per_base_quality()
        if not mean_qualities:
            print("Нет данных для графика")
            return
        positions = list(range(1, len(mean_qualities) + 1))
        sns.set_style("whitegrid")
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.axhspan(28, 42, facecolor='green', alpha=0.1)
        ax.axhspan(20, 28, facecolor='orange', alpha=0.1)
        ax.axhspan(0, 20, facecolor='red', alpha=0.1)
        ax.plot(positions, mean_qualities, linewidth=2, color='blue', label='Среднее качество')
        ax.set_xlabel('Позиция в прочтении (bp)', fontsize=12)
        ax.set_ylabel('Phred качество', fontsize=12)
        ax.set_title('Per Base Sequence Quality (Seaborn)', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 42)
        ax.legend()
        plt.tight_layout()
        output_path = os.path.join(self.graphs_dir, 'per_base_quality_seaborn.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"График сохранен: {output_path}")

    def plot_per_base_quality_plotly(self):
        """
        ИНТЕРАКТИВНЫЙ ГРАФИК КАЧЕСТВА (PLOTLY)
        ======================================
        Создает полностью интерактивный HTML график.
        
        Логика работы:
        - Создает multiple traces для среднего, медианы и IQR
        - Добавляет интерактивные элементы (hover, zoom, pan)
        - Использует заполненные области для визуализации распределения
        
        Выходные данные:
        - HTML файл с интерактивным графиком
        - Возможность детального изучения через браузер
        """
        mean_qualities, all_qualities = self.calculate_per_base_quality()
        if not mean_qualities:
            print("Нет данных для графика")
            return
        positions = list(range(1, len(mean_qualities) + 1))
        medians = [np.median(q) if q else 0 for q in all_qualities]
        q25 = [np.percentile(q, 25) if q else 0 for q in all_qualities]
        q75 = [np.percentile(q, 75) if q else 0 for q in all_qualities]
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=positions, y=mean_qualities, mode='lines', name='Среднее',
                                line=dict(color='blue', width=2)))
        fig.add_trace(go.Scatter(x=positions, y=medians, mode='lines', name='Медиана',
                                line=dict(color='red', width=2)))
        fig.add_trace(go.Scatter(x=positions + positions[::-1], y=q75 + q25[::-1],
                                fill='toself', fillcolor='rgba(255, 255, 0, 0.3)',
                                line=dict(color='rgba(255,255,255,0)'), name='IQR (25-75%)'))
        fig.add_hrect(y0=28, y1=42, fillcolor="green", opacity=0.1, line_width=0)
        fig.add_hrect(y0=20, y1=28, fillcolor="orange", opacity=0.1, line_width=0)
        fig.add_hrect(y0=0, y1=20, fillcolor="red", opacity=0.1, line_width=0)
        fig.update_layout(title='Per Base Sequence Quality (Interactive)',
                         xaxis_title='Позиция в прочтении (bp)',
                         yaxis_title='Phred качество',
                         yaxis=dict(range=[0, 42]),
                         template='plotly_white', height=600)
        output_path = os.path.join(self.graphs_dir, 'per_base_quality_plotly.html')
        fig.write_html(output_path)
        print(f"Интерактивный график сохранен: {output_path}")

    def plot_per_base_content(self):
        """
        ГРАФИК НУКЛЕОТИДНОГО СОСТАВА ПО ПОЗИЦИЯМ
        ========================================
        Визуализирует процентное содержание A,C,G,T вдоль последовательностей.
        
        Логика работы:
        - Рисует отдельные линии для каждого нуклеотида
        - Использует стандартные цвета для идентификации (A-зеленый, C-синий и т.д.)
        - Показывает динамику состава вдоль последовательности
        
        Выходные данные:
        - PNG файл с графиком нуклеотидного состава
        - Визуализация bias в составе последовательностей
        """
        content = self.calculate_per_base_content()
        if not content['A']:
            print("Нет данных для графика")
            return
        positions = list(range(1, len(content['A']) + 1))
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.plot(positions, content['A'], label='% A', color='green', linewidth=1.5)
        ax.plot(positions, content['C'], label='% C', color='blue', linewidth=1.5)
        ax.plot(positions, content['G'], label='% G', color='black', linewidth=1.5)
        ax.plot(positions, content['T'], label='% T', color='red', linewidth=1.5)
        ax.set_xlabel('Позиция в прочтении (bp)', fontsize=12)
        ax.set_ylabel('Процент (%)', fontsize=12)
        ax.set_title('Per Base Sequence Content', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        output_path = os.path.join(self.graphs_dir, 'per_base_content.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"График сохранен: {output_path}")

    def plot_sequence_length_distribution(self):
        """
        ГИСТОГРАММА РАСПРЕДЕЛЕНИЯ ДЛИН ПОСЛЕДОВАТЕЛЬНОСТЕЙ
        ==================================================
        Показывает вариабельность длин последовательностей.
        
        Логика работы:
        - Сортирует длины последовательностей
        - Создает столбчатую диаграмму с количеством для каждой длины
        - Добавляет сетку для удобства чтения
        
        Выходные данные:
        - PNG файл с гистограммой распределения длин
        - Визуализация модальной длины и вариабельности
        """
        length_dist = self.get_sequence_length_distribution()
        if not length_dist:
            print("Нет данных для графика")
            return
        lengths = sorted(length_dist.keys())
        counts = [length_dist[l] for l in lengths]
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(lengths, counts, color='steelblue', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Длина последовательности (bp)', fontsize=12)
        ax.set_ylabel('Количество', fontsize=12)
        ax.set_title('Sequence Length Distribution', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        output_path = os.path.join(self.graphs_dir, 'sequence_length_distribution.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"График сохранен: {output_path}")

if __name__ == "__main__":
    """
    ОСНОВНАЯ ПРОГРАММА - ДЕМОНСТРАЦИЯ ВОЗМОЖНОСТЕЙ
    ==============================================
    Запрашивает путь к FASTQ файлу и выполняет полный анализ.
    
    Логика работы:
    1. Запрос пути к файлу у пользователя
    2. Проверка существования файла
    3. Создание директории для графиков
    4. Вычисление базовой статистики
    5. Генерация всех типов графиков
    
    Выходные данные:
    - Статистика в консоли
    - Графические файлы в папке graphs
    - Подробный отчет о выполнении
    """
    print("Демоверсия анализа FASTQ файлов")
    FASTQ_FILE = input("Введите полный путь к FASTQ файлу: ").strip().strip('"').strip("'")
    if not os.path.isfile(FASTQ_FILE):
        print(f"\nОШИБКА: Файл не найден: {FASTQ_FILE}")
        exit(1)
    graphs_dir = os.path.join(os.path.dirname(os.path.abspath(FASTQ_FILE)), "graphs")
    if not os.path.exists(graphs_dir):
        os.makedirs(graphs_dir)
        print(f"\nСоздана папка для графиков: {graphs_dir}")
    analyzer = FastqAnalyzer(FASTQ_FILE, graphs_dir)
    print(f"\nВычисление базовой статистики...")
    print(f"  Количество последовательностей: {analyzer.get_sequence_count()}")
    print(f"  Средняя длина последовательности: {analyzer.get_average_sequence_length():.2f} bp")
    print("\nГенерация графиков...")
    analyzer.plot_per_base_quality_matplotlib()
    analyzer.plot_per_base_quality_seaborn()
    analyzer.plot_per_base_quality_plotly()
    analyzer.plot_per_base_content()
    analyzer.plot_sequence_length_distribution()
    print(f"\nАнализ завершен!\nВсе графики сохранены в папке: {graphs_dir}")