import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from typing import Iterator, List, Dict
import gzip
import os
from collections import defaultdict

# 1. БАЗОВЫЕ КЛАССЫ

class Record:
    # Одна запись FASTQ: seq_id, sequence, quality
    def __init__(self, seq_id: str, sequence: str, quality: str):
        self.seq_id = seq_id
        self.sequence = sequence
        self.quality = quality

# 2. КЛАСС ДЛЯ ЧТЕНИЯ FASTQ (+ генераторы)

class FastqReader:
    def __init__(self, filename: str):
        self.filename = filename
        self.file_handle = None
        self._is_gzipped = filename.endswith('.gz')

    def __enter__(self):
        if self._is_gzipped:
            self.file_handle = gzip.open(self.filename, 'rt')
        else:
            self.file_handle = open(self.filename, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def read(self) -> Iterator[Record]:
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
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None

    @staticmethod
    def quality_to_scores(quality_str: str, phred_offset: int = 33) -> List[int]:
        return [ord(char) - phred_offset for char in quality_str]

# 3. АНАЛИЗАТОР FASTQ (статистика и графики)

class FastqAnalyzer:
    def __init__(self, filename: str, graphs_dir: str):
        self.reader = FastqReader(filename)
        self.graphs_dir = graphs_dir

    def get_sequence_count(self) -> int:
        count = 0
        with self.reader:
            for _ in self.reader.read():
                count += 1
        return count

    def get_average_sequence_length(self) -> float:
        total_length = 0
        count = 0
        with self.reader:
            for record in self.reader.read():
                total_length += len(record.sequence)
                count += 1
        return total_length / count if count > 0 else 0.0

    def get_sequence_length_distribution(self) -> Dict[int, int]:
        length_dist = defaultdict(int)
        with self.reader:
            for record in self.reader.read():
                length_dist[len(record.sequence)] += 1
        return dict(length_dist)

    def calculate_per_base_quality(self):
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

    # 4. ГРАФИКИ

    def plot_per_base_quality_matplotlib(self):
        # Per Base Sequence Quality (matplotlib)
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
        # Per Base Sequence Quality (seaborn)
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
        # Per Base Sequence Quality (Plotly interactive)
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
        # Per Base Sequence Content
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
        # Sequence Length Distribution
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

# 4. ЗАПУСК

if __name__ == "__main__":
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