from typing import Iterator, List, Dict
import gzip
import os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

# 1. БАЗОВЫЕ КЛАССЫ

class Record:
    # Одна запись FASTQ: seq_id, sequence, quality
    def __init__(self, seq_id: str, sequence: str, quality: str):
        self.seq_id = seq_id
        self.sequence = sequence
        self.quality = quality

# 2. КЛАСС ДЛЯ ЧТЕНИЯ FASTQ + генераторы

class FastqReader:
    # Инициализация reader’а (учёт .gz)
    def __init__(self, filename: str):
        self.filename = filename
        self.file_handle = None
        self._is_gzipped = filename.endswith('.gz')

    def __enter__(self):
        # Открытие файла в текстовом режиме
        if self._is_gzipped:
            self.file_handle = gzip.open(self.filename, 'rt')
        else:
            self.file_handle = open(self.filename, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Закрытие файла
        self.close()

    def read(self) -> Iterator[Record]:
        # Генератор чтения: header, sequence, +, quality
        if self.file_handle is None:
            raise RuntimeError("Use 'with FastqReader(...)'")
        while True:
            header = self.file_handle.readline().strip()
            if not header:
                break
            seq = self.file_handle.readline().strip()
            self.file_handle.readline()  # плюс
            qual = self.file_handle.readline().strip()
            seq_id = header[1:] if header.startswith('@') else header
            yield Record(seq_id, seq, qual)

    def close(self):
        # Закрытие дескриптора
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None

    @staticmethod
    def _quality_string_to_scores(qs: str, phred_offset: int = 33) -> List[int]:
        # ASCII → Phred scores
        return [ord(c) - phred_offset for c in qs]

# 3. АНАЛИЗАТОР FASTQ (статистика)

class FastqAnalyzer:
    def __init__(self, filename: str, graphs_dir: str):
        self.reader = FastqReader(filename)
        self.graphs_dir = graphs_dir

    def get_sequence_count(self) -> int:
        # Количество прочтений
        cnt = 0
        with self.reader:
            for _ in self.reader.read():
                cnt += 1
        return cnt

    def get_average_sequence_length(self) -> float:
        # Средняя длина чтений
        total, cnt = 0, 0
        with self.reader:
            for r in self.reader.read():
                total += len(r.sequence)
                cnt += 1
        return total / cnt if cnt else 0.0

    def get_sequence_length_distribution(self) -> Dict[int, int]:
        # Распределение длин
        dist = defaultdict(int)
        with self.reader:
            for r in self.reader.read():
                dist[len(r.sequence)] += 1
        return dict(dist)

    def calculate_per_base_quality(self):
        # Сбор Phred баллов по позициям
        pos_q = defaultdict(list)
        with self.reader:
            for r in self.reader.read():
                scores = FastqReader._quality_string_to_scores(r.quality)
                for i, s in enumerate(scores):
                    pos_q[i].append(s)
        max_i = max(pos_q.keys()) if pos_q else -1
        mean_q, all_q = [], []
        for i in range(max_i + 1):
            qs = pos_q[i]
            mean_q.append(np.mean(qs) if qs else 0)
            all_q.append(qs)
        return mean_q, all_q

    def calculate_per_base_content(self) -> Dict[str, List[float]]:
        # Процент ACGT по позициям
        pos_c = defaultdict(lambda: {'A':0,'C':0,'G':0,'T':0,'total':0})
        with self.reader:
            for r in self.reader.read():
                for i, b in enumerate(r.sequence.upper()):
                    if b in 'ACGT':
                        pos_c[i][b] += 1
                        pos_c[i]['total'] += 1
        max_i = max(pos_c.keys()) if pos_c else -1
        res = {b: [] for b in 'ACGT'}
        for i in range(max_i + 1):
            cnt = pos_c[i]
            tot = cnt['total']
            for b in 'ACGT':
                res[b].append((cnt[b] / tot)*100 if tot else 0)
        return res

    # 4. ГРАФИКИ

    def plot_per_base_quality_matplotlib(self):
        # 4.1 Per Base Sequence Quality (matplotlib)
        mean_q, all_q = self.calculate_per_base_quality()
        if not mean_q:
            print("Нет данных")
            return
        pos = list(range(1, len(mean_q)+1))
        perc = []
        for qs in all_q:
            if qs:
                perc.append({
                    'q25': np.percentile(qs,25),
                    'med': np.median(qs),
                    'q75': np.percentile(qs,75),
                    'q10': np.percentile(qs,10),
                    'q90': np.percentile(qs,90),
                    'mean': np.mean(qs)})
            else:
                perc.append({k:0 for k in ['q25','med','q75','q10','q90','mean']})
        fig, ax = plt.subplots(figsize=(14,6))
        ax.axhspan(28,42,facecolor='green',alpha=0.1)
        ax.axhspan(20,28,facecolor='orange',alpha=0.1)
        ax.axhspan(0,20,facecolor='red',alpha=0.1)
        w = 0.8
        for i, p in zip(pos, perc):
            ax.add_patch(plt.Rectangle((i-w/2,p['q25']),w,p['q75']-p['q25'],facecolor='yellow',edgecolor='black'))
            ax.plot([i-w/2,i+w/2],[p['med'],p['med']],'r-')
            ax.plot([i,i],[p['q75'],p['q90']],'k-')
            ax.plot([i,i],[p['q25'],p['q10']],'k-')
        ax.plot(pos,[p['mean'] for p in perc],'b-',label='Mean')
        ax.set_xlabel('Position (bp)'); ax.set_ylabel('Phred score'); ax.set_ylim(0,42)
        ax.set_title('Per Base Sequence Quality'); ax.legend(); plt.tight_layout()
        path = os.path.join(self.graphs_dir,'per_base_quality_matplotlib.png')
        plt.savefig(path,dpi=300); plt.close(); print(f"Saved: {path}")

    def plot_per_base_quality_seaborn(self):
        # 4.2 Per Base Sequence Quality (seaborn)
        mean_q, _ = self.calculate_per_base_quality()
        if not mean_q:
            print("Нет данных"); return
        pos = list(range(1,len(mean_q)+1))
        sns.set_style('whitegrid')
        fig, ax = plt.subplots(figsize=(14,6))
        ax.axhspan(28,42,facecolor='green',alpha=0.1)
        ax.axhspan(20,28,facecolor='orange',alpha=0.1)
        ax.axhspan(0,20,facecolor='red',alpha=0.1)
        ax.plot(pos, mean_q, color='blue',linewidth=2, label='Mean')
        ax.set_xlabel('Position (bp)'); ax.set_ylabel('Phred score')
        ax.set_title('Per Base Sequence Quality (Seaborn)'); ax.set_ylim(0,42)
        ax.legend(); plt.tight_layout()
        path = os.path.join(self.graphs_dir,'per_base_quality_seaborn.png')
        plt.savefig(path,dpi=300); plt.close(); print(f"Saved: {path}")

    def plot_per_base_quality_plotly(self):
        # 4.3 Per Base Sequence Quality (interactive Plotly)
        mean_q, all_q = self.calculate_per_base_quality()
        if not mean_q:
            print("Нет данных"); return
        pos = list(range(1,len(mean_q)+1))
        med = [np.median(q) if q else 0 for q in all_q]
        q25 = [np.percentile(q,25) if q else 0 for q in all_q]
        q75 = [np.percentile(q,75) if q else 0 for q in all_q]
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=pos,y=mean_q,mode='lines',name='Mean'))
        fig.add_trace(go.Scatter(x=pos,y=med,mode='lines',name='Median'))
        fig.add_trace(go.Scatter(x=pos+pos[::-1],y=q75+q25[::-1],fill='toself',fillcolor='rgba(255,255,0,0.3)',name='IQR'))
        fig.add_hrect(y0=28,y1=42,fillcolor='green',opacity=0.1)
        fig.add_hrect(y0=20,y1=28,fillcolor='orange',opacity=0.1)
        fig.add_hrect(y0=0,y1=20,fillcolor='red',opacity=0.1)
        fig.update_layout(title='Per Base Sequence Quality (Interactive)',xaxis_title='Position (bp)',yaxis_title='Phred score',yaxis=dict(range=[0,42]))
        path = os.path.join(self.graphs_dir,'per_base_quality_plotly.html')
        fig.write_html(path); print(f"Saved: {path}")

    def plot_per_base_content(self):
        # 4.4 Per Base Sequence Content
        content = self.calculate_per_base_content()
        if not content['A']:
            print("Нет данных"); return
        pos = list(range(1,len(content['A'])+1))
        fig, ax = plt.subplots(figsize=(14,6))
        ax.plot(pos, content['A'], label='%A', color='green')
        ax.plot(pos, content['C'], label='%C', color='blue')
        ax.plot(pos, content['G'], label='%G', color='black')
        ax.plot(pos, content['T'], label='%T', color='red')
        ax.set_xlabel('Position (bp)'); ax.set_ylabel('Percent (%)'); ax.set_ylim(0,100)
        ax.set_title('Per Base Sequence Content'); ax.legend(); ax.grid(alpha=0.3)
        plt.tight_layout()
        path = os.path.join(self.graphs_dir,'per_base_content.png')
        plt.savefig(path,dpi=300); plt.close(); print(f"Saved: {path}")

    def plot_sequence_length_distribution(self):
        # 4.5 Sequence Length Distribution
        dist = self.get_sequence_length_distribution()
        if not dist:
            print("Нет данных"); return
        L = sorted(dist.keys()); C = [dist[i] for i in L]
        fig, ax = plt.subplots(figsize=(12,6))
        ax.bar(L, C, color='steelblue',edgecolor='black',alpha=0.7)
        ax.set_xlabel('Length (bp)'); ax.set_ylabel('Count')
        ax.set_title('Sequence Length Distribution'); ax.grid(axis='y',alpha=0.3)
        plt.tight_layout()
        path = os.path.join(self.graphs_dir,'sequence_length_distribution.png')
        plt.savefig(path,dpi=300); plt.close(); print(f"Saved: {path}")

# 5. ЗАПУСК

if __name__ == "__main__":
    print("Демоверсия анализа FASTQ файлов")
    FASTQ_FILE = "demo/test_data_SRR35776190.fastq"
    if not os.path.isfile(FASTQ_FILE):
        print(f"ERROR: not found {FASTQ_FILE}")
        exit(1)
    graphs_dir = os.path.join(os.path.dirname(os.path.abspath(FASTQ_FILE)), "graphs")
    if not os.path.exists(graphs_dir):
        os.makedirs(graphs_dir)
        print(f"Created graphs dir: {graphs_dir}")
    analyzer = FastqAnalyzer(FASTQ_FILE, graphs_dir)
    print(f"Count: {analyzer.get_sequence_count()}")
    print(f"Avg length: {analyzer.get_average_sequence_length():.2f} bp")
    analyzer.plot_per_base_quality_matplotlib()
    analyzer.plot_per_base_quality_seaborn()
    analyzer.plot_per_base_quality_plotly()
    analyzer.plot_per_base_content()
    analyzer.plot_sequence_length_distribution()
    print("Demo complete")