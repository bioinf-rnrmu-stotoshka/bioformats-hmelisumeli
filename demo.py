from fastacounterF import FastaAnalyzer
from fastq_reader import FastqAnalyzer
from sam_reader import SAMfile
from vcf_reader3 import Vcf_reader
import os


# === FASTA ===
def demo_fasta():
    print("\n=== Пример работы с FASTA ===")
    fasta_path = "example.fasta"

    # Создаём тестовый файл, если его нет
    if not os.path.exists(fasta_path):
        with open(fasta_path, "w") as f:
            f.write(">seq1\nATGCGTAGCTAG\n>seq2\nATTTGGGCCCAA")

    analyzer = FastaAnalyzer(fasta_path)
    avg_len, seq_count = analyzer.fasta_counter()
    print(f"Количество последовательностей: {seq_count}")
    print(f"Средняя длина последовательностей: {avg_len:.2f}")


# === FASTQ ===
def demo_fastq():
    print("\n=== Пример работы с FASTQ ===")
    fastq_path = "example.fastq"

    # Минимальный FASTQ для примера
    if not os.path.exists(fastq_path):
        with open(fastq_path, "w") as f:
            f.write(
                "@read1\nATCGT\n+\nIIIII\n"
                "@read2\nGGGTT\n+\nIIHHH\n"
            )

    graphs_dir = "graphs"
    os.makedirs(graphs_dir, exist_ok=True)

    analyzer = FastqAnalyzer(fastq_path, graphs_dir)
    print(f"Количество последовательностей: {analyzer.get_sequence_count()}")
    print(f"Средняя длина: {analyzer.get_average_sequence_length():.2f}")
    analyzer.plot_sequence_length_distribution()
    print("График сохранён в папке graphs/")


# === SAM ===
def demo_sam():
    print("\n=== Пример работы с SAM ===")
    sam_path = "example.sam"

    if not os.path.exists(sam_path):
        with open(sam_path, "w") as f:
            f.write(
                "@HD\tVN:1.0\tSO:unsorted\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "r001\t99\tchr1\t7\t30\t8M2I4M1D3M\t=\t37\t39\tACCTTGAACTGAC\t*\n"
            )

    headers, alignments = SAMfile.parse_sam_file(sam_path)
    print(f"Заголовков: {len(headers)}")
    print(f"Выравниваний: {SAMfile.count_alignments(alignments)}")

    stats = SAMfile.get_chromosome_stats(alignments)
    print("\nСтатистика по хромосомам:")
    print(stats)


# === VCF ===
def demo_vcf():
    print("\n=== Пример работы с VCF ===")
    vcf_path = "example.vcf"

    if not os.path.exists(vcf_path):
        with open(vcf_path, "w") as f:
            f.write(
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=35\n"
                "chr1\t500\t.\tG\tC\t99\tPASS\tDP=60\n"
            )

    reader = Vcf_reader(vcf_path)
    print(f"Количество вариантов: {reader.count()}")

    stats = reader.stats()
    print("\nПример статистики:")
    print(stats.head())


# === Главный запуск ===
if __name__ == "__main__":
    print("=== Демонстрация биоинформатических форматов ===")
    demo_fasta()
    demo_fastq()
    demo_sam()
    demo_vcf()
    print("\n Демонстрация завершена успешно!")
