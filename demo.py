
from fastacounterF import FastaAnalyzer
from fastq_reader import FastqAnalyzer
from sam_reader import SAMfile
from vcf_reader3 import Vcf_reader


def demo_fasta():
    print("\n=== Пример работы с FASTA ===")
    reader = FastaAnalyzer("example.fasta")
    for seq_id, seq in reader.read_sequences().items():
        print(f"ID: {seq_id}, длина: {len(seq)}")
    reader.close()


def demo_fastq():
    print("\n=== Пример работы с FASTQ ===")
    reader = FastqAnalyzer("example.fastq")
    for seq_id, scores in reader.get_quality_scores().items():
        avg = reader.get_average_quality(seq_id)
        print(f"ID: {seq_id}, среднее качество: {avg:.2f}")
    reader.close()


def demo_sam():
    print("\n=== Пример работы с SAM ===")
    reader = SAMfile("example.sam")
    for aln in reader.read_alignments():
        print(aln)
        break  # показываем только первую запись
    reader.close()


def demo_vcf():
    print("\n=== Пример работы с VCF ===")
    reader = Vcf_reader("example.vcf")
    for var in reader.read_variants():
        print(var)
        break  # показываем только первую запись
    reader.close()


if __name__ == "__main__":
    demo_fasta()
    demo_fastq()
    demo_sam()
    demo_vcf()
