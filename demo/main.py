# Программа-объединяющая скрипты для анализа FASTA, FASTQ, SAM, VCF файлов.
# Запрашивает путь к файлу, определяет тип и вызывает анализатор из папки-пакета с скриптами.


import os
import sys


def detect_file_type(filepath):
    ext = os.path.splitext(filepath.lower())[1]
    if filepath.lower().endswith('.gz'):
        base = os.path.splitext(filepath[:-3].lower())[1]
        if base in ['.fastq', '.fq']:
            return 'fastq'
        elif base in ['.fasta', '.fa']:
            return 'fasta'
    if ext in ['.fq', '.fastq']:
        return 'fastq'
    elif ext in ['.fa', '.fasta', '.fna']:
        return 'fasta'
    elif ext == '.sam':
        return 'sam'
    elif ext == '.vcf':
        return 'vcf'
    return None


def analyze_fasta(filepath):
    from code.fastacounterF import FastaAnalyzer
    print("\nАнализ FASTA")
    analyzer = FastaAnalyzer(filepath)
    avg_len, count = analyzer.fasta_counter()
    print(f"\nКоличество последовательностей: {count}")
    print(f"Средняя длина последовательности: {avg_len:.2f} bp")


def analyze_fastq(filepath):
    from code.fastq_reader import FastqAnalyzer
    print("\nАнализ FASTQ")
    graphs_dir = os.path.join(os.path.dirname(os.path.abspath(filepath)), "graphs")
    if not os.path.exists(graphs_dir):
        os.makedirs(graphs_dir)
    analyzer = FastqAnalyzer(filepath, graphs_dir)
    print("\nВычисление статистики:")
    seq_count = analyzer.get_sequence_count()
    avg_len = analyzer.get_average_sequence_length()
    print(f"Количество последовательностей: {seq_count}")
    print(f"Средняя длина последовательности: {avg_len:.2f} bp")
    print("Генерация графиков:")
    analyzer.plot_per_base_quality_matplotlib()
    analyzer.plot_per_base_quality_seaborn()
    analyzer.plot_per_base_quality_plotly()
    analyzer.plot_per_base_content()
    analyzer.plot_sequence_length_distribution()
    print(f"Графики сохранены в '{graphs_dir}'")


def analyze_vcf(filepath):
    from code.vcf_reader3 import Vcf_reader

    print("\nАнализ VCF")
    vcf = Vcf_reader(filepath)

    # Заголовки по группам
    print("\nЗаголовки общего уровня (##):")
    for line in vcf.title():
        print(line)
    print("\nГруппа INFO (##INFO):")
    for line in vcf.info():
        print(line)
    print("\nГруппа FILTER (##FILTER):")
    for line in vcf.filter():
        print(line)
    print("\nГруппа FORMAT (##FORMAT):")
    for line in vcf.format():
        print(line)
    print("\nГруппа ALT (##ALT):")
    for line in vcf.alt():
        print(line)
    print("\nГруппа CONTIG (##contig):")
    for line in vcf.contig():
        print(line)

    # Количество вариантов
    print(f"\nКоличество вариантов: {vcf.count()}")

    # Статистика по регионам
    stats = vcf.stats(region_size=1000)
    print("\nСтатистика по регионам (первые 10):")
    print(stats.head(10).to_string(index=False))

    # Интерактивный ввод для получения вариантов в заданном отрезке
    print("\nПолучение вариантов в определённом геномном отрезке")
    chrom = input("Введите хромосому (значение из столбца CHROM): ").strip()
    start = int(input("Введите начальную позицию (в пределах значений из столбца REGION)): ").strip())
    end = int(input("Введите конечную позицию (в пределах значений из столбца RREGION): ").strip())
    var_region = vcf.varregion(chrom, start, end)
    if var_region:
        print(f"\nНайдено {len(var_region)} вариантов в {chrom}:{start}-{end}")
        for i, var in enumerate(var_region[:10], 1):
            print(f"{i}. CHROM={var[0]}, POS={var[1]}, REF={var[3]}, ALT={var[4]}")
    else:
        print("\nВарианты в указанном регионе не найдены.")

def analyze_sam(filepath):
    from code.sam_reader import SAMfile

    print("\nАнализ SAM")
    sam = SAMfile()

    # Заголовки по группам
    headers, alignments = sam.parse_sam_file(filepath)
    groups = sam.get_header_groups(headers)
    print("\nЗаголовки по группам:")
    for group, lines in groups.items():
        print(f"{group} ({len(lines)}):")
        for line in lines:
            print(f"  {line}")

    # Общее количество выравниваний
    total = sam.count_alignments(alignments)
    print(f"\nОбщее количество выравниваний: {total}")

    # Расширенная статистика по хромосомам с мин/макс позициями
    stats = sam.get_chromosome_stats(alignments)
    print("\nСтатистика по хромосомам:")
    print(stats.to_string(index=False))

    # Интерактивный ввод для поиска выравниваний в регионе
    print("\nПоиск выравниваний в указанном геномном отрезке")
    chrom = input("Хромосома (значение из столбца Хромосома): ").strip()
    start = int(input("Начальная позиция (в пределах значений столбцов Мин_позиция и Макс_позиция): ").strip())
    end   = int(input("Конечная позиция (в пределах значений столбцов Мин_позиция и Макс_позиция): ").strip())
    results = sam.find_alignments_in_region(alignments, chrom, start, end)
    if not results.empty:
        print(f"\nНайдено {len(results)} выравниваний в {chrom}:{start}-{end}")
        print(results.to_string(index=False))
    else:
        print("\nВыравниваний не найдено.")


def main():
    print("Анализ файлов форматов: FASTA, FASTQ, SAM, VCF")
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = input("Введите путь к файлу: ").strip().strip('"').strip("'")


    if not os.path.isfile(filepath):
        print(f"Ошибка: файл не найден: {filepath}")
        sys.exit(1)


    ftype = detect_file_type(filepath)
    if ftype is None:
        print("Неизвестный формат файла. Поддерживаются .fasta, .fa, .fastq, .fq, .sam, .vcf (+ .gz)")
        sys.exit(1)


    print(f"Определён формат: {ftype.upper()}")
    if ftype == 'fasta':
        analyze_fasta(filepath)
    elif ftype == 'fastq':
        analyze_fastq(filepath)
    elif ftype == 'sam':
        analyze_sam(filepath)
    elif ftype == 'vcf':
        analyze_vcf(filepath)

    print("\nЕсли хотите проверить другой файл, запустите программу заново.")


if __name__ == "__main__":
    main()
