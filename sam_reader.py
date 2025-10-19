import pandas as pd


class SAMfile:
    def parse_sam_file(filename):  # Читаем файл и разделяем на заголовки и выравнивания
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

    def get_header_groups(
        headers,
    ):  # Группируем заголовки по типам (@HD, @SQ, @RG, @PG)
        groups = {}

        for header in headers:
            header_type = header.split("\t")[0]
            if header_type not in groups:
                groups[header_type] = []
            groups[header_type].append(header)

        return groups

    def print_headers(headers):  # Выводит информацию о заголовках
        groups = get_header_groups(headers)

        print("Заголовки sam-файла:\n")
        for header_type, lines in groups.items():
            print(f"{header_type} - найдено {len(lines)} записей:")
            for line in lines:
                print(f"  {line}")
            print()

    def count_alignments(alignments):  # Считаем кол-во выравниваний
        return len(alignments)

    def get_chromosome_stats(alignments):  # Делаем статистику по хромосомам
        chromosomes = []

        for alignment in alignments:
            fields = alignment.split("\t")
            chromosome = fields[2]
            chromosomes.append(chromosome)

        df = pd.DataFrame({"Хромосома": chromosomes})
        stats = df["Хромосома"].value_counts().reset_index()
        stats.columns = ["Хромосома", "Количество"]

        return stats

    def calculate_alignment_end(
        position, cigar
    ):  # Вычисляет конечную позицию выравнивания из CIGAR строки
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

    def find_alignments_in_region(
        alignments, chromosome, start, end
    ):  # Находим выравнивания в заданном геномном регионе
        results = []

        for alignment in alignments:
            fields = alignment.split("\t")

            read_name = fields[0]
            chrom = fields[2]
            position = int(fields[3])
            cigar = fields[5]

            if chrom != chromosome:
                continue

            align_end = calculate_alignment_end(position, cigar)

            # Проверяем пересечение с регионом
            if not (align_end < start or position > end):
                results.append(
                    {
                        "Название": read_name,
                        "Хромосома": chrom,
                        "Начало": position,
                        "Конец": align_end,
                        "CIGAR": cigar,
                    }
                )

        return pd.DataFrame(results)
