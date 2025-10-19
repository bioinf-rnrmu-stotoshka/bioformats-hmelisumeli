class FastaAnalyzer:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
    
    def fasta_sequence_generator(self):
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


#использование
#if __name__ == "__main__":
    #analyzer = FastaAnalyzer("example.fasta")
    
    #статистика
    #average_length, sequence_count = analyzer.fasta_counter()
    
    #print(f"Количество последовательностей: {sequence_count}")
    #print(f"Средняя длина последовательности: {average_length:.2f}")

