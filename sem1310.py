#Требования к fasta
# Получение количества последовательностей
# Получение средней длины последовательностей

def fasta_sequence_generator(fasta_file):
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
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
        
        if current_header is not None and current_sequence:
            full_sequence = ''.join(current_sequence)
            yield current_header, full_sequence

def fasta_counter(fasta_file):
    sequences = {}
    current_header = None
    current_sequence = []
    total_length = 0 
    sequence_count = 0  #n
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:  # пропускаем пустые строки
                continue
                
            if line.startswith('>'):
                if current_header is not None:
                    full_sequence = ''.join(current_sequence)
                    sequences[current_header] = full_sequence
                    total_length += len(full_sequence)  
                    sequence_count += 1
                
                current_header = line[1:]  #убираем '>'
                current_sequence = []
            else:
                current_sequence.append(line)
        
        #обрабатываем последнюю последовательность отдельно, потому что нет тригера >
        if current_header is not None:
            full_sequence = ''.join(current_sequence)
            sequences[current_header] = full_sequence
            total_length += len(full_sequence)
            sequence_count += 1
    
    #вычисляем среднюю длину если есть последовательности
    if sequence_count > 0:
        average_length = total_length / sequence_count
    else:
        average_length = 0
    
    return average_length, sequence_count