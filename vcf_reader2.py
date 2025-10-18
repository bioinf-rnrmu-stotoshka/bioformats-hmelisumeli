
import pandas as pd #для работы с данными в таблице
class Vcf_reader:
  # Получение заголовка и информации по отдельным группам заголовков
    def __init__(self, path:str):
        self.path = path

    def title(self):
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##')]
            
    def info(self):
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##INFO')]
            
    def filter(self):
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##FILTER')]
            
    def format(self):  #oпределяет, какие параметры могут быть указаны для каждого образца
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##FORMAT')]
    
    def alt(self): #описание альтернативных типов аллелей
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##ALT')]
            
    def contig(self): #информация о хромосомах
        with open(self.path, 'rt', encoding='latin-1') as f:
            return [x.strip() for x in f if x.startswith('##contig')]

    # Получение количества вариантов.    
    def count(self):
        count = 0
        with open(self.path, 'rt', encoding='latin-1') as f:
            for x in f:
                if not x.startswith('#'):
                    count += 1
        return count

# Получение статистики “количество выравниваний - регион.” (Используйте pandas)
    def stats(self, region_size=1000):
        d = pd.read_csv(self.path, sep='\t', comment='#', header=None) #header is none чтобы прочитать первую строку именно как данные а не заголовок
        d.columns = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'] #задаём эти заголовки сами
        d['DP'] = d['INFO'].str.extract(r'DP=(\d+)') #реджекс для создания отдельного столбца с числами если в инфо есть дп
        d['DP'] = d['DP'].fillna(0).astype(int) #fillna заменяет все пустые значения на 0. astype переводит значения в целые числа
        result = {}
        for index, row in d.iterrows(): #проходим по каждой строке
            region_start = (row['POS'] // region_size) * region_size #опр начало региона
            key = (row['CHROM'], region_start)
            if key not in result:
                result[key] = {'TOTAL_DEPTH': 0, 'VARIANT_COUNT': 0}
            result[key]['TOTAL_DEPTH'] += row['DP']
            result[key]['VARIANT_COUNT'] += 1
        #преобразуем словарь в таблицу
        data = []
        for k, v in result.items():
            data.append([k[0], k[1], v['TOTAL_DEPTH'], v['VARIANT_COUNT']])
        return pd.DataFrame(data, columns=['CHROM','REGION','TOTAL_DEPTH','VARIANT_COUNT'])

# Получение вариантов, лежащем в определенном геномном отрезке (аналог bedtools intersect).
    def varregion(self, chrom, start, end):
        variants_in_region = []
        with open(self.path, 'rt', encoding='latin-1') as f:  # лэтин для чтения спец символов в всф
            for x in f:  # проходим по каждой строке файла
                if x.startswith('#'):
                    continue # пропускаем заголовки
                columns = x.strip().split('\t')  #убирает пробелы и символы переноса строки в начале и конце + разделяет строку на столбцы 
                x_chrom = columns[0] 
                x_pos = int(columns[1])  #позицию переводим в число
                if x_chrom == chrom and start <= x_pos <= end: 
                    variants_in_region.append(columns)
        return variants_in_region
