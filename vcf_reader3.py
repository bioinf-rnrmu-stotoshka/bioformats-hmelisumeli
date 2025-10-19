import pandas as pd #для работы с данными в таблице
class Vcf_reader:
  # Получение заголовка и информации по отдельным группам заголовков
    def __init__(self, path:str):
        self.path = path

    def _lines(self): #служебный метод для вн функций
        with open(self.path, 'rt', encoding='latin-1') as f:
            for x in f:
                yield x.strip()

    def title(self):
        return (x for x in self._lines() if x.startswith('##'))
            
    def info(self):
        return (x for x in self._lines() if x.startswith('##INFO'))
            
    def filter(self):
        return (x for x in self._lines() if x.startswith('##FILTER'))
            
    def format(self):  #oпределяет, какие параметры могут быть указаны для каждого образца
        return (x for x in self._lines() if x.startswith('##FORMAT'))
    
    def alt(self): #описание альтернативных типов аллелей
        return (x for x in self._lines() if x.startswith('##ALT'))
            
    def contig(self): #информация о хромосомах
        return (x for x in self._lines() if x.startswith('##contig'))

    # Получение количества вариантов.   
    def count(self):
        count = 0
        for x in self._lines():
            if not x.startswith('#'):
                count += 1
        return count

    # Получение статистики “количество выравниваний - регион.” (Используйте pandas)
    def stats(self, region_size=1000):
        result = {}
        def varlines(): # читает только строки с вариантами без заголовков
            for x in self._lines():
                if not x.startswith('#'):
                    yield x

        for y in varlines(): #проходим по каждой строке
            columns = y.split('\t') #разделяем по столбцам
            if len(columns) < 8: 
                continue
            chrom = columns[0]
            pos = int(columns[1])
            info = columns[7]
            dp = 0
            for x in info.split(';'): #перебор строки и разделение по ;
                if x.startswith('DP='):
                    dp = int(x[3:])
                    break

            region_start = (pos // region_size) * region_size #опр начало региона
            key = (chrom, region_start)
            if key not in result:
                result[key] = {'TOTAL_DEPTH': 0, 'VARIANT_COUNT': 0}
            result[key]['TOTAL_DEPTH'] += dp
            result[key]['VARIANT_COUNT'] += 1

        #преобразуем словарь в таблицу
        data = ([k[0], k[1], v['TOTAL_DEPTH'], v['VARIANT_COUNT']] for k, v in result.items())
        return pd.DataFrame(data, columns=['CHROM','REGION','TOTAL_DEPTH','VARIANT_COUNT'])

    # Получение вариантов, лежащем в определенном геномном отрезке (аналог bedtools intersect).
    def varregion(self, chrom, start, end):
        def gen(): #нужен чтобы вернуть не генератор
            for x in self._lines(): 
                if x.startswith('#'):
                    continue # пропускаем заголовки
                columns = x.strip().split('\t')  #убирает пробелы и символы переноса строки в начале и конце + разделяет строку на столбцы 
                if len(columns) < 2:  # пропуск пустых или битых строк
                    continue
                x_chrom = columns[0] 
                x_pos = int(columns[1])  #позицию переводим в число
                if x_chrom == chrom and start <= x_pos <= end: 
                    yield columns
        return list(gen())
