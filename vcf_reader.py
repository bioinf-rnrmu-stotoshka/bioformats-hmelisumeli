# Получение заголовка и информации по отдельным группам заголовков
# Получение количества вариантов.
# Получение статистики “количество выравниваний - регион.” (Используйте pandas)
# Получение вариантов, лежащем в определенном геномном отрезке (аналог bedtools intersect).
import gzip

class Vcf_reader:
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
            
vcf = Vcf_reader(r"C:\Users\veram\.conda\OneDrive\Documents\labs\lab 14\homo_sapiens-chrMT.vcf\homo_sapiens-chrMT.vcf")
print(vcf.contig())
print(vcf.info())     