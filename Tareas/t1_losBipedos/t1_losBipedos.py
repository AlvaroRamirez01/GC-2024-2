import numpy as np
import math, random, re

# Solucion Ejercicio 1.2
lista_secuencias = [
    'ATATATACATACTGGTAATGGGCGCGCGTGTGTTAAGTTCTGTTGTAGGGGTGATTAGGGGCG',
    'GGCCCACACCCCACACCAATATATGTGGTGTGGGCTCCACTCTCTCGCGCTCGCGCTGGGGAT',
    'ATAAGGTGTGTGGGCGCGCCCCGCGCGCGCGTTTTTTCGCGCGCCCCCGCGCGCGCGCGCGCG',
    'GGCGCGGGACGCGGCGGCGGATCCCGATCCGTGCGTCAATACTATTATGGCCAGATAGAATAA',
    'GTGCTGCTGCGGCGCCCACACCTATTATCTCTCTCTCTCTGCCTCTCCACCTCGGGGCTTAAT',
    'GCGCTGCTGCTGGCTCGATGGGCGCGTGCGTCGTAGCTCGATGCTGGCTCGAGCTGTAATCTT',
    'GGCGCTCGCTCGGATGCGCGGCCGGGCTCTCTGCTCGCGCTCGCTTCGCGCTCGTGACCGCTG',
    'AATTGGTGCGCGCTCGCGCACACACAGAGAGAGGGTTTATATAGGATGATATATCCACATTGG',
    'ATGCTGCTGCTGGCTCTGCTTGCGCTCTGCTCGCTGGGGTGTGTGTGCCGCGCGCTGCTGCTC',
    'GCTGGGCTCGCTCGATGCGCGCGGGCGCGCGACCGCGGACGGCGTCGCTGCTAAATGGGCTTC']

def ejercicio_1_2(lista_secuencias):
    L = []
    r = re.compile('(ATG|TTG|GTG)([ACTG]{3})+(TAG)')
    L = [i for i, item in enumerate(lista_secuencias) if re.search(r, item)]
    return L

# Solucion Ejercicio 2.1
promotores = ['AGATAG', 'TGATAG', 'AGATAA', 'TGATAA']

def leer_archivo(archivo):
    with open(archivo, 'r') as file:
        for line in file:
            yield line.strip()
            
def ejercicio_2_1(archivo, promotores):
    L = []
    for line in leer_archivo(archivo):
        contador = 0
        identificador = re.split(r'\s+', line)[0]
        line = re.split(r'\s+', line)[1]
        for _ in re.finditer('|'.join(promotores), line):
            contador += 1
        L.append(identificador+" : " + str(contador))
    return L

# lista = ejercicio_2_1("promotores.txt", promotores)
# for item in lista:
#     print(item)

# Solucion Ejercicio 2.2
def calcula_pi(iteraciones: int) -> float:
    numeros_dentro = 1
    for _ in range(iteraciones):
        x = np.random.uniform(-1, 1)
        y = np.random.uniform(-1, 1)
        if np.sqrt((x**2) + (y**2)) <= 1:
            numeros_dentro += 1
    return 4 * (numeros_dentro / iteraciones)

# Solucion Ejercicio 4

# Importamos las bibliotecas necesarias
from Bio import SeqIO
import re, math

# cargando genoma de SARS-CoV2
genome = SeqIO.read("sequence.fasta", "fasta")
sequence = str(genome.seq)
r_sequence = str(genome.seq.reverse_complement())

# Numero de ORFs
patron = r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)'
orfs = re.findall(patron, sequence)
print(str(len(orfs)))

# ORFs con probabilidad < 0.05 de ser no espurios
k = math.log(0.05, 61/64)
orfs_5 = []
for orf in orfs:
    if len(orf) > int(k):
        orfs_5.append(orf)
print(str(len(orfs_5)))

# ORFs con probabilidad < 0.01 de no ser espurios
k = math.log(0.01, 61/64)
orfs_1 = []
for orf in orfs:
    if len(orf) > int(k):
        orfs_1.append(orf)
print(str(len(orfs_1)))