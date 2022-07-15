from table import Tabla

SNPS = Tabla("example_file.txt",context=False) # (EDITABLE)
SNPS.armar_visualizacion()

# Create output files

L = []
for linea in SNPS.V:
    l_s = "".join(str(e) for e in linea)
    L.append(l_s)

file1 = open("no_context_scheme.txt", "w") # SNPs file (EDITABLE)
for frase in L:
    file1.write(frase)
    file1.write('\n')
file1.close()

file2 = open("positions.txt", "w") # SNPs positions file (EDITABLE)
p = " ".join(str(e) for e in SNPS.posiciones)
file2.write(p)
file2.close()
