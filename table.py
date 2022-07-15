
from collections import deque


class Tabla:
    def __init__(self, path,context=True):
        self.path = path
        self.T = []
        self.V = []
        self.NC = []
        self.num_cepas = 0
        self.cepas_names = []
        self.posiciones = []
        self.context = context

        self.cargar()

    def cargar(self):
        with open(self.path, "r") as f:
            for line in f:
                l = line.split("\t")
                self.T.append(l)
            f.close()
        registro = []
        n = 0
        for fila in self.T:
            if fila[-1] not in registro:
                n += 1
                registro.append(fila[-1])
            fila.append(n)
        self.cepas_names = registro


        # Pasar a int la primera y última columna, ordenamiento y etiqueta:

        for linea in self.T:
            linea[0] = int(linea[0])
            # linea[14] = int(linea[14][:-1])
            if linea[1] != '.' and linea[2] != '.': # cambio de base
                linea.append(linea[2])
            elif linea[1] == '.': # inserción
                linea.append("I")
            else: # deleción
                linea.append("D")
        
        self.num_cepas = self.T[-1][15] # se declara el numero de cepas

        self.T.sort(key=lambda x: x[0], reverse=False)
        return
    

    def dibujar_entorno_NC(self, entorno):
        pos = []
        cepas = []
        inserts = self.check_inserts(entorno)
        Ref = []
        if inserts == ([0],[0]):
            # Crear Ref
            for snp in entorno:
                if snp[0] not in pos:
                    Ref.append(snp[8][10])
                    pos.append(snp[0])
        else:
            for snp in entorno:
                if snp[-1]=='I' and snp[0] in inserts[0]: # si es una inserción y está en inserts (no se ha usado)
                    pos_ins = inserts[0].popleft()
                    cant_ins = inserts[1].popleft()
                    for i in range(cant_ins):
                        pos.append(pos_ins)
                        Ref.append('-')
                elif snp[-1]!='I' and snp[0] not in pos: # si no es una inserción y es una posicion nueva
                    Ref.append(snp[8][10])
                    pos.append(snp[0])
        
        # Crear cepas
        snps_cepas = []
        for i in range(self.num_cepas):
            snps_cepas.append([[],[]])
        for fila in entorno:
            snps_cepas[fila[15]-1][0].append(fila[0]) # pos
            if fila[2] == '.':
                snps_cepas[fila[15]-1][1].append('-') # reemplazo
            else:
                snps_cepas[fila[15]-1][1].append(fila[2]) # reemplazo
        for i in range(self.num_cepas):
            cepas.append([])
            for coord, nuc in zip(pos,Ref):
                if coord not in snps_cepas[i][0]:
                    cepas[i].append(nuc)
                else:
                    indx = snps_cepas[i][0].index(coord)
                    cepas[i].append(snps_cepas[i][1][indx])
        return pos, Ref, cepas

    def crear_entorno(self, origen=0):
        """
        Devuelve el entorno de la posición indicada
        """
        rango = 20
        entorno = []
        entorno.append(self.T[origen]) # se une la fila inicial
        pos = self.T[origen][0]

        while True:
            if (origen + 1) >= len(self.T): # Final
                break
            pos = self.T[origen + 1][0] # siguiente posición
            if pos > self.T[origen][0] + rango: # si la sgte pos está fuera de rango, se acaba el ciclo
                break
            else: # si no está fuera de rango, se agrega al entorno, y se declara como nuevo origen
                entorno.append(self.T[origen + 1])
                origen += 1
        sgte_origen = origen + 1
        return entorno, sgte_origen
    
    def armar_visualizacion(self):
        start = 0
        end = len(self.T) # Completo
        # end = 80 # Pequeño (mejor display)
        for i in range(self.num_cepas + 1): # SIGUE DESDE AQUÍ
            self.V.append([])
        while start < end:
            entorno_actual, sgte = self.crear_entorno(origen=start)
            if self.context:
                P, R, C = self.dibujar_entorno(entorno_actual)
            else:
                P, R, C = self.dibujar_entorno_NC(entorno_actual)
            self.V[0] += R
            self.posiciones += P
            for i in range(self.num_cepas):
                self.V[i + 1] += C[i]
            # Quitamos los asteriscos
            # for i in range(self.num_cepas + 2):
            #     self.V[i] += ['*','*','*','*']
            start = sgte
        self.V.insert(0, ['>Ref'])
        counter = 2
        for name in self.cepas_names:
            clean_name = ['>' + name[:-10]]
            self.V.insert(counter, clean_name)
            counter += 2
        return


    def dibujar_entorno(self, entorno):
        inicio = entorno[0]
        pos = []
        cepas = []
        inserts = self.check_inserts(entorno)
        if inserts == ([0],[0]):
            # Se arma Ref y pos
            Ref = list(inicio[8])
            for i in range(inicio[0] - 10, inicio[0] + 11):
                pos.append(i)
            start = inicio[0]
            for snps in entorno:
                dist = snps[0] - start
                if dist > 0:
                    Ref += list(snps[8][-dist:])
                    last = pos[-1]
                    for i in range(dist):
                        pos.append(last + dist + 1)
                    start = snps[0]
            # Se arma cepas
            snps_cepas = []
            for i in range(self.num_cepas):
                snps_cepas.append([[],[]])
            for fila in entorno:
                snps_cepas[fila[15]-1][0].append(fila[0]) # pos
                if fila[2] == '.':
                    snps_cepas[fila[15]-1][1].append('-') # reemplazo
                else:
                    snps_cepas[fila[15]-1][1].append(fila[2]) # reemplazo
            for i in range(self.num_cepas):
                cepas.append([])
                for coord, nuc in zip(pos,Ref):
                    if coord not in snps_cepas[i][0]:
                        cepas[i].append(nuc)
                    else:
                        indx = snps_cepas[i][0].index(coord)
                        cepas[i].append(snps_cepas[i][1][indx])

        else:
            # Voy a asumir que solo hay un indel cada 21 nuc (sino no sería snp)
            # Se arma Ref y pos
            Ref = list(inicio[8][:10])
            for i in range(inicio[0] - 10, inserts[0][0]):
                pos.append(i)
            while len(inserts[0]) > 0:
                pos_ins = inserts[0].popleft()
                cant_ins = inserts[1].popleft()
                for i in range(cant_ins):
                    pos.append(pos_ins)
                    Ref.append('-')
                if len(inserts[0]) > 0:
                    for i in range(pos_ins + 1, inserts[0][0]):
                        pos.append(i)
            for i in range(pos_ins + 1, pos_ins + 11):
                pos.append(i)
            Ref += list(inicio[8][-10:])
            # Se arma cepas
            snps_cepas = []
            for i in range(self.num_cepas):
                snps_cepas.append([[],[]])
            for fila in entorno:
                snps_cepas[fila[15]-1][0].append(fila[0]) # pos
                snps_cepas[fila[15]-1][1].append(fila[2]) # reemplazo
            for i in range(self.num_cepas):
                cepas.append([])
                for coord, nuc in zip(pos,Ref):
                    if coord not in snps_cepas[i][0]:
                        cepas[i].append(nuc)
                    else:
                        indx = snps_cepas[i][0].index(coord)
                        cepas[i].append(snps_cepas[i][1][indx])
                        snps_cepas[i][0].pop(indx)
                        snps_cepas[i][1].pop(indx)
        return pos, Ref, cepas


    def check_inserts(self, entorno):
        """
        Retorna (deque) la posición(es) y largo de las inserciones en un entorno dado.
        """
        posiciones = dict()
        largo_max = deque()
        for snp in entorno:
            if snp[-1] == 'I' and snp[0] not in posiciones:
                contadores = []
                for i in range(self.num_cepas):
                    contadores.append(0)
                contadores[snp[-2] -1] = 1
                posiciones[snp[0]] = contadores
            elif snp[-1] == 'I':
                posiciones[snp[0]][snp[-2]-1] =  posiciones[snp[0]][snp[-2]-1] + 1
        f_pos = deque()
        if bool(posiciones) == False:
            return [0], [0]
        else:
            for p in posiciones:
                f_pos.append(p)
                largo_max.append(max(posiciones[p]))
            return f_pos, largo_max
