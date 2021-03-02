### Quantum Associative Memory by D. Ventura, T. Martinez
### Repository reference: https://gitlab.com/prince-ph0en1x/QaGs
#(by A. Sarkar)
## Importing libraries
#matplotlib inline
import qiskit
from qiskit import IBMQ
from qiskit import Aer
from qiskit import QuantumCircuit, ClassicalRegister,QuantumRegister, QiskitError
from qiskit.tools.visualization import circuit_drawer
from qiskit.tools.visualization import plot_histogram
from qiskit.providers.aer import noise
import random
import matplotlib
import math
from math import *

## Defining some auxiliary functions
def convertToNumEncoding(genome_str):
    bin_str = ""
    for i in range(0, len(genome_str)):
        if genome_str[i] == 'A':
            bin_str = bin_str + "0"
        elif genome_str[i] == 'C':
            bin_str = bin_str + "1"
        elif genome_str[i] == 'G':
            bin_str = bin_str + "2"
        elif genome_str[i] == 'T':
            bin_str = bin_str + "3"
    return bin_str

## Initializing global variables
# Alphabet set (0, 1, 2, 3) := {A, C, G, T} for DNA Nucleotide bases
AS = {'00', '01', '10', '11'}
A = len(AS)
genome_file = open("Sars_cov_2.txt", "r")
#genome_file = open("Carsonella_rudii.txt", "r")
#genome_file = open("nsp14.txt", "r")
RG = genome_file.read()
genome_file.close()
RG = convertToNumEncoding(RG)
RG = RG[0:160000]
N = len(RG)


# Short Read search string
SR = '23?'
M = len(SR)
Q_A = ceil(log2(A))
# Data Qubits
Q_D = Q_A * M
# Tag Qubits
Q_T = ceil(log2(N - M + 1))
# Ancilla Qubits
Q_anc = 8
# Ancilla qubits ids
anc = []
for qi in range(0, Q_anc):
    anc.append(Q_D + Q_T + qi)
# Total number of qubits
Q = Q_D + Q_T + Q_anc

#Start codon
Start = '032'
Xl = len(Start)

#Stop codon
Stop = ['302','320','300']

IBMQ.save_account('APIKEY', overwrite=True)
## Initialization of IBM QX
IBMQ.enable_account('bf333c31a0a3b9e5829fa0a0e0097216062f044a366ff7cce49dc0f004ad35da0d8067cb98c1f59e8d79704b055dac79de9e88c4e10306dc2cb14e20183f4e68')
#provider = IBMQ.get_provider()
# Pick an available backend
# If this isn’t available pick a backend whose name containes ’_qasm_simulator’ from the output above
#backend = provider.get_backend('ibmq_qasm_simulator')
backend = Aer.get_backend('qasm_simulator')

## Defining functions
def QAM():
    print('Reference genome:', RG)
    print('Chosen pattern for testing:', SR)
    print('Total number of qubits:', Q)
    print('Number of ancilla qubits:', Q_anc)
    qr = qiskit.QuantumRegister(Q)
    cr = qiskit.ClassicalRegister(Q_T)
    qc = qiskit.QuantumCircuit(qr, cr)
# Patterns are stored
    generateInitialState(qc, qr)
# Patterns are turned into Hamming Distances
    evolveToHammingDistances(qc, qr)
# Marking the zero Hamming Distance states
    markZeroHammingDistance(qc, qr)
    inversionAboutMean(qc, qr)
# Applying again this funcions turns back the Hamming Distances into the original patterns
    evolveToHammingDistances(qc, qr)
# Marking the stored patterns
    markStoredPatterns(qc, qr)
    # From patterns to Hamming Distances again
    evolveToHammingDistances(qc, qr)
    inversionAboutMean(qc, qr)
# Grover’s iterations
    it = 0
    for i in range(0, 1):
        markZeroHammingDistance(qc, qr)
        inversionAboutMean(qc, qr)
        it = it + 1
    print("Grover’s algorithm had {} iterations.".format(int(it)))
    finalGroverMeasurement(qc, qr, cr)
    return qc

# Initialization as suggested by L. C. L. Hollenberg
def generateInitialState(qc, qr):
    for qi in range(0, Q_T):
        qc.h(qr[qi])
    
    control_qubits = []
    for ci in range(0, Q_T):
        control_qubits.append(qr[ci])
    
    ancilla_qubits = []
    for qi in anc:
        ancilla_qubits.append(qr[qi])
    
    for qi in range(0, N - M + 1):
        qis = format(qi, '0' + str(Q_T) + 'b')
        for qisi in range(0, Q_T):
            if qis[qisi] == '0':
                qc.x(qr[qisi])
        wMi = RG[qi:qi + M]
        #print("Tag: {} - Data: {}".format(qis, wMi))

        for wisi in range(0, M):
            wisia = format(int(wMi[wisi]), '0' + str(Q_A) + 'b')
            for wisiai in range(0, Q_A):
                if wisia[wisiai] == '1':
                    qc.mct(control_qubits, qr[Q_T + wisi * Q_A + wisiai], ancilla_qubits)
        for qisi in range(0,Q_T):
            if qis[qisi] == '0':
                qc.x(qr[qisi])
                
## To compare given search string (GD)
        Seq = []
        counter = compare(wMi, Start)
        #print(counter, "yes")
        if counter == 3 :
            #print(Start, "found at ", qis)
            y = int(qis, 2)
            #print("y",y)
            for qi in range(y, N - M + 1, 3):
                qis = format(qi, '0' + str(Q_T) + 'b')
                for qisi in range(0, Q_T):
                    if qis[qisi] == '0':
                        qc.x(qr[qisi])
                wMi = RG[qi:qi + M]
                #print("Tag: {} - Data: {}".format(qis, wMi))
                
                Seq.append(wMi)
                Length=len(Seq)
                for i in range(0,3):
                    StopC = compare(wMi, Stop[i])
                #print(Length)
                if StopC == 3 and Length > 7096:
                    break
                
            print(Length, "\t codons long sequence found at \t", qis, "/",int(qis,2), ":\t",int(qis,2),"-",int(qis,2)+Length-1)
            #print(Seq)
    return

def compare(a,b):
    #print(a,b)
    c = 0
    for x, y in zip(a, b):
        if x == y:
            c+=1
            #print(x,y)
    return c


def evolveToHammingDistances(qc, qr):
    for pi in range(0, M):
        if SR[pi] == '?':
            continue
        ppi = format(int(SR[pi]), '0' + str(Q_A) + 'b')
        for ppii in range(0, Q_A):
            if ppi[ppii] == '1':
                qc.x(qr[Q_T + pi * Q_A + ppii])
    return

# Oracle to mark zero Hamming Distance
def markZeroHammingDistance(qc, qr):
    for qi in range(0, Q_D):
        qc.x(qr[Q_T + qi])
    control_qubits = []
    for mi in range(0, M):
        if SR[mi] == '?':
            for ai in range(0, Q_A):
                control_qubits.append(qr[Q_T + mi * Q_A + ai])
    
    ancilla_qubits = []
    for qi in anc:
        ancilla_qubits.append(qr[qi])
    target = control_qubits.pop()
    qc.h(target)
    qc.mct(control_qubits, target, ancilla_qubits)
    qc.h(target)
    for qi in range(0, Q_D):
        qc.x(qr[Q_T + qi])
    return

# Inversion about mean
def inversionAboutMean(qc, qr):
    for si in range(0, Q_D + Q_T):
        qc.h(qr[si])
        qc.x(qr[si])
    
    control_qubits = []
    for qi in range(1, Q_D + Q_T):
        control_qubits.append(qr[qi])
    
    ancilla_qubits = []
    for qi in anc:
        ancilla_qubits.append(qr[qi])
    qc.h(qr[0])
    qc.mct(control_qubits, qr[0], ancilla_qubits)
    qc.h(qr[0])
    for si in range(0, Q_D + Q_T):
        qc.x(qr[si])
        qc.h(qr[si])
    return

# Oracle to mark stored patterns
def markStoredPatterns(qc, qr):
    control_qubits = []
    for qi in range(1, Q_D + Q_T):
        control_qubits.append(qr[qi])
    
    ancilla_qubits = []
    for qi in anc:
        ancilla_qubits.append(qr[qi])
   
    for qi in range(0, N - M + 1):
        qis = format(qi, '0' + str(Q_T) + 'b')
        for qisi in range(0, Q_T):
            if qis[qisi] == '0':
                qc.x(qr[qisi])
        wMi = RG[qi:qi + M]
        for wisi in range(0, M):
            wisia = format(int(wMi[wisi]), '0' + str(Q_A) + 'b')
            for wisiai in range(0, Q_A):
                if wisia[wisiai] == '0':
                    qc.x(qr[Q_T + wisi * Q_A + wisiai])
        qc.h(qr[0])
        qc.mct(control_qubits, qr[0], ancilla_qubits)
        qc.h(qr[0])
        for wisi in range(0,M):
            wisia = format(int(wMi[wisi]), '0' + str(Q_A) + 'b')
            for wisiai in range(0, Q_A):
                if wisia[wisiai] == '0':
                    qc.x(qr[Q_T + wisi * Q_A + wisiai])
        for qisi in range(0, Q_T):
            if qis[qisi] == '0':
                qc.x(qr[qisi])
    return

# Final measurement
def finalGroverMeasurement(qc, qr, cr):
    for qi in range(0, Q_T):
        qc.measure(qr[qi], cr[qi])
    return

# Main function
if __name__ == '__main__':
# Printing some data for testing
    qc = QAM()
    print("Circuit depth: {}".format(qc.depth()))
# Total number of gates
    print("Number of gates: {}".format(len(qc.data)))
    gate_num = 1
    for item in qc.data:
        qb_list = ""
        for qb in item[1]:
            qb_list = qb_list + str(qb.index) + ", "
        qb_list = qb_list[:len(qb_list)-2]
        #print("#{}: {}, {}".format(gate_num, item[0].name, qb_list))
        gate_num = gate_num + 1

    
