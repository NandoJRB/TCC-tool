'''
    Nesse módulo irei usar plugins para o software Primer3 que será usado para
    projetar primers. Por enquanto, estou só usando para fazer testes mesmo. 
'''

import primer3 as _p3



a = _p3.designPrimers(seq_args={
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'
                             'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'
                             'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'
                             'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
                             'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'
                             'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'
                             'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'
                             'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'
                             'GATGTTTCCCTCTAGTAG',
        'SEQUENCE_INCLUDED_REGION': [36,342]
    },
    global_args={
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 5,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [
            [75,100],[100,125],[125,150],
            [150,175],[175,200],[200,225]
        ],
    })


# A função calcHeterodimer recebe duas sequências distintas e informa se existe dímero entre elas
b = _p3.calcHeterodimer('GCTTGCATGCCTGCAGGTCGACTCTCAGAGGATCCCCCTACATTTT', 'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG')

# A função calcHairpin verifica se a formação de um anel no primer
c = _p3.calcHairpin('GATGTTTCCCTCTAGTAG')

# A função verifica se ocorre homodímero, isto é, primer parea com outro primer igual
d = _p3.calcHomodimer('TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA')

# Calcula a temperatura de Melting
e = _p3.calcTm('TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA')
print(e)

