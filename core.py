from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os as _stdos


_TERMINAL_SIZE = int(_stdos.get_terminal_size().columns * (2 / 3))
_SEQEND_SIZE = int(_TERMINAL_SIZE / 4)


'''
    A classe BioBlock é derivada de SeqRecord de Biopython. Essa classe tem
    como propósito representar uma sequência biologica qualquer que pode ser
    unida com outros BioBlocks em um Assembly
'''
class BioBlock(SeqRecord):

    def __init__(self, seq, 
                       id="<unknown id>", 
                       name="<unknown name>", 
                       description="<unknown description>", 
                       dbxrefs=None, 
                       features=None, 
                       annotations=None, 
                       letter_annotations=None,
                       forward_primer=None, 
                       reverse_primer=None, 
                       circular=False, 
                       leftbord=0, 
                       rightbord=-1):
        if isinstance(seq, SeqRecord):
            self.__dict__.update(seq.__dict__)
            self.seq = Seq(self.seq) 
        elif isinstance(seq, (str, Seq)):
            super().__init__(Seq(seq), 
                             id, 
                             name, 
                             description, 
                             dbxrefs, 
                             features, 
                             annotations, 
                             letter_annotations)
        else:
            raise TypeError('constructor might receive a SeqRecord object or sequence string')
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.circular = circular
        self.leftbord = leftbord
        self.rightbord = rightbord
        if self.rightbord >= len(self.seq) or self.rightbord < (-len(self.seq)):
            raise ValueError('invalid value for segment right border')
        elif self.rightbord < 0:
            self.rightbord += len(self.seq)

    def __str__(self):
        slen = len(self.seq)
        strform = f'BioBlock {self.id}\n\tsize: {slen}\n'
        if self.circular:
            strform += f'\ttopology: circular\n'
        else:
            strform += f'\ttopology: linear\n'
        if slen <= _TERMINAL_SIZE:
            strform += f'\tseq: {self.seq}'
        else:
            strform += f'\tseq: {self.seq[:_SEQEND_SIZE]}...{self.seq[-_SEQEND_SIZE:]}'
        return strform
    
    def reverse_complement(self, 
                           id=True, 
                           name=True, 
                           description=False, 
                           features=True, 
                           annotations=False, 
                           letter_annotations=True, 
                           dbxrefs=False):
        return BioBlock(super().reverse_complement(id, 
                                                    name, 
                                                    description, 
                                                    features, 
                                                    annotations, 
                                                    letter_annotations, 
                                                    dbxrefs))

'''
    A classe Primer é derivada de SeqRecord de Biopython. Essa classe tem
    como propósito representar uma sequência biológica que sirva como primer
    ou outro adaptador utilizado para modificar BioBlocks.
'''
class Primer(SeqRecord):

    def __init__(self, 
                seq, 
                id="<unknown id>", 
                name="<unknown name>", 
                description="<unknown description>", 
                dbxrefs=None, 
                features=None, 
                annotations=None, 
                letter_annotations=None,
                melting_temperature=None):
        if isinstance(seq, SeqRecord):
            self.__dict__.update(seq.__dict__)
            self.seq = Seq(self.seq) 
        elif isinstance(seq, (str, Seq)):
            super().__init__(Seq(seq), 
                            id, 
                            name, 
                            description, 
                            dbxrefs, 
                            features, 
                            annotations, 
                            letter_annotations)
        else:
            raise TypeError('constructor might receive a SeqRecord object or sequence string')
        self.melting_temperature = melting_temperature
    
    def __str__(self):
        slen = len(self.seq)
        strform = f'Primer {self.id}\n\tsize: {slen}\n'
        if self.melting_temperature:
            strform += f'\ttemp-melting (°C): {self.melting_temperature}\n'
        if slen <= _TERMINAL_SIZE:
            strform += f'\tseq: {self.seq}'
        else:
            strform += f'\tseq: {self.seq[:_SEQEND_SIZE]}...{self.seq[-_SEQEND_SIZE:]}'
        return strform

    def reverse_complement(self, 
                            id=True, 
                            name=True, 
                            description=True, 
                            features=True, 
                            annotations=True, 
                            letter_annotations=True, 
                            dbxrefs=True):
        return Primer(super().reverse_complement(id, 
                                                name, 
                                                description, 
                                                features, 
                                                annotations, 
                                                letter_annotations, 
                                                dbxrefs))

'''
    A classe DnaAssembly é um genérico para representar um conjunto de blocos
    biológicos a serem unidos por uma técnica qualquer de DNA Assembly. Os
    blocos são unidos sequencialmente.
'''
class DnaAssembly:

    def __init__(self, records, circular=False):
        self.circular = circular
        self.records = []
        if hasattr(records, "__iter__"):
            for item in records:
                if isinstance(item, BioBlock):
                    self.records.append(item)
                elif isinstance(item, (SeqRecord, str)):
                    self.records.append(BioBlock(item))
                else:
                    raise TypeError('incompatible sequence object in given iterable')
        else:
            raise TypeError('constructor might receive iterable of sequence objects')



if __name__ == "__main__":

    a = BioBlock('ATTACCGAAACCTTGCATGC', 'AXRT43.554')
    b = SeqRecord('ATCAGCCACTGATGCAATTTCTGTGTTTCACTGACTGCACAGGTTAAGGAAGGCCTCCAAATGTGCATCGGAAACTACTACTGGTCTTCCAGGGAAAAACTGACTAC', 'AXTR553.6')
    c = BioBlock(b)
    # BioBlock pode ser contruído a partir de sequência ou de SeqRecord
    print(a)
    print(c.reverse_complement())

    d = Primer('GCTTCGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT')
    e = d.reverse_complement()
    print(d)
    print(e)
    f = 'AAAATGTAGGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGC'[::-1]
    print(f)

