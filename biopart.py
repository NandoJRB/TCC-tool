'''
    Nesse arquivo, vou implementar os objetos padrão para manipulação de sequências.
    A princípio, serão apenas os 3 objetos: BioBlocks, Primer e Assembly. Estão
    sujeitos a mudança e eu posso colocá-los em arquivos distintos adiante.
    
'''


from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio.Restriction import RestrictionBatch, CommOnly, BamHI, EcoRI, SpeI, SfaNI, EcoRV
from compstrings import end_overlaps, are_rotations
from typing import Union, Optional, Literal



'''
    A classe Biopart é derivada de SeqRecord de Biopython. Essa classe tem
    como propósito representar uma sequência biologica qualquer que pode ser
    unida com outros BioBlocks em um Assembly

    O que essa classe precisa fazer:
        * Ela precisa poder ser criada a partir de um objeto SeqRecord de
          Biopython. Ela precisa poder ser criada a partir de uma string

        * Ela deve conter tanto a sequência forward, quanto a sequência
          reverse simultaneamente
        
        * Ela deve ter métodos simples, tais como: __str__, __hash__, len,
          reverse_complement, left_border, right_border, find, etc. Todos
          os métodos relevantes para ela como em pydna, dnacauldron e Biopython

        * Ela precisa ter propertys. 

'''


class Biopart:

    empty_seqrecord = SeqRecord('')
    penzymes = RestrictionBatch([x for x in CommOnly if x.site == Seq(x.site).reverse_complement()])

    def __init__(self, 
                 forward : Union[str, bytes, bytearray, Seq, MutableSeq] = None,
                 reverse : Union[str, bytes, bytearray, Seq, MutableSeq] = None,
                 id : Optional[str] = None, 
                 name : Optional[str] = None, 
                 description : Optional[str] = None, 
                 circular : bool = False,
                 overhang : int = None):
        # 1) Checks for forward and reverse parameters and converts it to Bio.Seq objects
        if not forward and not reverse:
            raise TypeError('forward and/or reverse parameters are required')
        if forward:
            forward = Seq(f'{forward}'.upper())
        else:
            if overhang and not circular:
                raise TypeError('given overhang parameter without a forward parameter')
            forward = Seq(f'{reverse}'.upper()).complement()
        if reverse:    
            reverse = Seq(f'{reverse}'.upper())
        else:
            if overhang and not circular:
                raise TypeError('given overhang parameter without a reverse parameter')
            reverse = forward.complement()
        # 2) Check if forward, reverse, overhang and circular parameters are compatible
        self._circular = circular
        if self._circular:
            temp = are_rotations(forward, reverse.complement())
            if temp:
                self._overhang, self._forward = 0, forward
                self._reverse = Seq(f'{reverse}'[temp[1]:] + f'{reverse}'[:temp[1]])
            else:
                temp = are_rotations(forward, reverse.reverse_complement())
                if temp:
                    self._overhang, self._forward, reverse = 0, forward, f'{reverse}'[::-1]
                    self._reverse = Seq(f'{reverse}'[temp[1]:] + f'{reverse}'[:temp[1]])
                else:
                    raise ValueError('incompatible strands for circular topology')
        else:
            if isinstance(overhang, int):
                temp = end_overlaps(forward,reverse.complement(), 1)
                for x in temp:
                    if overhang == (x[1] - x[0]):
                        self._overhang, self._forward, self._reverse = overhang, forward, reverse
                        break
                else:
                    temp = end_overlaps(forward, reverse.reverse_complement(), 1)
                    for x in temp:
                        if overhang == (x[1] - x[0]):
                            self._overhang, self._forward = overhang, forward
                            self._reverse = Seq(f'{reverse}'[::-1])
                            break
                    else:
                        raise ValueError('incompatible overhang value for given strands')
            else:
                if forward == reverse.complement():
                    self._overhang, self._forward, self._reverse = 0, forward, reverse
                else:
                    # Tenho que corrigir isso para aceitar fragmentos não sticky-ends e
                    # usar longest common substrings

                    temp = end_overlaps(forward, reverse.complement(), 1)
                    if temp:
                        temp = max(temp, key=lambda x: x[2])
                        self._forward, self._reverse = forward, reverse
                        self._overhang = temp[1] - temp[0]
                    else:
                        temp = end_overlaps(forward, reverse.reverse_complement(), 1)
                        if temp:
                            temp = max(temp, key=lambda x: x[2])
                            self._forward, self._reverse = forward, Seq(f'{reverse}'[::-1])
                            self._overhang = temp[1] - temp[0]
                        else:
                            raise ValueError('forward and reverse strands are incompatible')
        if id is None:
            id = Biopart.empty_seqrecord.id
        if name is None:
            name = Biopart.empty_seqrecord.name
        if description is None:
            description = Biopart.empty_seqrecord.description
        self._id = id
        self._name = name
        self._description = description

    def __str__(self) -> str:
        return f'<Biopart: {self._id}; len: {len(self)}; topology:' + \
            f' {"circular" if self._circular else "linear"}>'

    def __len__(self) -> int:
        return len(self.global_forward())

    def __eq__(self, __o) -> bool:
        try:
            return self._forward == __o.forward and \
                   self._reverse == __o.reverse and \
                   self._overhang == __o.overhang
        except:
            return False

    @classmethod
    def from_seqrecord(cls, 
                       seqrecord : SeqRecord, 
                       circular : bool = False):
        return cls(forward=seqrecord.seq, 
                   id=seqrecord.id, 
                   name=seqrecord.name, 
                   description=seqrecord.description, 
                   circular=circular)

    @property
    def forward(self):
        return self._forward

    @property
    def reverse(self):
        return self._reverse      

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def description(self):
        return self._description

    @property
    def circular(self):
        return self._circular
    
    @property
    def overhang(self):
        return self._overhang

    def strands_alignment(self):
        if self._overhang < 0:
            return f'{self._forward}\n{" " * abs(self._overhang)}{self._reverse}'
        else:
            return f'{" " * self._overhang}{self._forward}\n{self._reverse}'
    
    def copy(self):
        return Biopart(forward=self._forward,
                       reverse=self._reverse,
                       id=self.id,
                       name=self.name,
                       description=self.description,
                       circular=self.circular
                    )

    def invert_strands(self):
        return Biopart(forward=self._reverse,
                       reverse=self._forward,
                       id=self.id,
                       name=self.name,
                       description=self.description,
                       circular=self.circular
                    )

    def global_forward(self):
        if self._overhang < 0:
            return self._forward[:abs(self._overhang)] + self._reverse.complement()
        else:
            return self._reverse.complement()[:self._overhang] + self._forward

    def global_reverse(self):
        return self.global_forward().complement()

    def right_end(self):
        # Não converte a fita reversa para o seus complemento
        if self._overhang < 0:
            return self._reverse[len(self._forward[abs(self._overhang):]):]
        else:
            return self._forward[len(self._reverse[self._overhang:]):]

    def left_end(self):
        # Não converte a fita reversa para o seus complemento
        if self._overhang < 0:
            return self._forward[:abs(self._overhang)]
        else:
            return self._reverse[:self._overhang]

    def double_strand(self):
        # Retorna apenas a porção da sequencia que possui fita dupla
        if self._overhang < 0:
            return self._forward[abs(self._overhang):]
        else:
            return self._reverse[self._overhang:].complement()
    
    def add_tag_forward(self, 
                        tag: Union[str, bytes, bytearray, Seq, MutableSeq], 
                        direction : Literal['right','left'] = 'right'):
        if direction == 'right':
            if self._overhang < 0:
                pass
            else:
                pass
        elif direction == 'left':
            if self._overhang < 0:
                pass
            else:
                pass
        else:
            raise ValueError('direction parameter value is not compatible')

    def add_tag_reverse(self, 
                        tag: Union[str, bytes, bytearray, Seq, MutableSeq], 
                        direction : Literal['right','left'] = 'right'):
        if direction == 'right':
            if self._overhang < 0:
                pass
            else:
                pass
        elif direction == 'left':
            if self._overhang < 0:
                pass
            else:
                pass
        else:
            raise ValueError('direction parameter value is not compatible')
    
    def has_blunt_ends(self) -> bool:
        return self.left_end() == Seq('') == self.right_end() and not self.circular

    def compatible_ends(self) -> bool:
        # Não checa se tem extremidades cegas
        if self._overhang < 0:
            return self.left_end() == self.right_end().complement()
        else:
            return self.left_end().complement() == self.right_end()

    def is_homologous(self, __o, fromleft=False) -> bool:
        # Isso não checa se os fragmentos tem blunt ends, retornará True nesse caso
        if fromleft:
            return __o.right_end() == self.left_end().complement()
        else:
            return self.right_end() == __o.left_end().complement()

    def join_ends(self):
        # Não tem efeito quando já é circular
        if not self._circular:
            if not self.compatible_ends():
                raise AssertionError('biopart ends are not compatible')
            self._circular, self._overhang = True, 0
            self._reverse = self._forward.complement()
    
    def assemble(self, __o):
        pass

    def search_for_all(self, minsites: Optional[int] = None, maxsites: Optional[int] = None):
        # considera apenas a porção da sequencia que possui dupla fita
        if not minsites:
            minsites = 0
        if not maxsites:
            maxsites = float('inf')
        enz_dict, gseq = dict(), self.double_strand()
        for enz in Biopart.penzymes:
            cutf = enz.search(gseq)
            if minsites <= len(cutf) <= maxsites:
                enz_dict[str(enz)] = cutf
        return enz_dict 
            
    def find_sites_enzyme(self, enzyme):
        # considera apenas a porção da sequencia que possui dupla fita
        cutenz = Biopart.penzymes.get(enzyme)
        if cutenz is not None:
            return cutenz.search(self.double_strand())
        else:
            raise TypeError('invalid enzyme name for search')

    def cut_with_enzyme(self, enzyme):
        # Eu deveria simplesmente cortar, mesmo que não fosse circular
        cutenz = Biopart.penzymes.get(enzyme)
        if cutenz is not None:
            seqdf = self.double_strand()
            catf = cutenz.catalyse(seqdf)
            if len(catf) == 1:
                return (self.copy(),)
            else:
                catr = cutenz.catalyse(seqdf.reverse_complement())
                retl = [Biopart(x, y[::-1]) for x, y in zip(catf, reversed(catr))]
                retl[0].add_tag_left(self.left_end())
                retl[-1].add_tag_right(self.right_end())
                return tuple(retl)
                

if __name__ == "__main__":

    fb = 'ATCAGCCACTGATGCAATTTCTGTGTTTCACTGACTGCACAGGTT'
    rb = 'TAGTCGGTGACTACGTTAAAGACACAAAGTGACTGACGTGTCCAA'
    # Adicionei um segundo sítio de BamHI e deu pau, implementar métodos add_tail
    bv = 'CAGCCACTGATGATATCGCAAGGATCCTTCTGTCCCGCAGGTTAAAGACGGATCCACAAACTACCGGGTTGTCCACGGG'
    bg = 'AATCCTGTAAAAGTCGGTGACTACTATAGCGTTCCTAGGAAGACAGGGCGTCCAATTTCTGCCTAGGTGTTTGATGGC'
    b1 = 'CCCCTTT'
    b2 = 'CCCCAAA'
    b3 = 'ATTGTTCCAACTCGAT'
    b4 = 'GGTTGAGCTATAACAA'
    b5 = b4[::-1]
    b6 = 'AATGGTGCCTACCGCCTTCCCCTGAACTACGGCA'
    b7 = 'CTACGGCAGGTCCAACTGTTCCACCTGGTAACTCCGGGACT'


    z = Biopart(bv, bg)
    print(z.strands_alignment())
    #print(z.overhang)

    for x in z.cut_with_enzyme(CommOnly.get('BamHI')):
        print()
        print(x.strands_alignment())
        print(x.overhang)



