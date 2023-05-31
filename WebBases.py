from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import xmltodict
import requests


def temporary_email():
    try:
        email = requests.post(
                    'https://www.1secmail.com/api/v1/?action=genRandomMailbox&count=1'
                ).json()[0]
    except:
        email = "514adm2s0c@wwjmp.com"
    return email


def from_genbank(record, email=temporary_email()):
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=record, rettype='fasta', retmode='text')
    return SeqIO.read(handle, 'fasta')


def from_ensembl(record):
    try:
        data = requests.get(
                        f"https://rest.ensembl.org/sequence/id/{record}", 
                        headers= {"Content-Type": "application/json"}
                    ).json()
        return SeqRecord(Seq(data['seq']), 
                         id=data['id'], 
                         name=data['query'], 
                         description=data['desc'])
    except:
        return None


def from_igem(record):
    try:
        data = xmltodict.parse(
                    requests.get(f"https://parts.igem.org/xml/part.{record}").text
                )['rsbpml']['part_list']['part']
        return SeqRecord(Seq(data['sequences']['seq_data']),
                         id=data['part_name'],
                         name=data['part_name'],
                         description=data['part_short_desc'])
    except:
        return None


if __name__ == '__main__':
    # Exemplos de teste
    print(from_genbank('LR134385'), end='\n\n')
    print(from_ensembl('ENSG00000157764'), end='\n\n')
    print(from_ensembl('ENSDARG00000024771'), end='\n\n')
    print(from_igem('BBa_J23100'), end='\n\n')
    print(from_igem('BBa_R0010'), end='\n\n')
    print(from_igem('BBa_B0034'), end='\n\n')
    print(from_igem('BBa_E0040'), end='\n\n')
    print(from_igem('BBa_I0462'), end='\n\n')
    print(from_igem('BBa_K1170001'), end='\n\n')
    print(from_igem('BBa_F2620'), end='\n\n')
    print(from_igem('BBa_C0040'), end='\n\n')
    print(from_igem('BBa_J04450'), end='\n\n')
    print(from_igem('BBa_E0840'), end='\n\n')

