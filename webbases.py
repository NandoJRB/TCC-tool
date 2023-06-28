from biopart import Biopart 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import xmltodict
import requests


def temporary_email() -> str:
    '''
        This function accesses the 1secmail API and generates a temporary and 
        disposable email. It's necessary for request nucleotide sequences from 
        Genbank if you does not want use your personal email.

        :return str: temporary email
    '''
    try:
        email = requests.post(
                    'https://www.1secmail.com/api/v1/?action=genRandomMailbox&count=1'
                ).json()[0]
    except:
        email = "514adm2s0c@wwjmp.com"
    return email


def from_genbank(record: str, email : str = temporary_email()):
    '''
        Request a nucleotide record from the Genbank database using biopython's 
        native Entrez module. It returns a Biopart object in case of a successful 
        search, otherwise it returns None.

        :param record: nucleotide sequence ID for search
        :param email: email for Entrez
        :return Biopart: requested nucleotide sequence
    '''
    try:
        Entrez.email = email
        handle = Entrez.efetch(db='nucleotide', id=record, rettype='fasta', retmode='text')
        return Biopart.from_seqrecord(SeqIO.read(handle, 'fasta'), circular=False)
    except:
        return None


def from_ensembl(record):
    '''
        Request a nucleotide record from the Ensembl database. It requires string 
        ids in Ensembl format. It returns a Biopart object in case of a successful 
        search, otherwise it returns None.

        :param record: nucleotide sequence ID for search
        :return Biopart: requested nucleotide sequence
    '''
    try:
        data = requests.get(
                        f"https://rest.ensembl.org/sequence/id/{record}", 
                        headers= {"Content-Type": "application/json"}
                    ).json()
        return Biopart(forward=Seq(data['seq']), 
                       id=data['id'], 
                       name=data['query'], 
                       description=data['desc'])
    except:
        return None


def from_igem(record):
    '''
        Request a nucleotide record from the iGem Registry of Standart Biological Parts. 
        This function is useful for request BioBricks parts from iGem. It returns a Biopart
        object in case of a successful search, otherwise it returns None.

        :param record: biobrick id for search (in iGem format)
        :return Biopart: requested nucleotide sequence
    '''
    try:
        data = xmltodict.parse(
                    requests.get(f"https://parts.igem.org/xml/part.{record}").text
                )['rsbpml']['part_list']['part']
        return Biopart(forward=Seq(data['sequences']['seq_data']),
                       id=data['part_name'],
                       name=data['part_name'],
                       description=data['part_short_desc'])
    except:
        return None


if __name__ == '__main__':

    # Test examples
    print(from_genbank('LR134385'))
    print(from_ensembl('ENSG00000157764'))    
    print(from_igem('BBa_E0840'))
