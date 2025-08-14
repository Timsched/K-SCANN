import numpy as np
from Bio.Align import substitution_matrices


def OneHotEncoding(sequence:str):
    """
    One hot encoding of a sequence
    :param sequence: Amino acid sequence as a string
    :return: One hot encoding of the sequence as a np.array with dimensions (sequence length*21)
    """
    AminoAcids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','_'] #_ is for stop
    AminoAcidsDict = {AminoAcids[i]:i for i in range(len(AminoAcids))}
    
    sequencelength=len(sequence)
    encoding=np.zeros((len(AminoAcids),sequencelength),np.int8)

    for i in range(sequencelength):
        encoding[AminoAcidsDict[sequence[i]],i]=1


    return encoding.flatten("F")

def BlossumEncoding(sequence:str,similarity:int=62):
    """
    BLOSUM encoding of a sequence
    :param sequence: Amino acid sequence as a string
    :param sim: Similarity score to use for encoding (e.g. 50, 62, 80)
    :return: BLOSUM encoding of the sequence as a np.array with dimensions (sequence length*21)
    """
    AminoAcids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','_'] #X is uknown aa, _ is stop codon, Z is Glutamic acid or Glutamine, B is Aspartic acid or Asparagine
    
    blosum_sub = substitution_matrices.load("BLOSUM"+str(similarity))

    #drop unwanted amino acids (B,Z,X)
    blosum_sub_thin=np.copy(blosum_sub)
    blosum_sub_thin=np.delete(blosum_sub_thin,[AminoAcids.index('B'),AminoAcids.index('Z'),AminoAcids.index('X')],0)
    blosum_sub_thin=np.delete(blosum_sub_thin,[AminoAcids.index('B'),AminoAcids.index('Z'),AminoAcids.index('X')],1)

    trimmedAA=AminoAcids[0:20]+[AminoAcids[23]]


    encoding=np.zeros((len(trimmedAA),len(sequence)),np.int16)

    for i,letter in enumerate(sequence):
        
        encoding[:,i]=blosum_sub_thin[trimmedAA.index(letter),:]
 
    return encoding.flatten("F")#flatten column wise

def PairwiseEncoding(sequence:str):
    """
    Encoding the Aminoacids as pairs rather than single amino acids
    :param sequence: Amino acid sequence as a string
    :return: Pairwise encoding of the sequence as a np.array with dimensions (sequence length-1*21*21)
    """

    AminoAcids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','_']

    pairwise=np.zeros((len(AminoAcids),len(AminoAcids),len(sequence)-1),np.int8)

    for i,letter in enumerate(sequence[:-1]):
        
        pairwise[AminoAcids.index(letter),AminoAcids.index(sequence[i+1]),i]=1

    return pairwise.flatten("F")

def FuzzyBlossumEncoding(sequence:str,similarity:int=62,mean=0,std=1):
    """
    Encode the amino acids using BLOSSUM encoding with noise.
    :param sequence: Amino acid sequence as a string
    :return: Fuzzy encoding of the sequence as a np.array with dimensions (sequence length*21)
    """
    AminoAcids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','_'] #X is uknown aa, _ is stop codon, Z is Glutamic acid or Glutamine, B is Aspartic acid or Asparagine
    
    blosum_sub = substitution_matrices.load("BLOSUM"+str(similarity))

    #drop unwanted amino acids (B,Z,X)
    blosum_sub_thin=np.copy(blosum_sub)
    blosum_sub_thin=np.delete(blosum_sub_thin,[AminoAcids.index('B'),AminoAcids.index('Z'),AminoAcids.index('X')],0)
    blosum_sub_thin=np.delete(blosum_sub_thin,[AminoAcids.index('B'),AminoAcids.index('Z'),AminoAcids.index('X')],1)

    trimmedAA=AminoAcids[0:20]+[AminoAcids[23]]


    encoding=np.zeros((len(trimmedAA),len(sequence)),np.float16)

    for i,letter in enumerate(sequence):
        
        encoding[:,i]=blosum_sub_thin[trimmedAA.index(letter),:]+np.random.normal(mean,std,21)

    return encoding.flatten("F")#flatten column wise

def RandomPairwiseEncoding(sequence):
    """
    Random encoding of a sequence
    :return: Random encoding of the sequence as a np.array with dimensions (sequence length-1*21*21)
    """
    
    return np.random.rand(len(sequence)-1,21,21).flatten("F")