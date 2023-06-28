from typing import Union, List, Tuple


def suffix_array(string: Union[str, bytes, bytearray]) -> List[int]:
    '''
        Creates a suffix array of a string or similar sequence of characters.

        :param string: string or similar sequence of characters
        :return List[int]: suffix array
    '''
    return sorted(range(len(string)), key=lambda x: string[x:])


def lcp_array(text: Union[str, bytes, bytearray], suffix: List[int]) -> List[int]:
    '''
        Creates Longest Common Prefix (LCP)-array of a string or similar sequence
        of characters from its suffix array.

        :param text: string or similar sequence of characters
        :param suffix: suffix array of given text
        :return List[int]: LCP-array
    '''
    rank, height = [0] * len(suffix), [0] * len(suffix)
    for i in range(len(suffix)):
        rank[suffix[i]] = i
    h = 0
    for i in range(len(suffix)):
        if rank[i] > 0:
            k = suffix[rank[i] - 1]
            try:
                while text[i + h] == text[k + h]:
                    h += 1
            except:
                pass
            finally:
                height[rank[i]] = h
                if h > 0:
                    h -= 1
    return height


def longest_overlap(a: Union[str, bytes, bytearray], 
                    b: Union[str, bytes, bytearray]) -> Tuple[int,int,int]:
    '''
        Returns the absolute longest common sub-string between two strings a and b.
        The longest match is a tuple of 3 values: index that substring starts at a, 
        index that substring starts at b and length of substring.

        :param a: first string
        :param b: second string
        :return Tuple[int,int,int]: match
    '''
    concat = f'{a}#{b}'
    suffix_c = suffix_array(concat)
    lcp_c = lcp_array(concat, suffix_c)
    max_inda, max_indb, max_size = None, None, 0
    for x in range(1, len(lcp_c)):
        if lcp_c[x] > max_size:
            lowind = min((suffix_c[x - 1], suffix_c[x]))
            highind = suffix_c[x - 1] + suffix_c[x] - lowind
            if lowind < len(a) < highind:
                highind = highind - len(a) - 1
                max_inda, max_indb, max_size = lowind, highind, lcp_c[x]
    return max_inda, max_indb, max_size


def all_overlaps(a: Union[str, bytes, bytearray], 
                 b: Union[str, bytes, bytearray], 
                 minsize: int =10) -> List[Tuple[int,int,int]]:
    '''
        Return a list of all distinct non-overlapping common substrings between
        two strings or similar sequence of characters. Repeated subsequences
        inside the same string are not counted. Substrings can be equal and have
        same lenght, as long as they don't overlap.
        Matchs are tuples of 3 int values: index of substring at A, index of
        substring at B and length of substring.

        :param a: first string 
        :param b: second string 
        :minsize: minimal size of substring
        :return List[Tuple[int,int,int]]: list of matchs
    '''
    concat = f'{a}#{b}'
    suffix_c = suffix_array(concat)
    lcp_c = lcp_array(concat, suffix_c)
    matchs = set()
    for x in range(len(lcp_c)):
        xp = x + 1
        while xp < len(lcp_c) and lcp_c[xp] >= minsize:
            lowind = min((suffix_c[x], suffix_c[xp]))
            highind = suffix_c[x] + suffix_c[xp] - lowind
            if lowind < len(a) < highind:
                matchs.add((lowind, highind - len(a) - 1, lcp_c[xp]))
            xp += 1
    return [x for x in matchs if (x[0] - 1, x[1] - 1, x[2] + 1) not in matchs]


def end_overlaps(a: Union[str, bytes, bytearray], 
                 b: Union[str, bytes, bytearray], 
                 minsize: int = 10) -> List[Tuple[int,int,int]]:
    '''
        Returns a list of all possible terminal overlaps between two strings or
        similar character sequences, a and b.

        :param a: first string 
        :param b: second string 
        :minsize: minimal size of terminal overlap
        :return List[Tuple[int,int,int]]: list of overlaps
    '''
    a, b = f'{a}', f'{b}'
    matchs = [(0, len(b) - x, x) for x in range(minsize, len(a) + 1) if a[:x] == b[-x:]]
    matchs.extend((len(a) - x, 0, x) for x in range(minsize, len(b) + 1) if b[:x] == a[-x:])
    return matchs


def are_rotations(a: Union[str, bytes, bytearray], 
                  b: Union[str, bytes, bytearray]) -> Tuple:
    '''
        Check if two strings a and b are rotations of each other or not.
        If false, this function returns a tuple that can be evaluated to False.
        If true, returns a tuple with two int values: index 0 for the first 
        character in string a, and index for the same character in string b. 
        This integer tuple can be valued to true.

        :param a: first string 
        :param b: second string 
        :return Tuple: if strings are rotations of each other or not
    '''
    a, b = f'{a}', f'{b}'
    if len(a) != len(b):
        return ()
    else:
        for x in range(len(b)):
            if b[x] == a[0] and b[x:] == a[:len(b)-x] and b[:x] == a[len(b)-x:]: 
                return 0, x
    return ()


if __name__ == '__main__':

    from Bio.Seq import Seq
    from timeit import timeit
    a = 'AAATGTCGTTAACAACCTGTTAAACACACTCTCTGGGGAAGGGGGGGG'
    b = 'GGGGGTCCCTGCCTATTACCCTAGAACTGCAAGGCTGGC'
    c = 'CCTGTTAAACACACTCTCTGGGGAAGGGGGGGGAAATGTCGTTAACAA'
    
    print(are_rotations(a, b))
    print(are_rotations(a, c))
    print(a[0:])
    print(c[33:])

    x = 'CCCACCTCGAAACTC'
    y = 'AAAGTTGTGTGAACCCACCTCGAAACTCCTGGTTCCACACTGTCA'
