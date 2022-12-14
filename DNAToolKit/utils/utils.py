"utils functions for DNAToolKit"

def reverse_complement(pattern) -> str:
    """Returns the reverse complement of a DNA string

    Parameters
    ----------
    pattern : str
        DNA string

    Returns
    -------
    str
        reverse complement of pattern
    """
    bases = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    
    complement = ""
    
    for i in pattern:
        complement += bases[i]
    
    reverse_complement_ = complement[::-1]
    return reverse_complement_

def number_to_symbol(index: int) -> str:
    """Converts a number into a DNA basa based on lexicographic order"""
    numb_to_symb = {
        0: "A",
        1: "C",
        2: "G",
        3: "T"
    }
    return numb_to_symb[index]

def symbol_to_number(symbol: str) -> int:
    """Converts a DNA base into a int based on lexicographic order"""
    values = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
    }
    return values[symbol.upper()]

def prefix(pattern: str) -> str:
    """Returns all but the last letter of a string"""
    return pattern[:-1]

def suffix(pattern: str) -> str:
    """Returns a string without the first letter"""
    return pattern[1:]

def first_symbol(pattern: str) -> str:
    """Returns the first letter of a string"""
    return pattern[0]

def last_symbol(pattern: str) -> str:
    """Returns the last letter of a string"""
    return pattern[-1]

def quotient(index: int, k: int) -> int:
    return index // k

def remainder(index: int, k: int) -> int:
    return index % k

def number_to_pattern(index: int, k: int) -> str:
    """Compute a DNA string (DNA pattern) given an integer that represents its 
    position when all patterns of length k are ordered lexicographically

    Parameters
    ----------
    index : int
        position that the pattern would occupy when ordered lexicographically
    k : int
        length of the pattern

    Returns
    -------
    str
        DNA pattern
    """
    if k == 1:
        return number_to_symbol(index)
    
    prefix_index = quotient(index, 4)
    remainder_ = remainder(index, 4)
    symbol = number_to_symbol(remainder_)
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    
    return prefix_pattern + symbol

def pattern_to_number(pattern: str) -> int:
    """Given a DNA string, convert it to a number corresponding to
    its position when ordered lexicographically

    Parameters
    ----------
    pattern : str
        DNA string to convert to a number

    Returns
    -------
    int
    
    """
    if not pattern:
        return 0
    
    prefix_ = prefix(pattern)
    symbol = last_symbol(pattern)
    return 4 * pattern_to_number(prefix_) + symbol_to_number(symbol)

def minimum_skew(genome: str) -> str:
    """Defines the skew of a DNA string, denoted skew(genome) as the difference
    between the total number of occurrences of "G" and "C" in Genome
    
    Let Prefix_i(genome) denote the prefix (i.e., initial subtring) of genome of 
    length i

    Parameters
    ----------
    genome : str
        A DNA string

    Returns
    -------
    str
        all integers i minimizing skew(prefix_i(genome)) over all values of i from
        0 to len(genome) in a string
    
    Example
    --------
    skew(prefix_i("CATGGGCATCGGCCATACGCC")) are:
    0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    """
    skew = [0] # first value is always 0, since its the skew for a subtring with
               # length 0
    c_count = 0
    g_count = 0
    for i in range(0, len(genome)):
        if genome[i] == "C":
            c_count += 1
        elif genome[i] == "G":
            g_count += 1
        skew.append(g_count-c_count)
        
    indices = [str(i) for i, x in enumerate(skew) if x == min(skew)]
    return " ".join(indices)