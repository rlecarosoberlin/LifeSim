"""
Baer Hodge
Jan 7, 2025
"""


import argparse
import copy
import os
import random
from argparse import Namespace

import matplotlib.pyplot as plt
import numpy as np

from food import Food
from organism import Organism
from collections import Counter

import time


def unit_interval(x) -> float:
    """
    A float between 0 and 1.

    Parameters
    ----------
    x
        The value whose datatype and range are being checked.

    Returns
    -------
    x : float
        The input converted into a float.

    Raises
    ------
    ValueError
        If `x` cannot be converted into float.
    argparse.ArgumentTypeError
        If `x` is not between 0 and 1.

    """
    # if `x` is not float, raise error
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{x} is not a valid variable type for this argument. Must be a "
            f"float between 0 and 1")

    # if `x` is not between 0 and 1, raise error
    if not 0 <= x <= 1:
        raise argparse.ArgumentTypeError(
            f"{x} is out of range for this argument. Must be a float between "
            f"0 and 1")

    return x


def pos_int(x) -> int:
    """
    An integer not less than 0.

    Parameters
    ----------
    x
        The value whose datatype and range are being checked.

    Returns
    -------
    x : int
        The input converted into an integer.

    Raises
    ------
    ValueError
        If `x` cannot be converted into integer.
    argparse.ArgumentTypeError
        If `x` is less than zero.
    """
    # if `x` is not integer, raise error
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{x} is not a valid variable type for this argument. Must be 0 "
            f"or a positive integer")

    # if `x` is less than 0, raise error
    if x < 0:
        raise argparse.ArgumentTypeError(
            f"{x} is out of range for this argument. Must be 0 or a positive "
            f"integer")

    return x


def two_comma_separated_ints(x) -> tuple[int, int]:
    """
    Convert a string of two integers separated by a comma into a tuple of
    those two integers.

    Parameters
    ----------
    x : str
        Value of which the format is being checked.

    Returns
    -------
    tuple[int, int]
        A tuple of the two integers from the inputted string.

    Raises
    ------
    argparse.ArgumentTypeError
        If `x` is not a string of two integers separated by a comma

    Examples
    --------
    # >>> two_comma_separated_ints("3,4")
    (3, 4)

    use three numbers
    # >>> two_comma_separated_ints("3,4,5")
    argparse.ArgumentTypeError

    use a float or a letter
    # >>> two_comma_separated_ints("3.1,Y")
    argparse.ArgumentTypeError

    """
    # Split the input string by the comma, stored as a list
    nums: list[str] = x.split(',')

    # if `nums` has two elements and both are convertible to integers
    if len(nums) == 2:
        if nums[0].isdigit() and nums[1].isdigit():
            return int(nums[0]), int(nums[1])
    else:
        raise argparse.ArgumentTypeError(
            f"{x} is not a valid input for this argument. Must be a "
            f"string of two integers separated by a comma with no "
            f"whitespace.")


def parse_through_args() -> Namespace:
    """
    Assign arguments from commandline to attributes of an object.

    Returns
    -------
    args : Namespace
        Object with all arguments from commandline stored in its attributes.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""\
    This program simulates a population of Organism objects living in 
    1-dimensional space. Organisms have genome with one or more genes that can 
    mutate. Organisms can replicate, die, eat, move, and undergo horizontal 
    gene transfer (HGT) via both conjugation and transformation. Please input 
    the proper values into the commandline as arguments for the parameters 
    specified below.\
        """)

    parser.add_argument('--outputDirectory', type=str, required=True,
                        help="""\
    Enter STRING in following format: /path/to/my/directory/ 
                            
    The directory path to which all output files are written. There is no need 
    to change this value across the multiple times you run the simulation with 
    different arguments. 
                        """)

    parser.add_argument('--runID', type=str, required=True,
                        help="""\
    Enter STRING.
                        
    The unique name to give to this specific run of the simulation. All output 
    files from this run will have this value in their filenames, and will all 
    be stored in a folder named with this value. 
                        """)

    parser.add_argument('--pGenSize', type=pos_int, required=True,
                        help="""\
    Enter POSITIVE INTEGER.

    The number of Organisms in the parent generation. These Organisms are 
    instantiated before the first timestep.
                        """)

    parser.add_argument('--censusTimes', type=pos_int, required=True,
                        nargs='+',
                        help="""\
    Enter one or many POSITIVE INTEGERS.

    The timesteps when the living population of Organisms is surveyed and the 
    data are analyzed. 
                        """)

    parser.add_argument('--mutationRate', type=unit_interval, required=True,
                        help="""\
    Enter FLOAT BETWEEN 0 AND 1.

    For all Organisms in parent generation, the probability of various 
    attributes of a newly instantiated child Organism to be different from 
    those of its parent. Also, the probability of each individual nucleotide 
    in an Organism's genome to get copied incorrectly when replicating or 
    spontaneously mutating. 
                        """)

    parser.add_argument('--hgtConjRate', type=unit_interval, required=True,
                        help="""\
    Enter FLOAT BETWEEN 0 AND 1.

    For all Organisms in parent generation, the probability of an Organism 
    undergoing horizontal gene transfer (HGT) via conjugation in a given 
    timestep, if criteria are met. This involves the Organism being the DONOR 
    of genetic material to another nearby Organism, who acts as the recipient 
    of genetic material.        
                        """)

    parser.add_argument('--hgtTranRate', type=unit_interval, required=True,
                        help="""\
    Enter FLOAT BETWEEN 0 AND 1.
                        
    For all Organisms in parent generation, the probability of an Organism 
    undergoing horizontal gene transfer (HGT) via transformation in a given 
    timestep, if criteria are met. This involves the Organism being the 
    RECIPIENT of genetic material from a dead Organism who has left its genome 
    loose in the environment.     
                        """)

    parser.add_argument('--motility', type=pos_int, default=2,
                        help="""\
    Enter POSITIVE INTEGER.
    (default = 2)

    For all Organisms in parent generation, the maximum by which an Organism's 
    xCoord can change in a given timestep. 
                        """)

    parser.add_argument('--hgtMaxDistance', type=pos_int, default=1,
                        help="""\
    Enter a POSITIVE INTEGER. 
    (default = 1)

    The maximum distance that can separate two Organisms while still 
    allowing for HGT conjugation. 
                             """)

    parser.add_argument('--replDeathRate', type=unit_interval, default=0.3,
                        help="""\
    Enter FLOAT BETWEEN 0 AND 1. 
    (default = 0.3)
    
    For all Organisms in parent generation, the probability of an Organism 
    giving rise to a new Organism in a given timestep, all of whose attributes 
    are derived from the parent. 
    
    Also the probability of an Organism dying in a 
    given timestep, and leaving its genome loose in the environment at the 
    location of death.
                        """)

    parser.add_argument('--spaceSize', type=pos_int, default=50,
                        help="""\
    Enter POSITIVE INTEGER. 
    (default = 50)

    The size of the 1-dimensional space of the simulation.
                        """)

    parser.add_argument('--foodSources', type=two_comma_separated_ints,
                        default='10,500 40,500', nargs='+',
                        help="""\
    Enter one or many COMMA-SEPARATED PAIRS OF INTEGERS. 
    (default = 10,500 40,500)
    
    The food sources in the simulation. The first number in a comma-separated 
    pair of integers is the location of the food source in space, and the 
    second number is how much food is at the food source.
                        """)

    parser.add_argument('--foodShelfLife', type=two_comma_separated_ints,
                        default='20,40',
                        help="""\
    Enter a COMMA-SEPARATED PAIR OF INTEGERS.
    (default = 20,40)
    
    The minimum and maximum age a Food can be before respawning and resetting
    age to 0.
                        """)

    parser.add_argument('--hungerTolerance', type=pos_int, default=1,
                        help="""\
    Enter a POSITIVE INTEGER.
    (default = 1)
    
    The number of timesteps an Organism can go without eating Food before it
    dies.
                        """)

    parser.add_argument('--fedReplBoost', type=unit_interval, required=True,
                        help="""\
    Enter FLOAT BETWEEN 0 and 1.
    
    The temporary bonus added to an Organism's replication rate when its 
    hunger is 0. 
                        """)

    parser.add_argument('--plotYLim', type=pos_int, default=None,
                        help="""\
    Enter a POSITIVE INTEGER.
    (default = None)
    
    If specified, the maximum value on the Y-axis of the Organisms and Food 
    plot. If not specified, there is no limit.
                        """)

    parser.add_argument('--geneLength', type=pos_int, default=1500,
                        help="""\
    Enter a POSITIVE INTEGER.
    (default = 1500)
    
    The number of nucleotides in each gene of the genome of each Organism. 
    Only called when not using a predetermined genome.
                        """)

    parser.add_argument('--geneNumber', type=pos_int, default=1,
                        help="""\
    Enter a POSITIVE INTEGER.
    (default = 1)

    The number of genes in the genome of each Organism. Only called when not 
    using a predetermined genome.
                        """)

    parser.add_argument('--replRateChangeMax', type=unit_interval, default=0,
                        help="""\
    Enter a FLOAT BETWEEN 0 AND 1.
    (default = 0)
    
    The maximum by which an Organism's replication rate can differ from that 
    of its parent's.
                        """)

    parser.add_argument('--deathRateChangeMax', type=unit_interval, default=0,
                        help="""\
    Enter a FLOAT BETWEEN 0 AND 1.
    (default = 0)
    
    The maximum by which an Organism's death rate can differ from that 
    of its parent's.
                        """)

    parser.add_argument('--motilityChangeMax', type=pos_int, default=0,
                        help="""\
    Enter a POSITIVE INTEGER.
    (default = 0)
    
    The maximum by which an Organism's motility can differ from that 
    of its parent's.
                        """)

    parser.add_argument('--mutationRateChangeMax', type=unit_interval,
                        default=0,
                        help="""\
    Enter a FLOAT BETWEEN 0 AND 1.
    (default = 0)
    
    The maximum by which an Organism's mutation rate can differ from that 
    of its parent's.
                        """)

    parser.add_argument('--hgtTranRateChangeMax', type=unit_interval,
                        default=0,
                        help="""\
                        
    Enter a FLOAT BETWEEN 0 AND 1.
    (default = 0)
    
    The maximum by which an Organism's HGT transformation rate can differ from 
    that of its parent's.
                        """)

    parser.add_argument('--hgtConjRateChangeMax', type=unit_interval,
                        default=0,
                        help="""\
    Enter a FLOAT BETWEEN 0 AND 1.
    (default = 0)
    
    The maximum by which an Organism's HGT conjugation rate can differ from 
    that of its parent's.
                        """)

    parser.add_argument('--foodGaussMean', type=float, default=0,
                        help="""\
    Enter a FLOAT.
    (default = 0)
    
    The mean of a Gaussian distribution, indicating the average distance Food 
    will random-walk in each timestep.
                        """)

    parser.add_argument('--foodGaussSigma', type=float, default=1,
                        help="""\
    Enter a FLOAT.
    (default = 0)

    The standard deviation of a Gaussian distribution, indicating the 
    variability in the distance Food will random-walk in each timestep.
                        """)

    parser.add_argument('--genomeType', type=str, required=True,
                        choices=['eColi16s', 'random', 'relatedRandom'],
                        help="""\
    Enter STRING equal to one of the CHOICES below.
    
    The manner in which the genomes of the Organisms in the parent generation
    are generated.
    
     CHOICES
     --------------------------------------------------------------------------
     eColi16s - All Organisms in parent generation have a singular gene, 
                which is E. coli 16s rRNA sequence from the literature.
                   
       random - All Organisms in parent generation have randomly generated 
                genomes, with gene length and number of genes specified from 
                commandline.
                
relatedRandom - All Organisms in parent generation are mutated versions
                of a randomly generated genome, resulting in genetic 
                diversity of the parent generation being proportional 
                to the mutation rate of the Organisms in the parent 
                generation.
                        """)

    # assign arguments from commandline to attributes of `args`
    args: Namespace = parser.parse_args()

    # some attributes of `args` must be changed in case user entered them
    # incorrectly.

    # Convert to set to remove duplicates. Convert to sorted list (the largest
    # element must be last to know when simulation should stop). Then convert
    # to tuple to ensure immutability.
    args.censusTimes = tuple(sorted(set(args.censusTimes)))

    # Convert to sorted list because this later acts as an interval.
    # Then convert back into tuple to ensure immutability.
    args.foodShelfLife = tuple(sorted(args.foodShelfLife))

    print(f'\n{args=}')

    return args


args = parse_through_args()


def initialize_food(foodSources) -> dict[int, Food]:
    """
    Make a dictionary that stores Food objects.

    Parameters
    ----------
    foodSources : list[tuple[int, int]]
        Encodes how much Food has what value as spawnPoint attribute.
        First number in each tuple is spawnPoint and second number is how much
        Food to instantiate.

    Returns
    -------
    foods : dict[int, Food]
        Dictionary in which Foods are stored as values.

    Examples
    --------
    # >>> initialize_food_dict([(10, 2), (40,3)])
    {1:Food(10), 2:Food(10), 3:Food(40), 4:Food(40), 5:Food(40)}

    """
    foodsKeySetter = 0
    foods = {}

    # all elements in `foodSources` are tuples of two integers.
    # tuple[0] sets spawnPoint attribute, while tuple[1] is how much Food to
    # instantiate with that spawnPoint attribute
    for spawnPoint, amount in foodSources:
        for _ in range(amount):
            foods[foodsKeySetter] = Food(spawnPoint)

            foodsKeySetter += 1
    return foods


def initialize_p_gen(genomeType=args.genomeType,
                     pGenSize=args.pGenSize,
                     mutationRate=args.mutationRate,
                     motility=args.motility,
                     replRate=args.replDeathRate,
                     deathRate=args.replDeathRate,
                     hgtConjRate=args.hgtConjRate,
                     hgtTranRate=args.hgtTranRate,
                     spaceSize=args.spaceSize,
                     geneNumber=args.geneNumber,
                     geneLength=args.geneLength)\
                     -> dict[int, Organism]:
    """
    A dictionary that stores Organism objects.

    Parameters
    ----------
    genomeType : str, optional
        The manner in which the genome attribute of Organisms are generated.
    pGenSize : int, optional
        The number of Organisms to instantiate.
    mutationRate : float, optional
        The value to assign to Organisms' mutationRate attribute. Affects
        genetic diversity.
    motility : int, optional
        The value to assign to Organisms' motility attribute. Affects movement
        of Organisms
    replRate : float, optional
        The value to assign to Organisms' replRate attribute. Affects frequency
        of Organisms giving rise to new Organisms.
    deathRate : float, optional
        The value to assign to Organisms' deathRate attribute. Affects
        frequency of Organism death.
    hgtConjRate : float, optional
        The value to assign to Organisms' hgtConjRate attribute. Affects
        frequency of donating genes via HGT conjugation.
    hgtTranRate : float, optional
        The value to assign to Organisms' hgtTranRate attribute. Affects
        frequency of receiving genes via HGT transformation.
    spaceSize : int, optional
        The maximum an Organism's xCoord attribute can be set to.
        Globally, the size of the 1-dimensional space of the simulation.
    geneNumber: int, optional
        The number of genes in the genome attribute of Organisms.
    geneLength: int, optional
        The number of nucleotides in each gene of the genome attribute of
        Organisms.

    Returns
    -------
    pGenOrgs : dict[int, Organism]
        Dictionary in which Organisms are stored as values.

    Notes
    -----
    `geneNumber` and `geneLength` do not affect genome attribute of Organisms
    if `genomeType` equals "eColi16s". In this case, each Organism has one only
    one gene in genome attribute with the sequence being that of the
    E. coli 16s rRNA.

    """
    pGenOrgs = {}
    # these will be the keys of pGenOrgs as well as the idNumbs of the orgs
    idNumbSetter = 0
    # This is the sequence for E. coli 16s rRNA from the literature.
    rRNA_eColi16sSeq = [
        list("AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAA"
             "GTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAA"
             "TGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCAT"
             "AACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATG"
             "GGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGA"
             "GGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGG"
             "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCT"
             "TCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATT"
             "GACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAG"
             "GGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA"
             "GATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTC"
             "GTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACC"
             "GGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCA"
             "AACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCC"
             "CTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCA"
             "AGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAAT"
             "TCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAG"
             "AATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGA"
             "AATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGC"
             "CGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTC"
             "ATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCG"
             "ACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAAC"
             "TCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGT"
             "TCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGT"
             "AGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAA"
             "CAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA")]

    # if genomeType=='relatedRandom', then all parent generation Organisms have
    # genomes that are replicated (and mutated) from this `grandparent`
    # Organism. Otherwise, it is not used.
    grandparent = Organism(rand_genome(geneNumber, geneLength),
                           mutationRate,
                           0,  # genNumb
                           random.randint(1, spaceSize),  # xCoord
                           motility,
                           0,  # birthday
                           0,  # idNumb
                           replRate,
                           deathRate,
                           hgtConjRate,
                           hgtTranRate)

    # make parent generation
    for _ in range(pGenSize):
        # generate genome for new org depending on genomeType
        if genomeType == 'eColi16s':
            genome = rRNA_eColi16sSeq
        elif genomeType == 'relatedRandom':
            genome = grandparent.replicate_genome()
        elif genomeType == 'random':
            genome = rand_genome(geneNumber, geneLength)
        else:
            raise ValueError("parameter `genomeType` must be in either "
                             "'eColi16s', 'relatedRandom', or 'random'")

        # instantiate new org
        pGenOrgs[idNumbSetter] = \
            Organism(genome,
                     mutationRate,
                     0,  # genNumb
                     random.randint(1, spaceSize),  # xCoord
                     motility,
                     0,  # birthday
                     idNumbSetter,
                     replRate,
                     deathRate,
                     hgtConjRate,
                     hgtTranRate)
        idNumbSetter += 1
    del grandparent
    return pGenOrgs


def food_diffusion(gaussMean=args.foodGaussMean,
                   gaussSigma=args.foodGaussSigma,
                   shelfLife=args.foodShelfLife,
                   spaceSize=args.spaceSize):
    """
    Simulate diffusion of Food objects with a Gaussian random-walk, by changing
    each Food's xCoord attribute. If Food gets old or goes out of bounds, reset
    Food.

    Parameters
    ----------
    spaceSize: int, optional
        The maximum a Food's xCoord attribute can be.
        Globally, the size of the 1-dimensional space of the simulation.
        The size of the 1-dimensional space of the simulation.
    shelfLife : tuple[int,int], optional
        The minimum and maximum a Food's age attribute can be before the Food
        object's age and xCoord attributes reset.
    gaussMean : float, optional
        The mean of a gaussian distribution. Affects average distance Food
        objects travel each timestep.
    gaussSigma : float, optional
        The standard deviation of a gaussian distribution. Affects variation of
        distance Food objects travel in each timestep.

    Notes
    -----
    `shelfLife` is two numbers acting as an interval instead of one constant
    value to avoid all Food resetting at once, which would be very deadly to
    Organisms not near a food source.

    """
    for food in foods.values():

        # if Food expires due to age attribute, reset age and xCoord attributes
        if food.age >= random.randint(shelfLife[0], shelfLife[1]):
            food.age = 0
            food.xCoord = food.spawnPoint

        # move the Food by changing xCoord attribute, and increase age
        step = round(np.random.normal(gaussMean, gaussSigma))
        food.xCoord += step
        food.age += 1

        # if Food is out of bounds, reset age and xCoord attribute
        if not (1 <= food.xCoord <= spaceSize):
            food.age = 0
            food.xCoord = food.spawnPoint


def mutate_stats(org: Organism,
                 mutationRateChangeMax=args.mutationRateChangeMax,
                 replRateChangeMax=args.replRateChangeMax,
                 deathRateChangeMax=args.deathRateChangeMax,
                 hgtConjRateChangeMax=args.hgtConjRateChangeMax,
                 hgtTranRateChangeMax=args.hgtTranRateChangeMax,
                 motilityChangeMax=args.motilityChangeMax)\
                 -> Organism:
    """
    Give several attributes of an Organism a random chance to mutate.

    Parameters
    ----------
    org : Organism
        Organism whose stats are mutating.
    mutationRateChangeMax : float, optional
        Maximum by which `org`'s mutationRate attribute can change.
    replRateChangeMax : float, optional
        Maximum by which `org`'s replRate attribute can change.
    deathRateChangeMax : float, optional
        Maximum by which `org`'s deathRate attribute can change.
    hgtConjRateChangeMax : float, optional
        Maximum by which `org`'s hgtConjRate attribute can change.
    hgtTranRateChangeMax : float, optional
        Maximum by which `org`'s hgtTranRate attribute can change.
    motilityChangeMax : int, optional
        Maximum by which `org`'s motility attribute can change.

    Returns
    -------
    org : Organism
        Organism whose stats were mutated.

    Notes
    -----
    Probability of mutation of a given attribute is determined by the
    `mutationRate` attribute of the `Organism`.
    `MutationRate` attribute is given the chance to mutate last so that other
    random chances are not affected.

    """
    def mutate_stat(stat, statChangeMax) -> float:
        """
        Give a number a random chance to mutate by a random amount.

        Parameters
        ----------
        stat : int or float
            The number mutating.
        statChangeMax : int or float
            Maximum by which `stat` can change, either up or down.

        Returns
        -------
        stat : float
            The (potentially) new value of `stat`

        """
        if random.random() <= org.mutationRate:
            # increment by random float between pos/neg `statChangeMax`
            stat += random.uniform(-statChangeMax, statChangeMax)
        return stat

    def set_unit_interval(stat) -> float:
        """
        float between 0 and 1.

        Parameters
        ----------
        stat : float
            The value whose range is being checked and altered.

        Returns
        -------
        stat : float
            The input `stat` set to be between 0 and 1 if not already so.

        """
        if stat >= 1:
            stat = 1
        elif stat <= 0:
            stat = 0
        return stat

    # mutate several attributes of `org`
    # round() ensures motility stays an integer
    # set_unit_interval() ensures other stats stay between 0 and 1.
    org.motility = round(mutate_stat(org.motility, motilityChangeMax))
    org.replRate = set_unit_interval(
        mutate_stat(org.replRate, replRateChangeMax))
    org.deathRate = set_unit_interval(
        mutate_stat(org.deathRate, deathRateChangeMax))
    org.hgtConjRate = set_unit_interval(
        mutate_stat(org.hgtConjRate, hgtConjRateChangeMax))
    org.hgtTranRate = set_unit_interval(
        mutate_stat(org.hgtTranRate, hgtTranRateChangeMax))
    org.mutationRate = set_unit_interval(
        mutate_stat(org.mutationRate, mutationRateChangeMax))
    return org


def orgs_replicate_and_die(time,
                           idNumbSetter,
                           fedReplBoost=args.fedReplBoost,
                           hungerTolerance=args.hungerTolerance)\
                           -> tuple[int, int, int]:
    """
    Give each Organism a chance to replicate and die.

    Parameters
    ----------
    time : int
        The value to assign to a new Organism's birthday attribute,
        representing the timestep in which an Organism is instantiated.
    idNumbSetter : int
        The value to assign to a new Organism's idNumb attribute, as well as
        the key used to access said Organism in aliveOrgs dictionary.
    fedReplBoost : float, optional
        The temporary bonus added to an Organism's replRate when its
        hunger attribute equals 0, affecting probability of that Organism
        replicating when it is not hungry.
    hungerTolerance : int, optional
        The maximum an Organism's hunger attribute can reach without that
        Organism dying.

    Returns
    -------
    tuple[int, int, int]
        A tuple containing:
        - replCount (int): The number of Organism replications.
        - deathCount (int): The number of Organism deaths.
        - idNumbSetter (int): The new `idNumbSetter` input, after being
                              incremented each time an Organism was
                              instantiated by replication.

    Notes
    -----
    A list of Organisms from aliveOrgs dictionary is initially made so that
    only the Organisms that existed before this function was called can
    replicate and die. Newborn Organisms are put into the original aliveOrgs
    dictionary and not the aforementioned list, meaning they can neither
    replicate nor die in the current call of this function. This keeps the
    population more stable.

    """
    aliveOrgsLst = list(aliveOrgs.values())

    def orgs_replicate(idNumbSetter) -> tuple[int, int]:
        """
        Give each Organism a chance to replicate.

        Parameters
        ----------
        idNumbSetter : int
            The value to assign to a new Organism's idNumb attribute, as well
            as the key used to access said Organism in aliveOrgs dictionary.

        Returns
        -------
         tuple[int, int]
            A tuple containing:
            - replCount (int): The number of Organism replications.
            - idNumbSetter (int): The new `idNumbSetter` input, after being
                                  incremented each time an Organism was
                                  instantiated by replication.

        """
        replCount = 0

        for org in aliveOrgsLst:
            # if org's hunger attribute is 0, boost probability of replication.
            # otherwise there's effectively no boost.
            replBoost = fedReplBoost if not org.hunger else 0
            if random.random() <= (org.replRate + replBoost):
                # instantiate new Organism, many of whose attributes are
                # derived from org. mutate its stats. add it to aliveOrgs.
                aliveOrgs[idNumbSetter] = mutate_stats(
                    Organism(org.replicate_genome(),  # mutated genome
                             org.mutationRate,
                             org.genNumb + 1,
                             org.xCoord,
                             org.motility,
                             time,  # birthday
                             idNumbSetter,  # idNumb
                             org.replRate,
                             org.deathRate,
                             org.hgtConjRate,
                             org.hgtTranRate))

                idNumbSetter += 1
                replCount += 1
        return replCount, idNumbSetter

    def orgs_die() -> int:
        """
        Give each Organism a chance to die.


        Returns
        -------
        deathCount : int
            The number of Organism deaths.

        """
        deathCount = 0
        for org in aliveOrgsLst:
            # if org's hunger goes above hungerTolerance, always die regardless
            # of deathRate attribute
            if org.hunger > hungerTolerance \
                    or random.random() <= org.deathRate:
                # add genome to deadSeqsSpace
                deadSeqsSpace[org.xCoord - 1].append(org.genome)

                # move org from aliveOrgs to deadOrgs
                deadOrgs[org.idNumb] = aliveOrgs.pop(org.idNumb, 0)

                deathCount += 1
        return deathCount

    replCount, idNumbSetter = orgs_replicate(idNumbSetter)
    deathCount = orgs_die()
    return replCount, deathCount, idNumbSetter


def orgs_move(spaceSize=args.spaceSize):
    """
    Change each Organism's xCoord attribute.

    Parameters
    ----------
    spaceSize : int, optional
        The maximum an Organism's xCoord attribute can be set to.
        Globally, the size of the 1-dimensional space of the simulation.

    """
    for org in aliveOrgs.values():
        # move org by changing xCoord attribute
        org.xCoord += random.randint(-org.motility, org.motility)

        # keep xCoord attribute in bounds
        if org.xCoord < 1:
            org.xCoord = 1
        if org.xCoord > spaceSize:
            org.xCoord = spaceSize


def orgs_hgt_conj(hgtMaxDistance=args.hgtMaxDistance) -> int:
    """
    Give each Organism a chance to undergo horizontal gene transfer via
    conjugation.

    Extended Summary
    ----------------
    Give each Organism the chance to donate a gene to a nearby recipient,
    resulting in both having identical corresponding genes, and the recipient's
    original gene being permanently removed.

    Parameters
    ----------
    hgtMaxDistance: int, optional
        The maximum distance that can separate two Organisms while still
        allowing for HGT conjugation.

    Returns
    -------
    hgtConjCount : int
        The number of conjugations.

    """
    hgtConjCount = 0

    # Create a list of Organisms and sort by xCoord attribute.
    # Shuffle the list before sorting to ensure the order of Organisms with
    # equal xCoord values is not influenced by their insertion order in the
    # original dictionary.
    aliveOrgsLst = list(aliveOrgs.values())
    random.shuffle(aliveOrgsLst)
    aliveOrgsLst.sort(key=lambda organism: organism.xCoord)

    for i, donor in enumerate(aliveOrgsLst):
        if random.random() <= donor.hgtConjRate:

            recipientCandidates = []

            # If the donor is not the first in the list and is within
            # `hgtMaxDistance` of the previous Organism, add the previous
            # Organism to the list of potential recipients.
            if recipientCandidates != aliveOrgsLst[0]:
                if abs(donor.xCoord - aliveOrgsLst[i - 1].xCoord) \
                        <= hgtMaxDistance:
                    recipientCandidates.append(aliveOrgsLst[i - 1])

            # If the donor is not the last in the list and is within
            # `hgtMaxDistance` of the next Organism, add the next Organism
            # to the list of potential recipients.
            if donor != aliveOrgsLst[-1]:
                if abs(donor.xCoord - aliveOrgsLst[i + 1].xCoord) \
                        <= hgtMaxDistance:
                    recipientCandidates.append(aliveOrgsLst[i + 1])

            # If there are eligible recipients, select one at random.
            # Replace a randomly chosen gene in the recipient's genome with a
            # copy of the corresponding gene from the donor's genome.
            if recipientCandidates:
                recipient = random.choice(recipientCandidates)
                randIndex = random.randint(0, len(donor.genome) - 1)
                recipient.genome[randIndex] = copy.deepcopy(
                    donor.genome[randIndex])

                hgtConjCount += 1

    return hgtConjCount


def orgs_hgt_tran() -> int:
    """
    Give each Organism a chance to undergo horizontal gene transfer via
    transformation.

    Extended Summary
    ----------------
    Give each Organism the chanve to receive a gene from the genome of a dead
    Organism at the same location, replacing the corresponding gene in the
    recipient's genome and permanently removing the original gene.

    Returns
    -------
    hgtTranCount : int
        The number of transformations.
    """
    hgtTranCount = 0
    for recipient in aliveOrgs.values():

        # if recipient is going to hgt_tran
        # and there is a deadSeq in deadSeqSpace at recipient's xCoord
        if random.random() <= recipient.hgtTranRate:
            if deadSeqsSpace[recipient.xCoord - 1]:
                # pick random index in genome
                # (with which to access a random gene)
                randIndex = random.randint(0, len(recipient.genome) - 1)
                # pick random genome from deadSeqSpace that is at recipient's
                # xCoord, and remove it from deadSeqSpace
                randDeadGenome = copy.deepcopy(
                    deadSeqsSpace[recipient.xCoord - 1].pop(
                        random.randint(0, len(
                            deadSeqsSpace[recipient.xCoord - 1]) - 1)))
                # replace recipient's gene with randomDeadGenome
                recipient.genome[randIndex] = copy.deepcopy(
                    randDeadGenome[randIndex])

                del randDeadGenome
                del randIndex
                hgtTranCount += 1
    return hgtTranCount


def orgs_eat(foodSpace) -> tuple[int, int]:
    """
    Give Organisms a chance to eat Food, affecting their hunger attribute.

    Extended Summary
    ----------------
    Each Food object at a given xCoord can only feed one Organism at that
    xCoord per call of this function. If an Organism eats, its hunger
    attribute is set to 0. Otherwise, its hunger increases by 1.

    Parameters
    ----------
    foodSpace : list[int]
        A list where the index of each element represents a location in
        1-dimensional space, and each element is how many Foods exist at that
        location.

    Returns
    -------
    tuple[int, int]
        A tuple containing:
            - eatCount (int): The number of Organisms who ate.
            - notEatCount (int): The number of Organisms who did not eat.

    Notes
    -----
    No food objects are involved in this function. Only `foodSpace`, a list of
    integers representing amounts of Food at specific locations, is being
    altered when an Organism eats. This means that globally, Food objects can
    sate an Organism once per timestep for multiple timesteps, moving in space
    each time, until they eventually reset and go back to their spawnPoints.

    """
    # make list of Organisms and shuffle to ensure no bias of which Organisms
    # eat first
    foodQueue = list(aliveOrgs.values())
    random.shuffle(foodQueue)

    eatCount = 0
    notEatCount = 0

    for org in foodQueue:
        # if food is available at org's xCoord, org eats
        if foodSpace[org.xCoord - 1]:
            # decrease food available at that xCoord for Organisms still
            # waiting in the foodQueue. set org's hunger to 0, as it is no
            # longer hungry.
            foodSpace[org.xCoord - 1] -= 1
            eatCount += 1
            org.hunger = 0

        # if no food available at org's xCoord, org does not eat.
        # increase org's hunger by 1
        else:
            notEatCount += 1
            org.hunger += 1

    return eatCount, notEatCount


def orgs_spont_mutate():
    """
    Give every Organism a chance to replace their genome with a mutated copy.
    """
    for org in aliveOrgs.values():
        if random.random() <= org.spontMutationRate:
            org.genome = org.replicate_genome()


def dct_to_space(dct, spaceSize=args.spaceSize) -> list[int]:
    """
    Make a list of how many objects of a given class are at each location in
    space.

    Parameters
    ----------
    dct : dict[int, Food or Organism]
        A dictionary of objects with an attribute `xCoord` stored as values.
    spaceSize : int, optional
        The maximum value that the `xCoord` attributes of objects can equal.

    Notes
    -----
    space : list[int]
        A list where the index of each element represents a location in
        1-dimensional space, and each element is how many objects exist at that
        location.

    """
    # count how many objects in `dct` have a given `xCoord` attribute,
    # for every possible value `xCoord` can equal
    xCoords = [obj.xCoord for obj in dct.values()]
    xCoordCounts = Counter(xCoords)
    space = [xCoordCounts[i] for i in range(1, spaceSize + 1)]
    del xCoords, xCoordCounts
    return space


def rand_genome(geneNumber, geneLength) -> list[list[str]]:
    """
    Generate a genome where each gene consists of randomly generated
    nucleotides.

    Parameters
    ----------
    geneNumber : int
        The number of genes in genome.
    geneLength : int
        The number of nucleotides in genes.

    Returns
    -------
    list[list[str]]
        A randomly generated genome with each nested element representing one
        of four nucleotides: 'A', 'C', 'G', or 'T'.

    Examples
    --------
    # >>>rand_genome(2,3)
    [['C', 'C', 'A'], ['A', 'T', 'C']]  # random

    """
    genome = []
    for i in range(geneNumber):
        genome.append([random.choice("ACGT") for _ in range(geneLength)])

    return genome


def display_dead_seqs_space():
    """
    Display `deadSeqSpace` vertically, where all the genomes in a given
    location are printed on a new line. This is easier to look at then printing
    the list directly.
    """
    print('deadSeqsSpace:')
    for seqLst in deadSeqsSpace:
        print(seqLst)


def census_csv(time,
               runID=args.runID,
               outputDirectory=args.outputDirectory):
    """
    Write a csv file containing the DNA sequences and some other attributes
    of all living Organisms.

    Parameters
    ----------
    time : int
        Used to name output files, denoting from which timestep the data was
        gathered.
    runID : str, optional
        Used to name output folders and files within the scope of
        `outputDirectory`, denoting from which run of the simulation the data
        was gathered.
    outputDirectory : str, optional
        Directory pathway in which data across all runs of the simulation are
        designed to be stored.

    """
    # start writing string that will become file by adding column for `idNumb`
    # attribute.
    seqCSV = 'idNumb,'

    # Look in `aliveOrgs` and `deadOrgs` for Organism whose dictionary key is
    # 0 (arbitrary number). Then add columns to csv for each gene in `genome`.
    try:
        for i, _ in enumerate(deadOrgs[0].genome):
            seqCSV += f'gene{i},'
    except KeyError:
        for i, _ in enumerate(aliveOrgs[0].genome):
            seqCSV += f'gene{i},'

    # add columns to csv for `xCoord` and `genNumb` attributes
    seqCSV += 'xCoord,genNumb\n'

    # start filling in Organism data under the labeled columns for every
    # living Organism
    for org in aliveOrgs:
        # add org's `idNumb`
        seqCSV += f'{aliveOrgs[org].idNumb},'

        # add org's gene(s) in `genome`
        for gene in aliveOrgs[org].genome:
            for nucleotide in gene:
                seqCSV += nucleotide
            seqCSV += ','

        # add org's `xCoord` and `genNumb`
        seqCSV += f"{aliveOrgs[org].xCoord},{aliveOrgs[org].genNumb}\n"

    # make the branch directory to which the csv file will be added
    outputDirectory = os.path.join(outputDirectory, runID,
                                   f"{runID}_sequenceData")
    os.makedirs(outputDirectory, exist_ok=True)

    # open, write, and (automatically) close the csv file
    filePath = os.path.join(outputDirectory, f'{runID}_timeStep{time}.csv')
    with open(filePath, 'w') as f:
        f.write(seqCSV)


def args_csv():
    """
    write a csv file of all the parameters and their inputted arguments from
    argparse.
    """

    # `args` is not iterable, but rather a `Namespace` object with attributes.
    # Make a dictionary (iterable) where dictionary keys are attribute names
    # and dictionary values are strings of attribute values.
    # (e.g., `argsDct['hungerTolerance']` == '1')
    argsDct = {}
    for arg in vars(args):
        attrValue = str(getattr(args, arg))
        # replace commas with vertical bar to not break formatting of csv file
        attrValue = attrValue.replace(',', '|')
        argsDct[arg] = attrValue

    # format string that will be turned into csv file, adding columns for
    # parameter and user-inputted argument for that parameter.
    argsCSV = 'parameter, argument\n'
    for arg in argsDct:
        argsCSV += f'{arg},{argsDct[arg]}\n'

    # make the branch directory to which the csv file will be added
    outputDirectory = os.path.join(argsDct['outputDirectory'],
                                   argsDct['runID'])
    os.makedirs(outputDirectory, exist_ok=True)

    # open, write, and (automatically) close the file
    filePath = os.path.join(outputDirectory, f'{argsDct["runID"]}_args.csv')
    with open(filePath, 'w') as f:
        f.write(argsCSV)


def plot_orgs_and_food(orgSpace, foodSpace, time,
                       spaceSize=args.spaceSize,
                       yLim=args.plotYLim):
    """
    Plot the number of Organism and Food objects over space in a given
    timestep.

    Parameters
    ----------
    orgSpace : list[int]
        A list where the index of each element represents a location in
        1-dimensional space, and each element is how many `Organism` objects
        exist with that `xCoord` attribute.
    foodSpace : list[int]
        A list where the index of each element represents a location in
        1-dimensional space, and each element is how many `Food` objects
        exist with that `xCoord` attribute.
    time : int
        Part of plot title, denoting current timestep.
    spaceSize : int, optional
            Size of the x-axis. Globally, the size of the 1-dimensional space
            of the simulation.
    yLim : int or None, optional
        If not `None`, the size of the y-axis.

    """
    # plot Food and Organisms
    plt.plot([_ for _ in range(spaceSize)], foodSpace, marker='o',
             linestyle='-', color='b', label='Food')
    plt.bar([_ for _ in range(spaceSize)], orgSpace, color='r',
            label='Organisms')

    # set x- and y-axis limits
    if yLim is not None:
        plt.ylim(0, yLim)
    plt.xlim(0, spaceSize)

    # labels and title
    plt.xlabel('Space')
    plt.ylabel('Number of Objects')
    plt.title(f'Number of Organisms and Food Across Space\n '
              f'Time = {time}\n '
              f'Population = {len(aliveOrgs)}')
    plt.legend(loc="upper left", frameon=False)

    # pause before clearing plot
    plt.pause(0.01)
    plt.clf()


def run_sim(censusTimes=args.censusTimes):
    """
    Simulate a population of `Organism` objects that interact with each other
    and their environment and undergo biological processes, for several
    timesteps.

    Extended Summary
    ----------------
    The following are the various ecological interactions and biological
    processes occurring in each timestep:
        (1) `Food`s diffuse through space, affecting their `xCoord` attributes.
        (2) Organisms might spontaneously mutate their `genome` attribute.
        (3) Organisms might replicate, instantiating new Organisms.
        (4) Organisms might die.
        (5) Organisms might donate genes via conjugation, affecting another's
            `genome` attribute.
        (6) Organisms might receive genes via transformation, affecting their
            `genome` attribute.
        (7) Organisms might eat, affecting their `hunger` attribute.
        (8) Organisms move, affecting their `xCoord` attribute.

    Parameters
    ----------
    censusTimes : tuple[int, ...], optional
        the timesteps when the census_csv() function is called. Affects when
        data are written to files during the simulation, as well as how many
        timesteps are simulated.

    Notes
    -----
    `time` starts at 1 instead of 0. This is to distinguish between
    instances of `Organism` in the parent generation and those instantiated
    in the 1st timestep. These two groups have a different `birthday`
    attribute: 0 and 1, respectively.

    """
    # used to make new `Organism` objects. Acts as key in dictionary, and as
    # their `idNumb` attribute upon instantiation. Unique to each `Organism`
    idNumbSetter = len(aliveOrgs)

    for time in range(1, censusTimes[-1] + 1):

        # stop simulation if population went extinct.
        if not aliveOrgs:
            print("population has gone extinct.")
            break

        print(f'\n* * * * * * * * * * *'
              f'\n{time = }')

        # diffuse the `Food` objects through the space via random walk
        food_diffusion()

        # plot the `Food` and `Organism` objects in space with matplotlib
        foodSpace = dct_to_space(foods)
        orgSpace = dct_to_space(aliveOrgs)

        plot_orgs_and_food(orgSpace, foodSpace, time)

        # `Organism` objects undergo various actions and biological processes
        orgs_spont_mutate()
        replCount, deathCount, idNumbSetter = \
            orgs_replicate_and_die(time, idNumbSetter)
        hgtConjCount = orgs_hgt_conj()
        hgtTranCount = orgs_hgt_tran()
        eatCount, notEatCount = orgs_eat(foodSpace)
        orgs_move()

        del orgSpace, foodSpace

        # print some updates on what has happened in this timestep.
        print(f"\npopulation size = {len(aliveOrgs)}"
              f"\n{replCount = }"
              f"\n{deathCount = }"
              f"\n{hgtConjCount = }"
              f"\n{hgtTranCount = }"
              f"\n{eatCount = }"
              f"\n{notEatCount = }")

        if time in censusTimes:
            print('\n----------------------------------'
                  '\n- TAKING CENSUS; OUTPUTTING FILE -'
                  '\n----------------------------------')
            census_csv(time)

    print(f'\n ---------------------'
          f'\n - END OF SIMULATION -'
          f'\n ---------------------'
          f'\n{len(aliveOrgs) = }'
          f'\n{len(deadOrgs) = }')


# Initialize `deadSeqSpace`. This is where `genome` attributes of dead
# `Organism` objects are stored that have not yet been uptaken by alive
# `Organism` objects via transformation.
# It will contain nucleotides in genes in genomes at locations

start = time.time()

deadSeqsSpace: list[list[list[list[str]]]] = \
    [[] for _ in range(args.spaceSize)]

# make food
foods = initialize_food(args.foodSources)

# make parent generation
aliveOrgs = initialize_p_gen(genomeType='eColi16s', geneNumber=3, geneLength=4)
deadOrgs = {}

# run the sim
run_sim()

# output a file of all the parameters from argparse
args_csv()

end = time.time()

print(f"duration is {end - start} seconds")
