"""
Baer Hodge
Jan 7, 2025
"""

import random
import copy


class Organism:
    """
    Object with genetic information and other attributes meant to simulate a
    pre-LUCA microbe.

    Attributes
    ----------
    genome : list[list[str]]
        One or more genes, each of whose elements are strings of
        'A', 'C', 'G' or 'T', representing nucleotides.
    mutationRate : float
        Probability of various attributes of a newly instantiated child
        Organism to be different from those of its parent. Also, the
        probability of each individual nucleotide in an Organism's genome to
        get copied incorrectly in the replicate_genome() method.
    genNumb : int
        How many generations an Organism is removed from the parent generation.
         This value is always 1 greater than that of the Organism's parent.
    xCoord : int
        Location in 1-Dimensional space of the simulation.
    motility : int
        Maximum by which an Organism's xCoord can change in a given timestep.
    birthday : int
        Timestep in which an Organism is instantiated.
    idNumb : int
        Unique number for every Organism. Also acts as a key in the
        dictionaries in which the Organisms are stored.
    replRate : float
        Probability of an Organism giving rise to a new Organism in a given
        timestep, all of whose attributes are derived from the parent.
    deathRate : float
        Probability of an Organism dying in a given timestep, and leaving its
        genome loose in the environment at the location of death.
    hgtConjRate : float
        Probability of an Organism undergoing horizontal gene transfer (HGT)
        via conjugation in a given timestep, if criteria are met. This involves
        the Organism being the donor of genetic material to another nearby
        Organism, who acts as the recipient of genetic material.
    hgtTranRate : float
        Probability of an Organism undergoing horizontal gene transfer (HGT)
        via transformation in a given timestep, if criteria are met. This
        involves the Organism being the recipient of genetic material from a
        dead Organism who has left its genome loose in the environment.
    spontMutationRate : float
        Probability that an Organism's genome will spontaneously mutate in a
        given timestep and be replaced by self.replicate_genome().
        Always one tenth (1/10) of replRate, so that 10% of mutations come from
        spontaneous mutations while 90% come from vertical gene transfer.
    hunger : int
        How hungry an Organism is. Organism dies if hunger gets too high.

    Methods
    -------
    replicate_genome()
        returns a mutated deepcopy of self's genome.
    """

    def __init__(self,
                 genome,
                 mutationRate,
                 genNumb,
                 xCoord,
                 motility,
                 birthday,
                 idNumb,
                 replRate,
                 deathRate,
                 hgtConjRate,
                 hgtTranRate):

        self.genome = genome
        self.mutationRate = mutationRate
        self.genNumb = genNumb
        self.xCoord = xCoord
        self.motility = motility
        self.birthday = birthday
        self.idNumb = idNumb
        self.replRate = replRate
        self.deathRate = deathRate
        self.hgtConjRate = hgtConjRate
        self.hgtTranRate = hgtTranRate
        self.spontMutationRate = replRate / 10
        self.hunger = 0

    def replicate_genome(self) -> list[list[str]]:
        """
        Simulate DNA replication by iterating through a deepcopy of the genome
        and giving each nucleotide a chance to mutate. Likelihood of mutation
        determined by mutationRate.

        Returns
        -------
        genomeReplica : list[list[str]]
            The mutated copy of the genome.

        """
        # create deepcopy of nested list
        genomeReplica = copy.deepcopy(self.genome)

        # iterate through each gene in genomeReplica
        for gene in genomeReplica:
            for i, nucleotide in enumerate(gene):

                # give each nucleotide of each gene a chance to mutate
                # to a nucleotide other than the one it already is
                if random.random() <= self.mutationRate:
                    gene[i] = random.choice("ACGT".replace(gene[i], ''))

        return genomeReplica
