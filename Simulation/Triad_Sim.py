from Bio import Align
from collections import defaultdict

class Seq:
    """
    Holds a genome sequence and a list of every edge it is in
    """
    def __init__(self, seq):
        self.seq=seq
        # self.edge_set=set()
        # self.test=test

    # For print debugging, ignore
    def __str__(self):
        return self.seq
    
    # For debug representation, ignore
    def __repr__(self):
        return self.seq
    
    
    # # Not all sequences are unique so this is just a bad idea
    # def __eq__(self, value):
    #     return self.seq==value.seq
    
    # def __hash__(self):
    #     return hash(self.seq)

    def comp(self, other):
        """
        Compares genome sequences of a given sequence object to it's own genome sequence

        Returns
        -------
        Comparison value
        """
        align= Align.PairwiseAligner()
        return align.score(self.seq, other.seq)
    
class Edge:
    """
    Edge class that takes in two sequences and optionally can take in a value of the percent similarity of the two sequences.
    Edge gets added to each vertex's edgelist.
    """
    def __init__(self, s1: Seq, s2: Seq, value= 0.9):
        self.vertices= frozenset({s1, s2})
        # s1.edge_set.add(self)
        # s2.edge_set.add(self)
        self.value=value

    # For print debugging, ignore
    def __str__(self):
        return str(self.vertices)
    
    def __eq__(self, edge):
        return self.vertices==edge.vertices
    
    def __hash__(self):
        return hash(self.vertices)
        
    # For 'in' keyword
    def __contains__(self, seq: Seq):
        return seq in self.vertices
        

    def build_triad(self, seq: Seq, threshold):
        """
        Builds a Triad object from a given sequence onto this edge if all sequences meet the threshold similiarity 

        Returns
        -------
        Triad if threshold is met and None otherwise
        """
        # If seqs in the edge pass the threshold, add them to the contained list
        contained= []
        for s in self.vertices:
            if (seq.comp(s))>=threshold:
                contained.append(s)

        # If both sequences pass the threshold, build a Triad
        if len(contained)==2:
            return Triad(contained[0], contained[1], seq)
        else:
            return None

        # if (seq.comp(self.vertices[0])>=threshold) and (seq.comp(self.vertices[1])>=threshold):
        #     return Triad(self.vertices[0], self.vertices[1], seq)
        # else:
        #     return None
    
class Triad:
    """
    Triad of three sequences that pass a similarity threshold
    """
    def __init__(self, seq1, seq2, seq3):
        # Creates an unordered set of sequences so that Triad1.seqs==Triad2.seqs if they contain the same sequences that are not necessarily in the same order
        self.seqs= frozenset({seq1, seq2, seq3})

        # Creates an unordered set of edges
        self.edges= set((Edge(seq1, seq2),Edge(seq2,seq3),Edge(seq1, seq3)))

        # All similar triads in grouping
        self.group=set()

    # For 'in' keyword
    def __contains__ (self, edge: Edge):
        return edge in self.edges

    # For print debugging, ignore
    def __str__(self):
        return self.seqs
    
    # For debugger representation, ignore
    def __repr__(self):
        seqs = []
        for seq in self.seqs:
            seqs.append(seq)
        
        ret= f"Triad({seqs[0]}  {seqs[1]}  {seqs[2]})"

        return ret

def read_data() -> list[Seq]:
    # open file
    file_name= "Simulation\\testRun_modified_timeStep100.csv"
    file= open(file_name, "r")
    lines= file.readlines()
    file.close()

    seq_list= []
    
    for line in lines[1:]:
        # remove the newline character from the end
        line = line.strip()

        # turn the line into a list of strings
        line = line.split(",")

        # add the sequences to list of sequence objects
        seq= Seq(line[1])
        seq_list.append(seq)

    return seq_list

def dfs (triad_lst: list[Triad]):
    """ 
    Depth First Search to group Triads together
    """
    # Uses a dictionary that maps each edge to a list of indices of triads that contain it
    edge_triad_map= defaultdict(list)
    for i, triad in enumerate(triad_lst):
        for edge in triad.edges:
            edge_triad_map[edge].append(i)

    # Build adjacency lists for triads by adding all triads with a shared edge to each other's adjacency list
    adjacency_list= defaultdict(list)
    for edge, triads in edge_triad_map.items():
        for i in range(len(triads)):
            for ii in range(i+1, len(triads)):
                t1= triads[i]
                t2= triads[ii]
                adjacency_list[t1].append(t2)
                adjacency_list[t2].append(t1)

    # Visited triads
    visited = [False] * len(triad_lst)
    # Groups of triads
    groupings = []

    # Iterate over each triad using their index
    for i in range(len(triad_lst)):
        # If Triad is not visited -> create a new grouping by searching through until no further shared edges are found, else make a new group with next unvisited triad
        if not visited[i]:
            # Initialise the group of Triads
            group = []
            # create a stack exploration list that holds the indices of the triads 
            stack = [i]
            # sets the current triad as explored
            visited[i] = True
            
            while len(stack)>0:
                # add current to the groupings list and remove from the exploration queue
                current = stack.pop()
                group.append(triad_lst[current])

                # add all unvisited neighbors to the exploration list and marks them as visited
                for neighbour in adjacency_list[current]:
                    if not visited[neighbour]:
                        visited[neighbour] = True
                        stack.append(neighbour)
                    
            # Return the groupings
            groupings.append(group)
    
    return groupings

def test():
    """"
    Test sequence grouping
    """
    a=Seq("A")
    b=Seq("B")
    c=Seq("C")
    d=Seq("D")
    e=Seq("E")
    f=Seq("F")
    g=Seq("G")
    j=Seq("J")
    k=Seq("K")
    l=Seq("L")
    x=Seq("X")
    y=Seq("Y")
    z=Seq("Z")
    
    triad_list = [Triad(a,b,c),Triad(a,b,d),Triad(b,c,d),Triad(d,e,f),Triad(e,f,g),Triad(x,y,z),Triad(a,b,e),Triad(j,k,l),Triad(b,c,j)]
    groupings = dfs(triad_list)

    print(groupings)

def main():
    test()

    seq_list= read_data()
    threshold= 1500.0
    edge_set= set()

    
    # Create all edges that meet threshold
    for i in range(len(seq_list)):
        # print(i)
        seqA = seq_list[i]
        for ii in range(i+1,len(seq_list)):
            seqB= seq_list[ii]
            if (seqA.comp(seqB))>=threshold:
                edge_set.add(Edge(seqA,seqB))
            
    print(len(edge_set))
    triad_lst = []
    edge_triad_map = defaultdict(list)

    # Create triads for every edge that meets the threshold
    for edge in edge_set:
        for seq in seq_list:
            triad = edge.build_triad(seq, threshold)
            if not (triad==None):
                triad_lst.append(triad)
                edge_triad_map[edge].append(triad)
    
    print(len(triad_lst))
    
    groupings = dfs(triad_lst)
    print(groupings)


if __name__=='__main__':
    main()
