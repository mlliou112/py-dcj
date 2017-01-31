import collections
import itertools


class AbstractMarker:
    """
    General Representation for a marker in DCJ model.
    """
    # TODO Should look into abstract/metaclasses in python
    def __init__(self):
        pass


class Telomere(AbstractMarker):
    """
    Representation for a telomere in DCJ model.
    """
    def __init__(self):
        self.reverse = False
        pass

    def __repr__(self):
        return "."

    def __eq__(self, other):
        return type(self) == type(other)

    def __hash__(self):
        return hash(".")


class Marker(AbstractMarker):
    """
    Representation for DNA fragment in DCJ model.

    Parameters
    ----------
    uid: int (required)
        unique id

    dna: str
        What dna fragment does this marker represent?

    reverse: boolean
        Is the DNA fragment reversed?
    """
    # Requirements:
    # Reversed or not
    # distinguishable and unique
    # telomere
    # Hash defined as uid, because that shouldn't change

    def __init__(self, uid, dna="", reverse=False):
        self.reverse = reverse
        self.dna = dna
        self.uid = uid

    def __eq__(self, other):
        if type(other) is type(self):
            return self.uid == other.uid
        return False

    def __repr__(self):
        return "-" + str(self.uid) if self.reverse else str(self.uid)

    def __hash__(self):
        return hash(self.uid)


class MarkerEnd():
    """Markers with designated head or tail end. Used in adjacencies"""
    def __init__(self, marker, head=True):
        self.marker = marker
        self.head = head

    def __repr__(self):
        return str(self.marker) if self.head else str(self.marker) + "*"

    def __eq__(self, other):
        return self.marker == other.marker and self.head == other.head

    def other_adjacency_end(self, adj):
        other_end_marker = adj.right_end_marker if adj.left_end_marker == self else adj.left_end_marker
        return other_end_marker

    def __hash__(self):
        return hash(self.marker) * hash(self.head)


class Chromosome(collections.UserList):
    """The abstraction of a chromosome in DCJ model, a list of telomeres and markers

    Parameters
    ----------
    markers: list, file, or string (req)
        The object to be converted into a chromosome. Currently only supports pilon.changes file.

    Fields
    -------
    content: list
        The list of markers that represent the chromosome

    """
    # Just a list of markers
    #TODO Validate that chromosome has only 0 or 2 telomeres

    string_d = {}
    uid_counter = 0

    def __init__(self, markers):
        # If it's already a list of markers, skip all this. Otherwise.
        # FIXME Skipping code based on if the first element is a marker is really hackish
        if isinstance(markers, list):
            content = markers
        else:
        # Otherwise loop through like it is a string
            content = []
            for s in markers:
                if s in self.string_d:
                    i = self.string_d[s]
                else:
                    self.string_d[s] = self.uid_counter
                    i = self.uid_counter
                    self.uid_counter += 1
                dna = s
                rev = s.isupper()
                telomere = s == "."
                if telomere:
                    content.append(Telomere())
                else:
                    content.append(Marker(i, dna, rev))
        super().__init__(content)




class Genome(collections.UserList):
    """The abstraction of a genome in DCJ model, a set of chromosomes"""
    # List of Chromosomes
    def __init__(self, *args):
        super().__init__(args)

    # def __repr__(self):
    #     return str(self.content)


class Adjacency(collections.UserList):
    # Requirements
    # label
    # TODO Should inherit from something tuple? for now list.
    def __init__(self, left_end_marker=None, right_end_marker=None, label=None):
        # List of markers, chromosome.
        self.label = label
        self.left_end_marker = left_end_marker
        self.right_end_marker = right_end_marker
        super().__init__([left_end_marker, right_end_marker])

    def __repr__(self):
        return "[" + str(self.left_end_marker) + ", " + str(self.label) + ", " + str(self.right_end_marker) + "]"


class AdjacencyGraph:
    """Adjacency information for a genome"""

    # Find common elements
    # Make Adjacency lists
    # Make two tables for index
    # Form Bipartite graph

    def __init__(self, A, B):
        # PseudoCode
        # Find A, B, G marker set
        # Create Adjacency lists

        # Finding Common Elements
        #FIXME: . does not belong in the unique marker set, should be in common?
        marker_set_a = set(itertools.chain.from_iterable(A))
        marker_set_b = set(itertools.chain.from_iterable(B))
        self.commonMarkers = marker_set_a.intersection(marker_set_b)
        marker_set_a = marker_set_a.difference(self.commonMarkers)
        marker_set_b = marker_set_b.difference(self.commonMarkers)

        # Make Adjacency lists
        self.adjA = []
        self.adjB = []
        adj = [None, None]
        # Label is a list of markers NOT endmarkers

        # Created Adjacencies for Genomes A and B
        adjacencies = [self.adjA, self.adjB]
        reference_A = {MarkerEnd(Telomere(), head=True): None}
        reference_B = {MarkerEnd(Telomere(), head=True): None}
        references = [reference_A, reference_B]
        for i, genome in enumerate([A, B]):
            for chromosome in genome:
                index = 0
                adj_index = 0
                current_marker = chromosome[index]
                while index < len(chromosome) - 1:
                    label = []

                    # Fill first adjancy postion
                    if current_marker == Telomere():
                        adj[0] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[0] = MarkerEnd(current_marker, head=not current_marker.reverse)

                    # While loop to add to label
                    next_marker = chromosome[index + 1]
                    while next_marker not in self.commonMarkers:
                        label.append(next_marker)
                        index += 1
                        next_marker = chromosome[index + 1]

                    # fill second adjacency position
                    if next_marker == Telomere():
                        adj[1] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[1] = MarkerEnd(next_marker, head=next_marker.reverse)
                    adjacencies[i].append(Adjacency(left_end_marker=adj[0], right_end_marker=adj[1], label=label))

                    # Create the reference for A and B
                    if adj[0].marker != Telomere():
                        references[i][adj[0]] = adj_index
                    if adj[1].marker != Telomere():
                        references[i][adj[1]] = adj_index

                    current_marker = next_marker
                    index += 1
                    adj_index += 1


        # Traverse the graph
        to_visitA_index = set(range(len(self.adjA)))
        to_visitB_index = set(range(len(self.adjB)))
        visitedA_index = set()
        visitedB_index = set()

        cycles = 0
        ab_paths = 0

        while to_visitA_index:
            a_side = True
            current_adj_index = to_visitA_index.pop()
            visitedA_index.add(current_adj_index)

            left_marker = self.adjA[current_adj_index].left_end_marker
            right_marker = self.adjA[current_adj_index].right_end_marker

            # Checks for AB paths
            paths_end_on_A = [True, True]
            for i, current_marker in enumerate([left_marker, right_marker]):
                next_adj_index = reference_A[current_marker]
                if next_adj_index in visitedB_index:
                    # If cycle on left_marker, no need to check right marker
                    continue
                while next_adj_index is not None:
                    current_adj_index = next_adj_index
                    a_side = not a_side

                    adj_side = self.adjA if a_side else self.adjB
                    reference_side = reference_A if a_side else reference_B

                    visitedA_index.add(current_adj_index) if a_side else visitedB_index.add(current_adj_index)
                    to_visitA_index.remove(current_adj_index) if a_side else to_visitB_index.remove(current_adj_index)


                    # Switch markers
                    if current_marker != adj_side[current_adj_index].left_end_marker:
                        current_marker = adj_side[current_adj_index].left_end_marker
                    else:
                        current_marker = adj_side[current_adj_index].right_end_marker


                    # Find next location
                    next_adj_index = reference_side[current_marker]


                    if next_adj_index is None:
                        # Path ends
                        paths_end_on_A[i] = a_side
                        break
                    elif (a_side and next_adj_index in visitedB_index) or (not a_side and next_adj_index in visitedA_index):
                        # cycle
                        cycles += 1
                        break
                    else:
                        # The path continues
                        pass

            # AB_path
            if paths_end_on_A[0] != paths_end_on_A[1]:
                ab_paths += 1

#
#
# def dcj_pilon_changes(f):
#     # Convention used for marker UID
#     #  common -> x%3 = 0
#     # Unique A -> x%3 = 1
#     # Unique B -> x%3 = 2
#     genome_A_indx = 1
#     genome_B_indx = 2
#     common_index = 3
#     genome_A = []
#     genome_B = []
#     genome_A += [Telomere(), Marker(uid=0)]
#     genome_B += [Telomere(), Marker(uid=0)]
#     with open(f, "r") as fi:
#         for line in fi:
#             ar = line.split()
#             if len(ar) != 4:
#                 raise ValueError("Unrecognized line in pilon.changes file")
#             if ar[2] != "." and ar[3] != ".":
#                 # If it's the reverse dna strand,
#                 if ar[2][::-1] == ar[3]:
#                     genome_A.append(Marker(uid=common_index))
#                     genome_B.append(Marker(uid=common_index, reverse=True))
#                     common_index += 3
#             if ar[2] != ".":
#                 genome_A.append(Marker(uid=genome_A_indx, dna=ar[2]))
#                 genome_A_indx += 3
#             if ar[3] != ".":
#                 genome_B.append(Marker(uid=genome_B_indx, dna=ar[3]))
#                 genome_B_indx += 3
#             genome_A.append(Marker(uid=common_index))
#             genome_B.append(Marker(uid=common_index))
#             common_index += 3
#     genome_A.append(Telomere())
#     genome_B.append(Telomere())
#
#     return Chromosome(genome_A), Chromosome(genome_B)



def calculate_distance(A, B):
    """Calculate the DCJ distance between genome A and B"""


if __name__ == '__main__':
    a = Chromosome(".acb.")
    b = Chromosome(".ab.")
    A = Genome(a)
    B = Genome(b)
    res = AdjacencyGraph(A, B)
    print(res.adjA)
    print(res.adjB)

