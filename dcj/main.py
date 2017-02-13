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


class MarkerEnd:
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
        self.marker_set = set()
        # If it's already a list of markers, skip all this. Otherwise.
        # FIXME Skipping code based on if the first element is a marker is really hackish
        if isinstance(markers, list):
            content = markers
        else:
        # Otherwise loop through like it is a string
            content = []
            # If first is telomere, last must be telomere or vice-versa
            if not (markers[0] == ".") == (markers[-1] == "."):
                raise ValueError("Linear Chromosome must start and end with telomeres.")
            for i, s in enumerate(markers):
                dna = s
                rev = s.isupper()
                telomere = s == "."
                if s.lower() in self.marker_set and not telomere:
                    raise ValueError("Duplicated markers are not allowed. ({0})".format(s.lower()))
                elif s.lower() in Chromosome.string_d:
                    uid = Chromosome.string_d[s.lower()]
                    self.marker_set.add(s.lower())
                else:
                    Chromosome.string_d[s.lower()] = Chromosome.uid_counter
                    self.marker_set.add(s.lower())
                    uid = Chromosome.uid_counter
                    Chromosome.uid_counter += 1
                if telomere and i not in [0, len(markers) - 1]:
                    raise ValueError("Telomere cannot appear in middle of chromosome.")
                elif telomere:
                    content.append(Telomere())
                else:
                    content.append(Marker(uid, dna, rev))
        super().__init__(content)

    @classmethod
    def clear(cls):
        """
        Resets the string_d and uid_counter.
        :return:
        """

        cls.string_d = {}
        cls.uid_counter = 0


class Genome(collections.UserList):
    """The abstraction of a genome in DCJ model, a set of chromosomes"""
    # List of Chromosomes
    def __init__(self, *args):
        for c in args:
            if not isinstance(c, Chromosome):
                raise TypeError("Argument list passed to Genome constructor must be Chromosome class instance")

        if len(args) > 1:
            duplicated_markers = set.intersection(*[c.marker_set for c in args])
        else:
            duplicated_markers = None

        if duplicated_markers:
            raise ValueError("Duplicated markers found across chromosomes. ({0})".format(duplicated_markers))
        super().__init__(args)

    # Must retain information about markers between chromosomes
    # Validate that arguments passed in are in fact chromosomes
    # Validate that the markerset in each chromo
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
        marker_set_a = set(itertools.chain.from_iterable(A))
        marker_set_b = set(itertools.chain.from_iterable(B))
        self.commonMarkers = marker_set_a.intersection(marker_set_b)
        if Telomere() in self.commonMarkers:
            self.commonMarkers.remove(Telomere())
        marker_set_a = marker_set_a.difference(self.commonMarkers)
        marker_set_b = marker_set_b.difference(self.commonMarkers)

        # Make Adjacency lists
        self.adjA = []
        self.adjB = []
        adj = [None, None]
        # Label is a list of markers NOT endmarkers

        # Created Adjacencies for Genomes A and B
        # FIXME: Does not work for circular genomes
        adjacencies = [self.adjA, self.adjB]
        reference_A = {MarkerEnd(Telomere(), head=True): None}
        reference_B = {MarkerEnd(Telomere(), head=True): None}
        references = [reference_A, reference_B]
        for i, genome in enumerate([A, B]):
            adjacency_length = 0
            for chromosome in genome:
                index = 0
                adj_index = 0
                current_marker = chromosome[index]
                chromosome_markers = set([x for x in chromosome])
                adjacency_length += len(self.commonMarkers.intersection(chromosome_markers))
                if chromosome[0] == Telomere():
                    adjacency_length += 1
                while len(adjacencies[i]) < adjacency_length:
                    label = []

                    # Fill first adjacency position
                    if current_marker == Telomere():
                        adj[0] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[0] = MarkerEnd(current_marker, head=not current_marker.reverse)

                    # While loop to add to label
                    next_marker = chromosome[0] if index >= len(chromosome) - 1 else chromosome[index+1]
                    while next_marker not in self.commonMarkers:
                        # Fixes linear chromosomes
                        if next_marker == Telomere():
                            break
                        label.append(next_marker)
                        index += 1
                        # Circular genome, must wrap around
                        next_marker = chromosome[0] if index >= len(chromosome) - 1 else chromosome[index + 1]

                    # fill second adjacency position
                    if next_marker == Telomere():
                        adj[1] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[1] = MarkerEnd(next_marker, head=next_marker.reverse)
                    adjacencies[i].append(Adjacency(left_end_marker=adj[0], right_end_marker=adj[1], label=label))

                    # Create reference matrix for A and B
                    if adj[0].marker != Telomere():
                        references[i][adj[0]] = len(adjacencies[i]) - 1
                    if adj[1].marker != Telomere():
                        references[i][adj[1]] = len(adjacencies[i]) - 1

                    current_marker = next_marker
                    index += 1
                    adj_index += 1

        # print("self.adjA: ", self.adjA)
        # print("self.adjB: ", self.adjB)
        # Traverse the graph
        to_visit_a_index = set(range(len(self.adjA)))
        to_visit_b_index = set(range(len(self.adjB)))
        visited_a_index = set()
        visited_b_index = set()

        cycles = 0
        ab_paths = 0
        a_runs = 0
        on_a_run = False
        b_runs = 0
        on_b_run = False

        while to_visit_a_index:
            current_adj_index = to_visit_a_index.pop()
            visited_a_index.add(current_adj_index)

            left_marker = self.adjA[current_adj_index].left_end_marker
            right_marker = self.adjA[current_adj_index].right_end_marker

            # Check label for A_run
            if self.adjA[current_adj_index].label:
                a_runs += 1
                on_a_run = True

            # Checks for AB paths
            paths_end_on_a = [True, True]
            for i, current_marker in enumerate([left_marker, right_marker]):
                a_side = True
                next_adj_index = reference_B[current_marker]
                if next_adj_index in visited_b_index:
                    # If cycle on left_marker, no need to check right marker
                    continue
                while next_adj_index is not None:
                    current_adj_index = next_adj_index
                    a_side = not a_side
                    adj_side = self.adjA if a_side else self.adjB
                    current_adj = adj_side[current_adj_index]
                    reference_side = reference_B if a_side else reference_A

                    visited_a_index.add(current_adj_index) if a_side else visited_b_index.add(current_adj_index)
                    to_visit_a_index.remove(current_adj_index) if a_side else to_visit_b_index.remove(current_adj_index)

                    # Check label, and update runs
                    if current_adj.label:
                        if not a_side and on_a_run:
                            on_b_run = True
                            on_a_run = False
                            b_runs += 1
                        elif not a_side and not on_b_run:
                            on_b_run = True
                            b_runs += 1

                        elif a_side and on_b_run:
                            on_a_run = True
                            on_b_run = False
                            a_runs += 1
                        elif a_side and not on_a_run:
                            on_a_run = True
                            a_runs += 1


                    # Switch markers

                    if current_marker != current_adj.left_end_marker:
                        current_marker = current_adj.left_end_marker
                    else:
                        current_marker = current_adj.right_end_marker

                    # Find next location
                    next_adj_index = reference_side[current_marker]

                    if next_adj_index is None:
                        # Path ends
                        paths_end_on_a[i] = a_side
                        break
                    elif (a_side and next_adj_index in visited_b_index)\
                            or (not a_side and next_adj_index in visited_a_index):
                        # cycle
                        cycles += 1
                        break
                    else:
                        # The path continues
                        pass

            # AB_path
            if paths_end_on_a[0] != paths_end_on_a[1]:
                ab_paths += 1

        self.cycles = cycles
        self.ab_paths = ab_paths
        self.a_runs = a_runs
        self.b_runs = b_runs
        self.run_potential = a_runs + b_runs
        self.indel_potential = (self.run_potential + 1) // 2 + ((self.run_potential // 2) % 2) if self.run_potential > 0 else 0


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


def calculate_distance(A, B, method="dcj"):
    """Calculate the DCJ distance between genome A and B"""
    try:
        ["dcj", "indel"].index(method)
    except ValueError:
        print('method must be "dcj"(default) or "indel')

    # DCJ distance
    # DCJ(A, B) = |G| - c - b/2
    # Common markers in formula excludes telomeres
    # Clear Dictionary for creating dictionaries
    Chromosome.clear()
    ag = AdjacencyGraph(A, B)
    dcj_distance = len(ag.commonMarkers) - ag.cycles - ag.ab_paths / 2
    dcj_indel_distance = dcj_distance + ag.indel_potential

    if method == "indel":
        return dcj_indel_distance
    return dcj_distance


def validate_genome(genome):
    """
    Validates the form of the genome. Acceptable forms are


    [".abcdef.", ".skjkjfd."]
    :param genome:
    :return:
    """


    #Restraints, no duplicated markers within one genome.
    #Check Genome is a listsa

    print(Genome(*[Chromosome(c) for c in genome]))



    # If starts with telomere, must end with telomere
    #[".abcde.", "rtskli"]

if __name__ == '__main__':
    # b = Genome(Chromosome("ab"), Chromosome(".cd."), Chromosome(".e."), Chromosome("fg"))
    # a = Genome(Chromosome(".acD."), Chromosome("be"), Chromosome(".fg."))
    a = ["aa"]
    # res = AdjacencyGraph(a, b)
    # print(res.ab_paths)
    # print(res.commonMarkers)
    # print(res.cycles)
    # print(calculate_distance(a, b, method="dcj"))
    # A = Chromosome(".ab.")
    # markers = set([x for x in A])
    # print(markers)
    validate_genome(a)


