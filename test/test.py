import unittest
from dcj.main import Telomere, Marker, Genome, Chromosome, calculate_distance, transform_genome


class TestMarkers(unittest.TestCase):
    def test_telomere_type(self):
        a = Telomere()
        self.assertTrue(type(a) is Telomere)

    def test_forward_reverse_marker(self):
        self.assertEquals(Marker(uid=0, dna="d", reverse=True), Marker(uid=0, dna="d", reverse=False))

    def test_telomere_repr(self):
        self.assertEquals(str(Telomere()), ".")

    def test_reverse_equality(self):
        self.assertNotEqual(Marker(uid=0, dna="D", reverse=True), Marker(uid=1, dna="D", reverse=True))


class TestAdjacencyGraphLinear(unittest.TestCase):
    def test_dcj_single_a_dcj(self):
        a = Genome(Chromosome(".ac."))
        b = Genome(Chromosome(".abc."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 0)

    def test_dcj_single_a_indel(self):
        a = Genome(Chromosome(".ac."))
        b = Genome(Chromosome(".abc."))
        self.assertEquals(calculate_distance(a, b, method="indel"), 1)

    def test_dcj_single_swap_dcj(self):
        a = Genome(Chromosome(".acb."))
        b = Genome(Chromosome(".abc."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 2)

    def test_dcj_intermediate_dcj(self):
        a = Genome(Chromosome(".cbAdEGf."))
        b = Genome(Chromosome(".abcdefg."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 6)

    def test_dcj_intermediate_2_dcj(self):
        a = Genome(Chromosome(".aBdCE."))
        b = Genome(Chromosome(".abcde."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 4)


class TestAdjacencyGraphCircular(unittest.TestCase):
    def test_dcj_backwards(self):
        a = Genome(Chromosome("abc"))
        b = Genome(Chromosome("cba"))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 2)

    def test_dcj_single_split(self):
        a = Genome(Chromosome(".ab."))
        b = Genome(Chromosome("ab"))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 1)

    def test_dcj_circular_excision(self):
        a = Genome(Chromosome("ab"), Chromosome(".cd."))
        b = Genome(Chromosome(".abcd."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 1)

    def test_dcj_default_example(self):
        a = Genome(Chromosome("ab"), Chromosome(".cd."), Chromosome(".e."), Chromosome("fg"))
        b = Genome(Chromosome(".acD."), Chromosome("be"), Chromosome(".fg."))
        self.assertEquals(calculate_distance(a, b, method="dcj"), 5)


class TestValidateGenome(unittest.TestCase):
    def test_duplicated_markers(self):
        self.assertRaises(ValueError, transform_genome, ["aa"])

    def test_duplicated_marker_reversed(self):
        self.assertRaises(ValueError, transform_genome, ["Aa"])

    def test_duplicated_marker_between_chromosomes(self):
        self.assertRaises(ValueError, transform_genome, ["abcde", ".fghiC."])

    def test_missing_telomere(self):
        self.assertRaises(ValueError, transform_genome, [".abcd"])

    def test_missing_telomere(self):
        self.assertRaises(ValueError, transform_genome, ["abcd."])

    def test_rogue_telomere(self):
        self.assertRaises(ValueError, transform_genome, [".abc.de."])

    def test_genome_without_chromosomes(self):
        self.assertRaises(TypeError, transform_genome, "abcde", ".fg.")

    def test_chromosome_list_not_markers(self):
        self.assertRaises(TypeError, transform_genome, [["a", "b", "c"]])


if __name__ == '__main__':
    unittest.main()
