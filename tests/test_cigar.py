import parasail
from unittest import TestCase, main


def print_cigar_attributes(cigar):
    print(cigar)
    print(cigar.seq)
    print(cigar.len)
    print(cigar.beg_query)
    print(cigar.beg_ref)
    print(cigar.decode)


def print_traceback_attributes(traceback):
    print(traceback)
    print(traceback.query)
    print(traceback.comp)
    print(traceback.ref)


class Tests(TestCase):

    def test0(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.blosum62)
        with self.assertRaises(AttributeError):
            print_cigar_attributes(result.cigar)

    def test1(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.blosum62)
        with self.assertRaises(AttributeError):
            print_traceback_attributes(result.traceback)

    def test2(self):
        result = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        print_cigar_attributes(result.cigar)

    def test3(self):
        result = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        print_traceback_attributes(result.traceback)

    def test4(self):
        """ test simple error free alignment """
        ref = "ATCGATGTGC"
        seq = "ATCGATGTGC"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(ref) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "10=")

    def test5(self):
        """ test simple two base deletion """
        ref = "ATCGATGTGC"
        seq = "ATCGATGT"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(seq) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "%s=" % len(seq))

    def test6(self):
        """ test simple 2 base insertion """
        ref = "ATCGATGTGC"
        seq = "ATCGTTATGTGC"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(ref) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "4=2I6=")

    def test7(self):
        """ test offset start """
        ref = "ATCGATGTGC"
        seq =    "GATGTGC"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(ref) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "3D7=")

    def test8(self):
        """ test offset end """
        ref = "ATCGATGTGC"
        seq = "ATCGATG"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(seq) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "7=")

    def test9(self):
        """ test additional start """
        ref = "ATCGATGTGC"
        seq = "TTTTATCGATGTGC"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(ref) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "4S10=")

    def test10(self):
        """ test additional end """
        ref = "ATCGATGTGC"
        seq = "ATCGATGTGCTTTT"
        result = parasail.sw_trace_striped_16(seq, ref, 8, 4, parasail.dnafull)
        self.assertEqual((result.cigar.beg_ref, result.end_ref), (0, len(ref) - 1))
        self.assertEqual((result.cigar.beg_query, result.end_query), (0, len(seq) - 1))
        self.assertEqual(result.cigar.decode.decode(), "10=4S")


if __name__ == '__main__':
    main()
