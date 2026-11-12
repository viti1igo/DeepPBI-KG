"""Run from anywhere: python tests/test_go_extract.py

Small inline fixtures only, so this runs without the repo's data. The DIAMOND test runs only
when DIAMOND is installed.
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import go_extract
from go_extract import best_hits, read_emapper, read_interproscan, sprot_fasta

SEQ_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACRKCDRK"
SEQ_B = "MSLNEQQLAQWLTDSGLTPAQAAELIRLHSEGWSPRQIAEHFGLSRSTVHRALKAL"
FLAT = f"""ID   A_ECOLI Reviewed; {len(SEQ_A)} AA.
AC   P00001; Q00009;
OS   Escherichia coli.
OC   Bacteria.
OX   NCBI_TaxID=562;
DR   GO; GO:0005737; C:cytoplasm; IEA:X.
DR   GO; GO:0003677; F:DNA binding; IDA:X.
SQ   SEQUENCE   {len(SEQ_A)} AA;
     {SEQ_A[:60]}
     {SEQ_A[60:]}
//
ID   B_ECOLI Reviewed; {len(SEQ_B)} AA.
AC   P00002;
OS   Escherichia coli.
OC   Bacteria.
OX   NCBI_TaxID=562;
DR   GO; GO:0016020; C:membrane; IEA:X.
DR   GO; GO:9999999; C:retired; IEA:X.
SQ   SEQUENCE   {len(SEQ_B)} AA;
     {SEQ_B}
//
"""


def hit(q, s, pident, cov, evalue, bits):
    return {"qseqid": q, "sseqid": s, "pident": pident, "qcovhsp": cov, "evalue": evalue, "bitscore": bits}


class ExtractTest(unittest.TestCase):
    def test_sprot_fasta_reads_first_accession_and_whole_sequence(self):
        self.assertEqual(list(sprot_fasta(FLAT.splitlines(True))), [("P00001", SEQ_A), ("P00002", SEQ_B)])

    def test_best_hit_passes_every_threshold_then_wins_on_bit_score(self):
        rows = [hit("q1", "low_id", 20, 90, 1e-30, 500), hit("q1", "ok", 45, 80, 1e-20, 200),
                hit("q1", "better", 60, 70, 1e-25, 300), hit("q2", "short", 90, 30, 1e-40, 400),
                hit("q3", "weak", 90, 90, 1e-3, 100)]
        best = best_hits(rows, min_identity=30, min_coverage=50, max_evalue=1e-5)
        self.assertEqual({q: h["sseqid"] for q, h in best.items()}, {"q1": "better"})

    def test_readers(self):
        with tempfile.TemporaryDirectory() as folder:
            em, ips = os.path.join(folder, "e.tsv"), os.path.join(folder, "i.tsv")
            with open(em, "w") as handle:
                handle.write("## eggNOG-mapper\n#query\tseed_ortholog\tGOs\n"
                             "g1\tx\tGO:0005737,GO:0003677\ng2\tx\t-\n## done\n")
            with open(ips, "w") as handle:
                cols = ["p1", "md5", "100", "Pfam", "PF1", "d", "1", "50", "1e-9", "T", "date", "IPR1", "d"]
                handle.write("\t".join(cols + ["GO:0005524(InterPro)|GO:0016887(PANTHER)"]) + "\n")
                handle.write("\t".join(cols + ["GO:0005524"]) + "\n")
                handle.write("\t".join(["p2"] + cols[1:]) + "\n")
            self.assertEqual(read_emapper(em), {"g1": {"GO:0005737", "GO:0003677"}, "g2": set()})
            self.assertEqual(read_interproscan(ips), {"p1": {"GO:0005524", "GO:0016887"}, "p2": set()})

    @unittest.skipUnless(os.path.exists(os.path.join(os.path.dirname(sys.executable), "diamond")),
                         "DIAMOND not installed")
    def test_fasta_to_pairs_end_to_end(self):
        obo = os.path.join(ROOT, "tests", "fixtures", "mini.obo")
        with tempfile.TemporaryDirectory() as folder:
            sprot = os.path.join(folder, "sprot.dat.gz")
            with gzip.open(sprot, "wt") as handle:
                handle.write(FLAT)
            faa = os.path.join(folder, "q.faa")
            with open(faa, "w") as handle:
                handle.write(f">mine_a\n{SEQ_A[5:]}\n>mine_b\n{SEQ_B}\n>nothing\n{'W' * 80}\n")
            script = os.path.join(ROOT, "scripts", "go_extract.py")
            subprocess.run([sys.executable, script, "build-db", sprot, os.path.join(folder, "db")],
                           check=True, capture_output=True)
            subprocess.run([sys.executable, script, "run", faa, os.path.join(folder, "out"), "--db",
                            os.path.join(folder, "db"), "--sprot", sprot, "--obo", obo],
                           check=True, capture_output=True)
            with open(os.path.join(folder, "out", "pairs.tsv")) as handle:
                pairs = sorted(line.split() for line in handle)
            with open(os.path.join(folder, "out", "provenance.json")) as handle:
                prov = json.load(handle)
        self.assertEqual(pairs, [["mine_a", "GO:0003677"], ["mine_a", "GO:0005737"], ["mine_b", "GO:0016020"], ["mine_b", "GO:9999999"]])
        self.assertEqual((prov["proteins"], prov["proteins_with_hit"]), (3, 2))
        self.assertEqual(prov["database"]["entries"], 2)
        self.assertEqual((prov["pairs_dropped_retired_or_unknown"], prov["genes_left_without_go"]), (1, 0))


if __name__ == "__main__":
    unittest.main()
