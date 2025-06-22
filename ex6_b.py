from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

# Αφαίρεση των 'X' και εύρεση πρώτου/τελευταίου γνωστού κομματιού
raw_mystery_sequence = """
GATGAPGIAGAPGFPGARGAPGPQGPSGAPGPKXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGVQGPPGPQGPR
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGSAGPPGATGFP
GAAGRXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXGVVGLPGQR
""".replace('\n', '')

matches = list(re.finditer(r'[^X]+', raw_mystery_sequence))
first_known = matches[0].group()
last_known = matches[-1].group()
total_length = len(raw_mystery_sequence.replace('\n', ''))

query_sequence = first_known + last_known
query_record = SeqRecord(Seq(query_sequence), id="MysteryEnds", description="First and last known parts")

def run_blastp(seq_record, db="nr", hits=20):
    print("Running BLASTP...")
    result_handle = NCBIWWW.qblast("blastp", db, seq_record.format("fasta"), hitlist_size=hits, format_type="XML")
    with open("blast_results_partial.xml", "w") as f:
        f.write(result_handle.read())
    result_handle.close()
    print("Αποθηκεύτηκαν τα αποτελέσματα στο blast_results_partial.xml.")

run_blastp(query_record)

# Διάβασμα και φιλτράρισμα αποτελεσμάτων, αποθήκευση σε αρχείο
with open("blast_results_partial.xml") as blast_file, open("filtered_hits.txt", "w", encoding="utf-8") as out_file:
    blast_records = NCBIXML.parse(blast_file)
    found = False
    for record in blast_records:
        for alignment in record.alignments:
            org_name = None
            if "[" in alignment.title and "]" in alignment.title:
                org_name = alignment.title.split("[")[-1].split("]")[0]
            # Έλεγχος αν κάποια λέξη του ονόματος ξεκινά με "Tyr"
            if org_name and any(word.startswith("Tyr") for word in org_name.split()) and abs(alignment.length - total_length) < 30:
                found = True
                result_str = (
                    f"\nΒρέθηκε hit με οργανισμό που περιέχει λέξη που ξεκινά με 'Tyr' και κατάλληλο μήκος:\n"
                    f"ID: {alignment.accession}\n"
                    f"Organism: {org_name}\n"
                    f"Description: {alignment.title}\n"
                    f"Length: {alignment.length}\n"
                    f"Query length: {total_length}\n"
                )
                print(result_str)
                out_file.write(result_str)
                for hsp in alignment.hsps:
                    hsp_str = (
                        f"E-value: {hsp.expect}\n"
                        f"Identity: {hsp.identities}/{hsp.align_length} ({100 * hsp.identities / hsp.align_length:.2f}%)\n"
                        f"Query: {hsp.query}\n"
                        f"Match: {hsp.match}\n"
                        f"Subject: {hsp.sbjct}\n"
                        + "-" * 60 + "\n"
                    )
                    print(hsp_str)
                    out_file.write(hsp_str)
    if not found:
        msg = "Δεν βρέθηκε οργανισμός που να περιέχει λέξη που να ξεκινά με 'Tyr' και κατάλληλο μήκος.\n"
        print(msg)
