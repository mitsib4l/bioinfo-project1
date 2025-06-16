from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_blastp(seq_record, db="nr", hits=10):
    print("Running BLASTP...")
    result_handle = NCBIWWW.qblast("blastp", db, seq_record.format("fasta"), hitlist_size=hits, format_type="XML")
    with open("blast_results.xml", "w") as f:
        f.write(result_handle.read())
    result_handle.close()
    print("Αποθηκεύτηκαν τα αποτελέσματα στο blast_results.xml.")

# Αφαίρεση των 'X' από την ακολουθία
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

mystery_sequence = raw_mystery_sequence.replace('X', '')
mystery_record = SeqRecord(Seq(mystery_sequence), id="Mystery", description="Mystery protein (no X)")
query_length = len(mystery_sequence)

# Εκτέλεση BLAST και αποθήκευση αποτελεσμάτων
run_blastp(mystery_record)

# Διάβασμα του αρχείου αποτελεσμάτων και εκτύπωση μόνο για οργανισμό που ξεκινά με "Tyr"
with open("blast_results.xml") as blast_file:
    blast_records = NCBIXML.parse(blast_file)
    found = False
    for record in blast_records:
        for alignment in record.alignments:
            org_name = None
            if "[" in alignment.title and "]" in alignment.title:
                org_name = alignment.title.split("[")[-1].split("]")[0]
            if org_name and org_name.startswith("Tyr"):
                found = True
                print("\nΒρέθηκε hit με οργανισμό που ξεκινά με 'Tyr':")
                print(f"ID: {alignment.accession}")
                print(f"Organism: {org_name}")
                print(f"Description: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"Query length: {query_length}")
                for hsp in alignment.hsps:
                    print(f"E-value: {hsp.expect}")
                    print(f"Identity: {hsp.identities}/{hsp.align_length} ({100 * hsp.identities / hsp.align_length:.2f}%)")
                    print(f"Query: {hsp.query}")
                    print(f"Match: {hsp.match}")
                    print(f"Subject: {hsp.sbjct}")
                print("-" * 60)
    if not found:
        print("Δεν βρέθηκε οργανισμός που να ξεκινά με 'Tyr'.")






"""
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_blastp(seq_record, db="nr", hits=10):
    print("Running BLASTP...")
    result = NCBIWWW.qblast("blastp", db, seq_record.format("fasta"), hitlist_size=hits)
    blast_records = NCBIXML.parse(result)
    blast_record = next(blast_records)
    hits_list = []
    for alignment in blast_record.alignments:
        org_name = None
        if "[" in alignment.title and "]" in alignment.title:
            org_name = alignment.title.split("[")[-1].split("]")[0]
        for hsp in alignment.hsps:
            hits_list.append({
                "id": alignment.accession,
                "desc": alignment.title,
                "evalue": hsp.expect,
                "seq": hsp.sbjct,
                "organism": org_name
            })
            break
    return hits_list

# Αφαίρεση των 'X' από την ακολουθία
raw_mystery_sequence = """
#GATGAPGIAGAPGFPGARGAPGPQGPSGAPGPKXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGVQGPPGPQGPR
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGSAGPPGATGFP
#GAAGRXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXGVVGLPGQR
""".replace('\n', '')

mystery_sequence = raw_mystery_sequence.replace('X', '')
mystery_record = SeqRecord(Seq(mystery_sequence), id="Mystery", description="Mystery protein (no X)")

blast_results = run_blastp(mystery_record)

print("Όλα τα hits:")
for hit in blast_results:
    print(f"ID: {hit['id']}, Organism: {hit['organism']}, Desc: {hit['desc']}, E-value: {hit['evalue']}")

tyr_hits = [hit for hit in blast_results if hit["organism"] and hit["organism"].startswith("Tyr")]

if tyr_hits:
    print("\nHits με οργανισμό που ξεκινά με 'Tyr':")
    for hit in tyr_hits:
        print(f"ID: {hit['id']}, Organism: {hit['organism']}, Desc: {hit['desc']}, E-value: {hit['evalue']}")
else:
    print("\nΔεν βρέθηκε οργανισμός που να ξεκινά με 'Tyr'.")
"""