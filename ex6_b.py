from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

# Ορισμός της άγνωστης αλληλουχίας (αγνοώντας τα 'X')
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

mystery_seq = raw_mystery_sequence.replace('X', '')

# Υποβολή BLASTP ερωτήματος
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=mystery_seq,
    expect=0.1,  # Αυστηρό κατώφλι
    composition_based_statistics=True,
    hitlist_size=50,  # Περιορισμός στον αριθμό των hits
    format_type="XML"  # Επιστροφή αποτελεσμάτων σε XML μορφή
)

# Αποθήκευση των αποτελεσμάτων
with open("blast_results.xml", "w") as f:
    f.write(result_handle.read())
result_handle.close()

query_length = 234  

with open("blast_results.xml") as blast_file:
    blast_records = NCBIXML.parse(blast_file)
    for record in blast_records:
        for alignment in record.alignments:
            hit_length = alignment.length
            if hit_length == query_length:  # Φιλτράρισμα μήκους
                for hsp in alignment.hsps:
                    if hsp.expect < 0.1:  # Επιπλέον φιλτράρισμα E-value
                        print(f"Ομόλογη πρωτεΐνη: {alignment.title}")
                        print(f"Οργανισμός: {alignment.hit_def}")
                        print(f"Μήκος: {hit_length}")
                        print(f"E-value: {hsp.expect}")
                        print(f"Αναλογία: {hsp.identities / hsp.align_length * 100:.2f}%")
                        print("------")






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