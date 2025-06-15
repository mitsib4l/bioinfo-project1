from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Χρησιμοποιήστε ΜΟΝΟ τα γνωστά αμινοξέα (χωρίς X)
mystery_sequence = "GATGAPGIAGAPGFPGARGAPGPQGPSGAPGPKGVQGPPGPQGPRGSAGPPGATGFPGAAGRGVVGLPGQR"

print("Εκτέλεση BLASTP...")
blast_results_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=mystery_sequence,
    expect=1e-10,
    hitlist_size=5
)

# Αποθήκευση και ανάλυση αποτελεσμάτων
with open("blast_results.xml", "w") as out_handle:
    out_handle.write(blast_results_handle.read())

with open("blast_results.xml") as result_handle:
    print("\nΑποτελέσματα BLASTP:")
    blast_records = NCBIXML.parse(result_handle)
    found = False
    for record in blast_records:
        for alignment in record.alignments:
            found = True
            print("\n--- Ομόλογη Πρωτεΐνη ---")
            print(f"Τίτλος: {alignment.title}")
            print(f"Οργανισμός: {alignment.hit_def}")
            print(f"Μήκος: {alignment.length} αμινοξέα")
            for hsp in alignment.hsps:
                print(f"\nE-value: {hsp.expect}")
                print(f"Identity: {hsp.identities}/{hsp.align_length} ({100 * hsp.identities / hsp.align_length:.2f}%)")
                print(f"Query: {hsp.query}")
                print(f"Match: {hsp.match}")
                print(f"Subject: {hsp.sbjct}")
    if not found:
        print("Δεν βρέθηκαν ομόλογες πρωτεΐνες. Δοκιμάστε να αφαιρέσετε τα 'X' ή να αλλάξετε το φίλτρο οργανισμού.")

print("\nΗ ανάλυση ολοκληρώθηκε!")