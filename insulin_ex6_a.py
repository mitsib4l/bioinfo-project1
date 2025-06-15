from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import subprocess
import os
from Bio.SeqRecord import SeqRecord

Entrez.email = "dimitaaaa09@gmail.com"  

# Λήψη της πρωτεΐνης INS (P01308) από UniProt
def fetch_protein(uniprot_id):
    handle = Entrez.efetch(db="protein", id=uniprot_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

human_insulin = fetch_protein("P01308")
print(f"Human Insulin: {human_insulin.id}, {len(human_insulin.seq)} amino acids")

# BLASTP αναζήτηση για ομόλογες ακολουθίες
def run_blastp(seq_record, db="nr", hits=8):
    print("Running BLASTP...")
    result = NCBIWWW.qblast("blastp", db, seq_record.format("fasta"), hitlist_size=hits)
    blast_records = NCBIXML.parse(result)
    blast_record = next(blast_records)  # Πάρτε το πρώτο αποτέλεσμα
    hits_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            hits_list.append({
                "id": alignment.accession,
                "desc": alignment.title,
                "evalue": hsp.expect,
                "seq": hsp.sbjct
            })
            break  # Πάρτε μόνο το top HSP για κάθε alignment
    return hits_list[:hits]

blast_results = run_blastp(human_insulin)
print(f"Found {len(blast_results)} homologs.")

# Αποθήκευση των BLAST hits σε ένα αρχείο FASTA
with open("insulin_homologs.fasta", "w") as f:
    SeqIO.write(human_insulin, f, "fasta")
    for hit in blast_results:
        hit_record = SeqRecord(seq=hit["seq"], id=hit["id"], description=hit["desc"])
        SeqIO.write(hit_record, f, "fasta")

# Pairwise alignments (Human Insulin vs κάθε ομόλογη ακολουθία)
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

with open("pairwise_alignments.aln", "w") as f:
    for hit in blast_results:
        target_seq = hit["seq"]
        alignment = aligner.align(human_insulin.seq, target_seq)[0]  # Πάρτε το top alignment
        f.write(f">Human_INS vs {hit['id']}\n")
        f.write(str(alignment) + "\n\n")

# Δημιουργία αντικειμένου MultipleSeqAlignment με την ανθρώπινη ινσουλίνη και όλους τους ομόλογους στοιχισμένους σε ζεύγη
msa_records = [SeqRecord(human_insulin.seq, id=human_insulin.id, description="Human Insulin")]
for hit in blast_results:
    alignment = aligner.align(human_insulin.seq, hit["seq"])[0]
    # Εξαγωγή στοιχισμένων ακολουθιών από το αντικείμενο στοιχοθέτησης
    aligned_human = alignment.aligned[0]
    aligned_hit = alignment.aligned[1]
    # Για απλότητα, χρησιμοποιήστε τις στοιχισμένες ακολουθίες ως συμβολοσειρές
    msa_records.append(SeqRecord(alignment.target, id=hit["id"], description=hit["desc"]))

msa = MultipleSeqAlignment(msa_records)

# Save all sequences (human + homologs) to FASTA for Clustal Omega
with open("insulin_homologs.fasta", "w") as f:
    SeqIO.write(human_insulin, f, "fasta")
    for hit in blast_results:
        hit_record = SeqRecord(seq=hit["seq"], id=hit["id"], description=hit["desc"])
        SeqIO.write(hit_record, f, "fasta")

# Run Clustal Omega for multiple sequence alignment
subprocess.run([
    "C:\Program Files\clustal-omega-1.2.2-win64\clustalo.exe", "-i", "insulin_homologs.fasta", "-o", "multiple_alignment.aln", "--outfmt=clu", "--force"
], check=True)

print("Ολοκληρώθηκαν οι στοιχίσεις! Αποθηκεύτηκαν τα αρχεία:")
print("- pairwise_alignments.aln (στοιχίσεις ζευγών)")
print("- multiple_alignment.aln (πολλαπλή στοίχιση)")