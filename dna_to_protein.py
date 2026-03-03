# Create a sequence object for the coding DNA strand
from Bio.Seq import Seq

cDNA = Seq("ATGGCCCTATAGTGTCTAAGCTAG")
print("DNA Sequence:", cDNA)
print("Sequence Length:", len(cDNA))

# Transcription of coding DNA strand to mRNA

mRNA = cDNA.transcribe()
print("mRNA Sequence:", mRNA)

# Translation of mRNA to a protein sequence

protein_full = mRNA.translate(to_stop=True)
print("Protein Sequence:", protein_full)


#Stop at First Stop Codon

protein = mRNA.translate(to_stop=True)
print("Protein(to stop):", protein)


#Translate directly from DNA

protein_full = cDNA.translate()
print("Protein Sequence:", protein_full)


#Reverse Complement (Coding → Template Strand)

template_strand = cDNA.reverse_complement()
print("Template Strand:", template_strand)


print("-" * 50)
print("CDS Validation Check")


# -----------------------------
# CDS Validation Checks
# -----------------------------

valid = True

# 1️⃣ Length divisible by 3
if len(cDNA) % 3 != 0:
    print("Invalid CDS: Length not divisible by 3")
    valid = False

# 2️⃣ Start codon check
if not cDNA.startswith("ATG"):
    print("Invalid CDS: Missing start codon (ATG)")
    valid = False

# 3️⃣ Stop codon check
if cDNA[-3:] not in ["TAA", "TAG", "TGA"]:
    print("Invalid CDS: Missing stop codon")
    valid = False

# 4️⃣ Internal stop codon check
if "*" in protein_full[:-1]:
    print("Invalid CDS: Internal stop codon detected")
    valid = False


# Final Result
if valid:
    print("Sequence is a valid Coding DNA Sequence (CDS)")
else:
    print("Sequence failed CDS validation")
