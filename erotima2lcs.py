def read_fasta(filename):
    """
    Διαβάζει μια ακολουθία από αρχείο FASTA (αγνοεί το header).
    Επιστρέφει την ακολουθία ως string.
    """
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue  # Αγνόησε το header
            seq.append(line.strip())
    return ''.join(seq)

class LongestCommonSubstring:
    """
    Κλάση για εύρεση της μέγιστης κοινής υποσυμβολοσειράς (Longest Common Substring)
    μεταξύ δύο ακολουθιών S1 και S2.
    """

    def __init__(self, S1, S2):
        self.S1 = S1
        self.S2 = S2

    def longest_common_substring(self):
        """
        Υπολογίζει το μήκος και το ίδιο το string της μέγιστης κοινής υποσυμβολοσειράς.
        Επιστρέφει: (μήκος, substring)
        """
        n, m = len(self.S1), len(self.S2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        max_len = 0
        end_pos = 0  # Τέλος της υποσυμβολοσειράς στο S1

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if self.S1[i - 1] == self.S2[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1] + 1
                    if dp[i][j] > max_len:
                        max_len = dp[i][j]
                        end_pos = i
                else:
                    dp[i][j] = 0  # Δεν υπάρχει κοινή υποσυμβολοσειρά εδώ

        return max_len, self.S1[end_pos - max_len:end_pos]

if __name__ == "__main__":
    # Διάβασε τις ακολουθίες από τα αρχεία
    seq1 = read_fasta("A0A023SFE5.txt")
    seq2 = read_fasta("NC_0455122.txt")

    # Εφάρμοσε Longest Common Substring (υποσυμβολοσειρά)
    lcs_substr = LongestCommonSubstring(seq1, seq2)
    substr_len, substr = lcs_substr.longest_common_substring()
    print(f"Longest Common Substring length: {substr_len}")
    print(f"Longest Common Substring: {substr}")