import re

def read_fasta(filename):
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue  # Αγνόησε το header
            # Κράτα μόνο γράμματα (A-Z, a-z)
            clean_line = re.sub(r'[^A-Za-z]', '', line)
            seq.append(clean_line.upper())
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

class SequenceAlignment:
    def __init__(self, S1, S2):
        self.S1 = S1
        self.S2 = S2
        self.n = len(S1)
        self.m = len(S2)

    def lcs(self):
        n, m = self.n, self.m
        S1, S2 = self.S1, self.S2
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if S1[i - 1] == S2[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1] + 1
                else:
                    dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
        lcs_str = []
        i, j = n, m
        while i > 0 and j > 0:
            if S1[i - 1] == S2[j - 1]:
                lcs_str.append(S1[i - 1])
                i -= 1
                j -= 1
            elif dp[i - 1][j] > dp[i][j - 1]:
                i -= 1
            else:
                j -= 1
        lcs_str.reverse()
        return dp[n][m], ''.join(lcs_str)

    def edit_distance_weighted(self, d=1, r=2, e=0):
        n, m = self.n, self.m
        S1, S2 = self.S1, self.S2
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            dp[i][0] = i * d
        for j in range(m + 1):
            dp[0][j] = j * d
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if S1[i - 1] == S2[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1] + e
                else:
                    dp[i][j] = min(
                        dp[i - 1][j] + d,
                        dp[i][j - 1] + d,
                        dp[i - 1][j - 1] + r
                    )
        return dp[n][m]

if __name__ == "__main__":
    # Διάβασε τις ακολουθίες από τα αρχεία με σωστό preprocessing
    seq1 = read_fasta("A0A023SFE5.txt")
    seq2 = read_fasta("NC_0455122.txt")

    print(f"Μήκος ακολουθίας 1: {len(seq1)}")
    print(f"Μήκος ακολουθίας 2: {len(seq2)}")

    # Υπολογισμός και εμφάνιση Longest Common Substring
    lcs_substr = LongestCommonSubstring(seq1, seq2)
    substr_len, substr = lcs_substr.longest_common_substring()
    print(f"Longest Common Substring length: {substr_len}")
    print(f"Longest Common Substring: {substr}")

    # Υπολογισμός και εμφάνιση LCS και edit distance
    aligner = SequenceAlignment(seq1, seq2)
    lcs_len, lcs_str = aligner.lcs()
    print(f"LCS length: {lcs_len}")

    edit_dist = aligner.edit_distance_weighted()
    print(f"Weighted Edit Distance (d=1, r=2, e=0): {edit_dist}")