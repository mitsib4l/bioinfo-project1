def global_alignment(S1, S2, ρ, σ):
    n, m = len(S1), len(S2)
    DP = [[0] * (m + 1) for _ in range(n + 1)]
    C = [[0] * (m + 1) for _ in range(n + 1)]  # Τρέχον μήκος ασυμφωνιών

    # Αρχικοποίηση
    for i in range(1, n + 1):
        DP[i][0] = DP[i-1][0] - ρ
        C[i][0] = 0
    for j in range(1, m + 1):
        DP[0][j] = DP[0][j-1] - ρ
        C[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Match ή mismatch
            if S1[i-1] == S2[j-1]:
                match_score = DP[i-1][j-1] + 1
                mismatch_score = -float('inf')
                k = 0
            else:
                k = C[i-1][j-1] + 1
                mismatch_score = DP[i-1][j-1] - (ρ + σ * k)
                match_score = -float('inf')

            # Εισαγωγή ή διαγραφή
            insert_score = DP[i-1][j] - ρ
            delete_score = DP[i][j-1] - ρ

            # Ενημέρωση DP[i][j] και C[i][j]
            max_score = max(match_score, mismatch_score, insert_score, delete_score)
            DP[i][j] = max_score

            if max_score == match_score:
                C[i][j] = 0
            elif max_score == mismatch_score:
                C[i][j] = k
            else:
                C[i][j] = 0  # Εισαγωγή/διαγραφή μηδενίζει τις ασυμφωνίες

    return DP[n][m]