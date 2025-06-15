def build_generalized_suffix_tree(strings):
    # Υλοποίηση κατασκευής GST (π.χ., με Ukkonen's algorithm)
    pass

def compute_lcp_pairs(gst):
    # Επιστρέφει το a (συνολικό μήκος όλων των LCP)
    pass

def max_common_prefix(strings):
    gst = build_generalized_suffix_tree(strings)
    alpha = compute_lcp_pairs(gst)
    return alpha
