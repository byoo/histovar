import numpy as np

class EMParalogAssigner:
    def __init__(self, loci, max_iter=100, tol=1e-4):
        self.loci = loci
        self.max_iter = max_iter
        self.tol = tol

    def fit_predict(self, reads_info):
        """
        reads_info: list of dicts, each dict:
            {"read_id": "r1", "H3C14": 0.9, "H3C15": 0.1, ...}
        """
        if not reads_info:
            return {}, {}

        # --- Build matrix of mapping probabilities ---
        read_ids = [r["read_id"] for r in reads_info]
        n_reads = len(reads_info)
        n_loci = len(self.loci)

        locus_to_idx = {locus: idx for idx, locus in enumerate(self.loci)}
        # X[j, k] = mapping probability of read j to locus k
        X = np.zeros((n_reads, n_loci), dtype=float)
        for j, r in enumerate(reads_info):
            for locus, val in r.items():
                if locus in locus_to_idx:
                    X[j, locus_to_idx[locus]] = val

        # --- Initialize theta uniformly ---
        theta = np.full(n_loci, 1.0 / n_loci, dtype=float)

        # --- EM iterations ---
        for _ in range(self.max_iter):
            # E-step: responsibilities z[j, k]
            # denom[j] = sum_k X[j,k] * theta[k]
            denom = X @ theta  # shape (n_reads,)
            # Avoid division by zero: where denom == 0, keep zeros
            z = np.zeros_like(X)
            mask = denom > 0
            z[mask] = X[mask] * theta  # broadcasting theta over loci
            z[mask] /= denom[mask][:, None]

            # M-step: update theta
            new_theta = z.sum(axis=0) / n_reads

            # Convergence check
            diff = np.abs(new_theta - theta).sum()
            theta = new_theta
            if diff < self.tol:
                break

        # --- Posterior-based assignments per read ---
        # best_idx[j] = index of locus with max posterior for read j
        best_idx = np.argmax(z, axis=1)
        best_prob = z[np.arange(n_reads), best_idx]

        # For alternative locus: set best to -inf and take next max
        z_alt = z.copy()
        z_alt[np.arange(n_reads), best_idx] = -np.inf
        alt_idx = np.argmax(z_alt, axis=1)
        alt_prob = z_alt[np.arange(n_reads), alt_idx]
        alt_prob[alt_prob == -np.inf] = 0.0  # if only one locus has nonzero prob

        posterior_gap = best_prob - alt_prob

        assignments = {}
        for j, read_id in enumerate(read_ids):
            assignments[read_id] = {
                "locus": self.loci[best_idx[j]],
                "prob": float(best_prob[j]),
                "alt_locus": self.loci[alt_idx[j]] if alt_prob[j] > 0 else None,
                "alt_prob": float(alt_prob[j]),
                "gap": float(posterior_gap[j]),
            }

        # Return theta as dict for consistency with original version
        theta_dict = {locus: float(theta[k]) for k, locus in enumerate(self.loci)}
        return theta_dict, assignments

