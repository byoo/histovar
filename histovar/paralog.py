import numpy as np

class EMParalogAssigner:
    def __init__(self, loci, max_iter=100, tol=1e-4):
        self.loci = loci
        self.max_iter = max_iter
        self.tol = tol
        
    def fit_predict(self, reads_info):
        """
        reads_info: list of dicts, where each dict represents a read
        and contains mapping probabilities for each locus:
        [
            {"read_id": "r1", "H3C14": 0.9, "H3C15": 0.1},
            {"read_id": "r2", "H3C14": 0.5, "H3C15": 0.5}, ...
        ]
        """
        if not reads_info:
            return {}, {}

        # Initialize theta (locus abundances) uniformly
        theta = {locus: 1.0 / len(self.loci) for locus in self.loci}
        
        for iteration in range(self.max_iter):
            # E-step: calculate expected assignments z_{jk}
            z = []
            for r in reads_info:
                denom = sum(r.get(locus, 0) * theta[locus] for locus in self.loci)
                if denom > 0:
                    z_j = {locus: (r.get(locus, 0) * theta[locus]) / denom for locus in self.loci}
                else:
                    z_j = {locus: 0 for locus in self.loci}
                z.append(z_j)
            
            # M-step: update theta
            new_theta = {locus: 0.0 for locus in self.loci}
            for z_j in z:
                for locus in self.loci:
                    new_theta[locus] += z_j[locus]
            
            # Normalize new_theta
            total_reads = len(reads_info)
            if total_reads > 0:
                for locus in self.loci:
                    new_theta[locus] /= total_reads
            
            # Check convergence
            diff = sum(abs(new_theta[locus] - theta[locus]) for locus in self.loci)
            theta = new_theta
            if diff < self.tol:
                break
                
        # Calculate posterior gap (confidence) for each read
        assignments = {}
        for i, r in enumerate(reads_info):
            z_j = z[i]
            sorted_probs = sorted(z_j.items(), key=lambda x: x[1], reverse=True)
            best_locus, best_prob = sorted_probs[0]
            if len(sorted_probs) > 1:
                alt_locus, alt_prob = sorted_probs[1]
            else:
                alt_locus, alt_prob = None, 0.0
            posterior_gap = best_prob - alt_prob
            assignments[r["read_id"]] = {
                "locus": best_locus,
                "prob": best_prob,
                "alt_locus": alt_locus,
                "alt_prob": alt_prob,
                "gap": posterior_gap
            }
            
        return theta, assignments

