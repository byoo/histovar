import pytest
from histovar.paralog import EMParalogAssigner
from histovar.reference import ReferenceBuilder, format_codon_aa_change
from histovar.ml import MLScorer

def test_em_algorithm():
    assigner = EMParalogAssigner(loci=["H3C14", "H3C15"])
    reads = [
        {"read_id": "r1", "H3C14": 0.9, "H3C15": 0.1},
        {"read_id": "r2", "H3C14": 0.8, "H3C15": 0.2},
        {"read_id": "r3", "H3C14": 0.4, "H3C15": 0.6},
    ]
    theta, assignments = assigner.fit_predict(reads)
    
    assert "H3C14" in theta
    assert "H3C15" in theta
    assert abs(sum(theta.values()) - 1.0) < 1e-5
    
    # Check that r1 is strongly assigned to H3C14
    assert assignments["r1"]["locus"] == "H3C14"
    assert assignments["r1"]["prob"] > 0.8
    assert assignments["r1"]["gap"] > 0.6

def test_reference_builder():
    ref = ReferenceBuilder(build="hg38")
    assert ref.get_canonical_name("H3F3A") == "H3-3A"
    assert ref.get_canonical_name("HIST1H3B") == "H3C2"
    
    import os
    with open("dummy.gtf", "w") as f:
        f.write('chr1\tHAVANA\tCDS\t1000\t1500\t.\t+\t0\tgene_name "H3F3A"; transcript_id "ENST001";\n')
        
    ref.parse_gtf("dummy.gtf")
    coords = ref.get_hotspot_coordinates("H3-3A", "K27")
    assert coords is not None
    assert coords["chrom"] == "chr1"
    
    os.remove("dummy.gtf")

def test_deep_learning_scorer():
    scorer = MLScorer(model_type="dl")
    
    # Mock data
    X = [
        {"depth": 100, "af": 0.5, "mapq": 60, "strand_bias": 0.5, "pos_bias": 0.5, "posterior_gap": 0.9},
        {"depth": 20, "af": 0.1, "mapq": 10, "strand_bias": 0.9, "pos_bias": 0.1, "posterior_gap": 0.1}
    ]
    y = [1, 0]
    
    # Test training loop
    scorer.train(X, y)
    assert scorer.is_trained
    
    # Test inference
    probas = scorer.predict_proba(X)
    assert len(probas) == 2
    assert 0.0 <= probas[0] <= 1.0
    assert 0.0 <= probas[1] <= 1.0

def test_codon_aa_translation():
    # Plus strand: AAG (K) -> ATG (M), AGG (R), AAA (K, synonymous), AAT (N)
    assert format_codon_aa_change(["A", "A", "G"], "+", 1, "T", 28, "K") == "K28M"
    assert format_codon_aa_change(["A", "A", "G"], "+", 1, "G", 28, "K") == "K28R"
    assert format_codon_aa_change(["A", "A", "G"], "+", 2, "A", 28, "K") == "K28K"
    assert format_codon_aa_change(["A", "A", "G"], "+", 2, "T", 28, "K") == "K28N"

    # Minus strand: coding AAG (K) has genomic bases ['T', 'T', 'C'] at codon_coords [pos0, pos1, pos2]
    # Mutating middle genomic base T -> A gives coding ATG (M)
    assert format_codon_aa_change(["T", "T", "C"], "-", 1, "A", 28, "K") == "K28M"
    # Mutating third coding base (codon_idx=2, genomic C -> A) gives coding AAT (N) -> K28N, not K28M
    assert format_codon_aa_change(["T", "T", "C"], "-", 2, "A", 28, "K") == "K28N"

    # Indels
    assert format_codon_aa_change(["A", "A", "G"], "+", 1, "A+TG", 28, "K") == "K28ins"
    assert format_codon_aa_change(["A", "A", "G"], "+", 1, "A-3", 28, "K") == "K28del"
