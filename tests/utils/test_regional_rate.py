"""Tests for the regional mutation rate (rr) functions in ``gnomad_constraint.utils.regional_rate``."""

import hail as hl
import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from gnomad_constraint.resources.constants import (
    ADJ_FREQ_META,
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_SCALES,
    METHYLATION_CORRELATED_FEATURES,
    RR_CONTEXTS,
)
from gnomad_constraint.utils.regional_rate import (
    choose_rr_penalty,
    compute_rr_ht,
    compute_rr_normalization_factors,
    create_rr_training_ht,
    fit_rr_context_model,
    fit_rr_models,
    get_chrom_fold,
    get_rr_model_features,
    get_syn_mane_select_ht,
    get_syn_sites_ht,
    normalize_rr_ht,
    remove_qc_flagged_genes,
    tune_rr_penalty,
)

INTERVAL_TYPE = hl.tinterval(hl.tlocus("GRCh38"))


def _locus(pos: int) -> hl.Locus:
    """Return a GRCh38 chr1 locus."""
    return hl.Locus("chr1", pos, "GRCh38")


def _window(start: int, end: int) -> hl.Interval:
    """Return the chr1 interval ``[start, end)``."""
    return hl.Interval(_locus(start), _locus(end))


def _training_df(n: int = 4000, seed: int = 0) -> pd.DataFrame:
    """Return training rows whose labels follow a logistic model of features a and b."""
    rng = np.random.default_rng(seed)
    a = rng.normal(size=n)
    b = 3 * rng.normal(size=n) + 5
    y = rng.random(n) < expit(-2 + 0.8 * a - 0.5 * (b - 5) / 3)
    return pd.DataFrame({"y": y.astype("int64"), "a": a, "b": b})


class TestGetRrModelFeatures:
    """Test the get_rr_model_features function."""

    def test_non_cpg_context(self) -> None:
        """Test that a non-CpG context uses every feature at every scale."""
        assert sorted(get_rr_model_features("AAA")) == sorted(GENOMIC_FEATURE_COLUMNS)

    def test_cpg_context(self) -> None:
        """Test that CpG contexts leave out the methylation-correlated features."""
        dropped = {
            f"{f}_{s}"
            for f in METHYLATION_CORRELATED_FEATURES
            for s in GENOMIC_FEATURE_SCALES
        }
        features = get_rr_model_features("ACG")
        assert set(features) == set(GENOMIC_FEATURE_COLUMNS) - dropped
        assert len(features) == len(GENOMIC_FEATURE_COLUMNS) - len(dropped)


class TestFitRrContextModel:
    """Test the fit_rr_context_model function."""

    def test_minimizes_the_penalized_log_loss(self) -> None:
        """
        Test that the fit minimizes mean log loss + penalty / 2 * sum(coef ** 2).

        The gradient of that objective, with the intercept unpenalized and the
        features standardized with the stored ``mean`` and ``sd``, is 0 at the fit.
        This checks the per-site penalty passed to sklearn and that the stored
        standardization is the one the model was fit with.
        """
        df = _training_df()
        penalty = 0.3
        model = fit_rr_context_model(df, ["a", "b"], ridge_penalty=penalty)
        z = (df[["a", "b"]].to_numpy() - model["mean"]) / model["sd"]
        coefs = np.array(model["coefs"])
        residual = expit(model["intercept"] + z @ coefs) - df["y"].to_numpy()
        assert residual.mean() == pytest.approx(0, abs=1e-4)
        assert z.T @ residual / len(df) + penalty * coefs == pytest.approx(
            [0, 0], abs=1e-4
        )
        assert model["sd"] == pytest.approx(list(df[["a", "b"]].std(ddof=0)))

    def test_ave(self) -> None:
        """Test that ave is the predicted probability with every feature at its mean."""
        model = fit_rr_context_model(_training_df(), ["a", "b"], 0.3)
        assert model["ave"] == pytest.approx(expit(model["intercept"]))

    def test_constant_feature_is_left_out(self) -> None:
        """
        Test that a feature with a single value is left out of the model.

        Regression test: z-scoring a constant feature gave NaN, so
        ``--train-rr-models --test`` crashed in contexts where a feature was constant
        on chr20.
        """
        df = _training_df().assign(c=0.0)
        model = fit_rr_context_model(df, ["a", "c", "b"], 0.3)
        assert model["features"] == ["a", "b"]
        assert len(model["mean"]) == len(model["sd"]) == len(model["coefs"]) == 2
        assert np.isfinite(model["coefs"]).all()

    def test_rows_with_missing_features_are_dropped(self) -> None:
        """Test that sites missing any feature aren't counted or used."""
        df = _training_df()
        df.loc[:9, "a"] = np.nan
        model = fit_rr_context_model(df, ["a", "b"], 0.3)
        assert model["n_dnm"] + model["n_background"] == len(df) - 10
        assert model["n_dnm"] == int(df["y"][10:].sum())

    @pytest.mark.parametrize("case", ["one_class", "no_varying_feature"])
    def test_no_model(self, case: str) -> None:
        """Test that there's no model without both classes or a varying feature."""
        df = _training_df()
        if case == "one_class":
            df = df[df["y"] == 0]
        else:
            df = df.assign(a=1.0, b=2.0)
        assert fit_rr_context_model(df, ["a", "b"], 0.3) is None


class TestFitRrModels:
    """Test the fit_rr_models function."""

    def test_models_per_context(self) -> None:
        """Test that each context with training data gets a model with its features."""
        rng = np.random.default_rng(1)
        n = 600
        df = pd.DataFrame(rng.normal(size=(n, len(GENOMIC_FEATURE_COLUMNS))))
        df.columns = GENOMIC_FEATURE_COLUMNS
        df["y"] = (rng.random(n) < 0.2).astype("int64")
        df["context"] = ["AAA"] * (n // 2) + ["ACG"] * (n // 2)
        df.loc[df["context"] == "ACG", "dist2telo_1k"] = 0.0
        models = {m["context"]: m for m in fit_rr_models(df, 0.3, contexts=RR_CONTEXTS)}
        assert set(models) == {"AAA", "ACG"}
        assert models["AAA"]["features"] == get_rr_model_features("AAA")
        assert models["ACG"]["features"] == [
            f for f in get_rr_model_features("ACG") if f != "dist2telo_1k"
        ]


class TestGetChromFold:
    """Test the get_chrom_fold function."""

    @pytest.mark.parametrize(
        "contig, fold",
        [("chr1", "odd"), ("chr2", "even"), ("chr21", "odd"), ("chrX", None)],
    )
    def test_fold(self, contig: str, fold: str) -> None:
        """Test that autosomes are split by parity and other contigs have no fold."""
        assert get_chrom_fold(contig) == fold


class TestTuneRrPenalty:
    """Test tune_rr_penalty and choose_rr_penalty."""

    @staticmethod
    def _df() -> pd.DataFrame:
        """
        Return training rows of one context with a real effect, on chr1-chr4.

        Two features carry the effect; the others are constant, so the models leave
        them out.
        """
        df = _training_df(n=8000, seed=2)
        df = df.rename(columns={"a": "dist2telo_1k", "b": "GC_content_1k"})
        for column in GENOMIC_FEATURE_COLUMNS:
            if column not in df:
                df[column] = 0.0
        df["contig"] = [f"chr{1 + i % 4}" for i in range(len(df))]
        df["context"] = "AAA"
        return df

    def test_held_out_log_loss(self) -> None:
        """
        Test that each fold's model is scored on the other fold's sites.

        The summed log loss of the odd fold matches that of a model fit on the even
        chromosomes, computed by hand.
        """
        df = self._df()
        tuning_df = tune_rr_penalty(df, [0.3], contexts=["AAA"])
        assert sorted(tuning_df["fold"]) == ["even", "odd"]
        odd = tuning_df.set_index("fold").loc["odd"]
        is_odd = df["contig"].isin(["chr1", "chr3"])
        assert odd["n_test"] == is_odd.sum()

        model = fit_rr_context_model(df[~is_odd], get_rr_model_features("AAA"), 0.3)
        assert model["features"] == ["dist2telo_1k", "GC_content_1k"]
        test_df = df[is_odd]
        z = (test_df[model["features"]].to_numpy() - model["mean"]) / model["sd"]
        p = expit(model["intercept"] + z @ np.array(model["coefs"]))
        y = test_df["y"].to_numpy()
        expected = -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
        assert odd["log_loss"] == pytest.approx(expected, rel=1e-6)

    def test_choose_penalty(self) -> None:
        """Test that the penalty with the lowest summed held-out loss is chosen."""
        tuning_df = pd.DataFrame(
            {
                "penalty": [0.1, 0.1, 0.3, 0.3],
                "context": ["AAA", "ACG", "AAA", "ACG"],
                "fold": ["odd", "odd", "odd", "odd"],
                "n_test": [10, 10, 10, 10],
                "log_loss": [5.0, 3.0, 4.0, 3.5],
            }
        )
        assert choose_rr_penalty(tuning_df) == 0.3

    def test_no_held_out_sites(self) -> None:
        """Test that the default penalty is used when there's only one chromosome."""
        df = self._df().assign(contig="chr20")
        tuning_df = tune_rr_penalty(df, [0.01, 1000.0], contexts=["AAA"])
        assert tuning_df.empty
        assert choose_rr_penalty(tuning_df, default=0.3) == 0.3

    def test_too_strong_penalty_loses(self) -> None:
        """Test that a penalty that flattens a real effect has a higher held-out loss."""
        tuning_df = tune_rr_penalty(self._df(), [0.01, 1000.0], contexts=["AAA"])
        assert choose_rr_penalty(tuning_df) == 0.01


class TestComputeRrHt:
    """Test the compute_rr_ht function."""

    MODELS = [
        hl.Struct(
            context="AAA",
            features=["f1"],
            mean=[0.0],
            sd=[1.0],
            coefs=[1.0],
            intercept=0.0,
            ave=0.5,
        ),
        hl.Struct(
            context="ACG",
            features=["f2"],
            mean=[0.0],
            sd=[2.0],
            coefs=[-1.0],
            intercept=0.0,
            ave=0.5,
        ),
    ]

    def test_rr(self) -> None:
        """
        Test rr per window, the 1.0 fill, and which windows are kept.

        rr = sigmoid(eta) / ave. A window missing a context's feature gets 1.0 for
        that context, contexts without a model get 1.0, and a window with no rr for
        any modeled context is removed.
        """
        features_ht = hl.Table.parallelize(
            [
                {"interval": _window(1, 1001), "f1": 0.0, "f2": 2.0},
                {"interval": _window(1001, 2001), "f1": 1.0, "f2": None},
                {"interval": _window(2001, 3001), "f1": None, "f2": None},
            ],
            hl.tstruct(interval=INTERVAL_TYPE, f1=hl.tfloat64, f2=hl.tfloat64),
            key="interval",
        )
        rows = {
            r.interval.start.position: r.adj_r
            for r in compute_rr_ht(features_ht, self.MODELS).collect()
        }
        assert set(rows) == {1, 1001}
        assert set(rows[1]) == set(RR_CONTEXTS)
        assert rows[1]["AAA"] == pytest.approx(1.0)
        assert rows[1]["ACG"] == pytest.approx(expit(-1.0) / 0.5)
        assert rows[1]["CCC"] == 1.0
        assert rows[1001]["AAA"] == pytest.approx(expit(1.0) / 0.5)
        assert rows[1001]["ACG"] == 1.0


class TestNormalizeRr:
    """Test compute_rr_normalization_factors and normalize_rr_ht."""

    @staticmethod
    def _rr_ht() -> hl.Table:
        """Return two windows with rr for AAA, ACG, and CCC (no synonymous sites)."""
        return hl.Table.parallelize(
            [
                {
                    "interval": _window(1, 100),
                    "adj_r": {"AAA": 2.0, "ACG": 0.5, "CCC": 3.0},
                },
                {
                    "interval": _window(100, 200),
                    "adj_r": {"AAA": 1.0, "ACG": 1.0, "CCC": 1.0},
                },
            ],
            hl.tstruct(interval=INTERVAL_TYPE, adj_r=hl.tdict(hl.tstr, hl.tfloat64)),
            key="interval",
        )

    @staticmethod
    def _syn_ht() -> hl.Table:
        """Return synonymous sites, including one outside the windows and one with no exp."""
        ht = hl.Table.parallelize(
            [
                {"pos": 10, "context": "AAA", "exp": 1.0},
                {"pos": 150, "context": "AAA", "exp": 3.0},
                {"pos": 20, "context": "ACG", "exp": 2.0},
                {"pos": 500, "context": "AAA", "exp": 10.0},
                {"pos": 30, "context": "AAA", "exp": None},
            ],
            hl.tstruct(pos=hl.tint32, context=hl.tstr, exp=hl.tfloat64),
        )
        return ht.key_by(
            locus=hl.locus("chr1", ht.pos, "GRCh38"), alleles=hl.literal(["A", "G"])
        )

    def test_factors(self) -> None:
        """Test that a factor is the exp-weighted mean rr of the sites in a window."""
        factors = compute_rr_normalization_factors(self._syn_ht(), self._rr_ht())
        assert factors == pytest.approx({"AAA": (1 * 2 + 3 * 1) / 4, "ACG": 0.5})

    def test_normalized_mean_is_one(self) -> None:
        """Test that normalized rr averages 1, and contexts without a factor don't change."""
        syn_ht = self._syn_ht()
        rr_ht = self._rr_ht()
        rr_ht = normalize_rr_ht(rr_ht, compute_rr_normalization_factors(syn_ht, rr_ht))
        rr = {r.interval.start.position: r.adj_r for r in rr_ht.collect()}
        assert (1 * rr[1]["AAA"] + 3 * rr[100]["AAA"]) / 4 == pytest.approx(1.0)
        assert rr[1]["ACG"] == pytest.approx(1.0)
        assert rr[1]["CCC"] == 3.0


class TestGetSynManeSelectHt:
    """Test the get_syn_mane_select_ht function."""

    def test_rows_and_exp(self) -> None:
        """Test that only synonymous Ensembl MANE Select rows are kept, with adj exp."""
        rows = [
            ("synonymous_variant", True, "ENST1", 1),
            ("synonymous_variant", True, "NM_1", 2),
            ("synonymous_variant", False, "ENST2", 3),
            ("missense_variant", True, "ENST1", 4),
        ]
        ht = hl.Table.parallelize(
            [
                {
                    "pos": pos,
                    "annotation": csq,
                    "mane_select": mane,
                    "transcript": transcript,
                    "context": "AAA",
                    "expected_variants": [0.1 * pos, float(pos)],
                }
                for csq, mane, transcript, pos in rows
            ],
            hl.tstruct(
                pos=hl.tint32,
                annotation=hl.tstr,
                mane_select=hl.tbool,
                transcript=hl.tstr,
                context=hl.tstr,
                expected_variants=hl.tarray(hl.tfloat64),
            ),
        )
        ht = ht.key_by(
            locus=hl.locus("chr1", ht.pos, "GRCh38"), alleles=hl.literal(["A", "G"])
        )
        ht = ht.annotate_globals(
            exomes_freq_meta=[{"group": "adj", "downsampling": "10"}, ADJ_FREQ_META]
        )
        out = get_syn_mane_select_ht(ht).collect()
        assert [(r.locus.position, r.context, r.exp) for r in out] == [(1, "AAA", 1.0)]


class TestCreateRrTrainingHt:
    """Test get_syn_sites_ht and create_rr_training_ht."""

    def test_training_sites(self) -> None:
        """
        Test the DNM and background sites of the training set.

        DNMs are the synonymous rows of the DNM Table. The background (all of it, with
        a fraction of 1) is the synonymous sites that aren't a DNM and aren't in the
        exomes sites Table, one per locus. Every site gets its window's features.
        """
        context_rows = [
            (10, ["A", "G"], "AAA", "synonymous_variant"),
            (20, ["C", "T"], "ACG", "synonymous_variant"),
            (30, ["A", "G"], "AAA", "synonymous_variant"),
            (30, ["A", "C"], "AAA", "synonymous_variant"),
            (40, ["A", "G"], "AAA", "synonymous_variant"),
            (50, ["C", "T"], "ACG", "missense_variant"),
        ]
        context_ht = hl.Table.parallelize(
            [
                {
                    "locus": _locus(pos),
                    "alleles": alleles,
                    "context": context,
                    "ref": alleles[0],
                    "alt": alleles[1],
                    "methylation_level": 0,
                    "vep": {"most_severe_consequence": csq},
                }
                for pos, alleles, context, csq in context_rows
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                vep=hl.tstruct(most_severe_consequence=hl.tstr),
            ),
            key=["locus", "alleles"],
        )
        dnm_ht = context_ht.filter(
            hl.literal({10, 50}).contains(context_ht.locus.position)
        )
        dnm_ht = dnm_ht.select(
            "context", most_severe_consequence=dnm_ht.vep.most_severe_consequence
        )
        exomes_ht = context_ht.filter(context_ht.locus.position == 40).select()
        features_ht = hl.Table.parallelize(
            [
                {
                    "interval": _window(1, 1000),
                    **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 1.0),
                }
            ],
            hl.tstruct(
                interval=INTERVAL_TYPE,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
            key="interval",
        )
        ht = create_rr_training_ht(
            get_syn_sites_ht(context_ht), dnm_ht, exomes_ht, features_ht, 1.0, 0
        )
        rows = ht.collect()
        assert sorted((r.locus.position, r.y) for r in rows) == [
            (10, 1),
            (20, 0),
            (30, 0),
        ]
        assert all(r.dist2telo_1k == 1.0 for r in rows)


class TestRemoveQcFlaggedGenes:
    """Test the remove_qc_flagged_genes function."""

    def test_flagged_ensembl_mane_select_spans(self) -> None:
        """Test that only the spans of flagged Ensembl MANE Select transcripts are removed."""
        ht = hl.Table.parallelize(
            [{"locus": _locus(p)} for p in [100, 200, 300]],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
            key="locus",
        )
        metrics = [
            ("ENST1", True, {"low_exome_coverage"}, 150, 250),
            ("ENST2", True, set(), 50, 120),
            ("NM_1", True, {"low_exome_coverage"}, 280, 320),
            ("ENST3", False, {"low_exome_coverage"}, 90, 110),
        ]
        metrics_ht = hl.Table.parallelize(
            [
                {
                    "transcript": transcript,
                    "mane_select": mane_select,
                    "gene_flags": flags,
                    "chromosome": "chr1",
                    "start_position": start,
                    "end_position": end,
                }
                for transcript, mane_select, flags, start, end in metrics
            ],
            hl.tstruct(
                transcript=hl.tstr,
                mane_select=hl.tbool,
                gene_flags=hl.tset(hl.tstr),
                chromosome=hl.tstr,
                start_position=hl.tint32,
                end_position=hl.tint32,
            ),
        )
        out = remove_qc_flagged_genes(ht, metrics_ht).collect()
        assert [r.locus.position for r in out] == [100, 300]
