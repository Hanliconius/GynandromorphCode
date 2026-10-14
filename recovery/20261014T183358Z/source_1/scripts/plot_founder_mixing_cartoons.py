from pathlib import Path
from itertools import product
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from scipy.stats import betabinom

OUT = Path("analysis/figures/founder_mixing_cartoons.pdf")
OUT.parent.mkdir(parents=True, exist_ok=True)

DEPTH = 40
CONCENTRATION = 22
ERROR = 0.001
ALLELE_FREQUENCY = 0.5
STATES = list(product((0., .5, 1.), repeat=2))
MIXES = [(0., "#555555"), (.2, "#0072B2"), (.4, "#D55E00")]
K = np.arange(DEPTH + 1)
AF = K / DEPTH

plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})

def allele_probability(bits):
    return np.prod([
        ALLELE_FREQUENCY if b else 1 - ALLELE_FREQUENCY
        for b in bits
    ])

def state_weights(kind):
    weights = {state: 0. for state in STATES}

    if kind == "three":
        # One shared haplotype: Left=A+B; Right=A+C.
        for a, b, c in product((0, 1), repeat=3):
            state = ((a + b) / 2, (a + c) / 2)
            weights[state] += allele_probability((a, b, c))

    elif kind == "four":
        # Four independent haplotypes: Left=A+B; Right=C+D.
        for a, b, c, d in product((0, 1), repeat=4):
            state = ((a + b) / 2, (c + d) / 2)
            weights[state] += allele_probability((a, b, c, d))

    elif kind == "siblings":
        # Two unrelated diploid parents; independent transmissions
        # to two offspring at this locus.
        for parents in product((0, 1), repeat=4):
            parent_probability = allele_probability(parents)
            for i, j, k, m in product((0, 1), repeat=4):
                left = (parents[i] + parents[2 + j]) / 2
                right = (parents[k] + parents[2 + m]) / 2
                weights[(left, right)] += parent_probability / 16
    else:
        raise ValueError(kind)

    assert np.isclose(sum(weights.values()), 1)
    return weights

def pmf(q):
    mean = ERROR + (1 - 2 * ERROR) * q
    values = betabinom.pmf(
        K, DEPTH,
        mean * CONCENTRATION,
        (1 - mean) * CONCENTRATION,
    )
    return values / values.sum()

def component(state, mixing):
    left, right = state
    ql = (1 - mixing) * left + mixing * right
    qr = mixing * left + (1 - mixing) * right
    return ql, qr, np.outer(pmf(qr), pmf(ql))

def threshold(density, mass=.5):
    ordered = np.sort(density.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered) - 1)]

def style(axis, title):
    axis.set(
        xlim=(-.025, 1.025), ylim=(-.025, 1.025),
        aspect="equal", title=title,
        xlabel="Left ALT-read fraction",
        ylabel="Right ALT-read fraction",
    )
    axis.plot([0, 1], [0, 1], color="#DDDDDD", lw=.8, zorder=0)
    axis.spines[["top", "right"]].set_visible(False)

def draw(axis, kind, shaded):
    weights = state_weights(kind)

    for mixing, color in MIXES:
        density = np.zeros((DEPTH + 1, DEPTH + 1))

        for state, weight in weights.items():
            if weight == 0:
                continue
            ql, qr, conditional = component(state, mixing)
            density += weight * conditional

            if not shaded:
                # Stationary contours need drawing only once.
                if state[0] != state[1] or mixing == 0:
                    level = threshold(conditional)
                    axis.contour(
                        AF, AF, conditional, levels=[level],
                        colors=[color], linewidths=1.1, alpha=.75,
                    )

            if state[0] != state[1] or mixing == 0:
                axis.scatter(
                    ql, qr,
                    s=(450 * weight if shaded else 28),
                    color=color, edgecolors="white",
                    linewidths=.5, zorder=5,
                )

        assert np.isclose(density.sum(), 1)

        if shaded:
            level = threshold(density)
            axis.contourf(
                AF, AF, density,
                levels=[level, density.max() + 1e-12],
                colors=[color], alpha=.13,
            )
            axis.contour(
                AF, AF, density, levels=[level],
                colors=[color], linewidths=1.2,
            )

def page(pdf, title, models, shaded, explanation):
    fig, axes = plt.subplots(1, 2, figsize=(12, 7.5))
    fig.subplots_adjust(
        left=.075, right=.975, bottom=.25, top=.82, wspace=.23
    )

    for axis, (kind, label) in zip(axes, models):
        style(axis, label)
        draw(axis, kind, shaded)

    fig.suptitle(title, fontsize=19, weight="bold", y=.95)
    fig.legend(
        handles=[
            Line2D([0], [0], color=color, marker="o",
                   label=f"{mixing:.0%} mixing per side")
            for mixing, color in MIXES
        ],
        loc="upper center", bbox_to_anchor=(.5, .90),
        ncol=3, frameon=False,
    )
    fig.text(.075, .16, explanation, fontsize=10.5, linespacing=1.5)
    fig.text(
        .075, .045,
        "ILLUSTRATIVE: allele frequency 0.5; depth 40/side; "
        "beta-binomial concentration 22; error 0.001; no bias.\n"
        "Contours describe predictive read-fraction regions, not confidence "
        "in fitted centroids. Discrete boundaries are interpolated.\n"
        "No empirical sibling calibration or observed data are used.",
        fontsize=9, color="#555555", linespacing=1.5,
    )
    pdf.savefig(fig)
    plt.close(fig)

with PdfPages(OUT) as pdf:
    page(
        pdf,
        "Founder sharing changes the available genotype states",
        [
            ("three", "Three-founder: one shared haplotype"),
            ("four", "Four-founder: independent haplotypes"),
        ],
        False,
        "Contours: 50% predictive region conditional on each genotype pair.\n"
        "Three-founder inheritance excludes opposing homozygotes; "
        "the other centroid positions are shared between models.",
    )
    page(
        pdf,
        "Full-sibling inheritance redistributes predictive density",
        [
            ("four", "Unrelated diploid pair"),
            ("siblings", "Full siblings: Mendelian inheritance"),
        ],
        True,
        "Shading: 50% highest-probability region of each complete genotype mixture.\n"
        "Point area represents genotype-state probability. Relatedness changes "
        "the weights, while centroid positions and conditional spread stay fixed.",
    )

print(OUT.resolve())
