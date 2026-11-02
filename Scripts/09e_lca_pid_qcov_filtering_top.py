#!/usr/bin/env python3

import sys
from collections import defaultdict

# Rank order from finest to coarsest, with depth = number of lineage
# components to keep and default PID / qcov thresholds.
RANKS = [
    {"name": "species",      "depth": 8, "pidThreshold": 95, "qcovThreshold": 90},
    {"name": "genus",        "depth": 7, "pidThreshold": 92, "qcovThreshold": 0},
    {"name": "family",       "depth": 6, "pidThreshold": 80, "qcovThreshold": 0},
    {"name": "order",        "depth": 5, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "class",        "depth": 4, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "phylum",       "depth": 3, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "kingdom",      "depth": 2, "pidThreshold": 0,  "qcovThreshold": 0},
    {"name": "superkingdom", "depth": 1, "pidThreshold": 0,  "qcovThreshold": 0},
]

#updated ranks for VSEARCH fish database
RANKS = [
    {"name": "species",      "depth": 7, "pidThreshold": 95, "qcovThreshold": 90},
    {"name": "genus",        "depth": 5, "pidThreshold": 92, "qcovThreshold": 0},
    {"name": "family",       "depth": 5, "pidThreshold": 80, "qcovThreshold": 0},
    {"name": "order",        "depth": 4, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "class",        "depth": 3, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "phylum",       "depth": 2, "pidThreshold": 70, "qcovThreshold": 0},
    {"name": "kingdom",      "depth": 1, "pidThreshold": 0,  "qcovThreshold": 0}
]

# Rank prefixes from depth 1 (superkingdom) to depth 8 (species).
#Removed sk rank as thats what our database had and it meant everything else was being assigned at the wrong level
PREFIXES = ["sk__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
PREFIXES = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

def parse_threshold_arg(arg):
    """
    Parse:
    species=95,genus=92,family=85
    """
    overrides = {}

    if not arg:
        return overrides

    for pair in arg.split(","):
        if "=" not in pair:
            continue

        name, value = [x.strip() for x in pair.split("=", 1)]

        try:
            overrides[name.lower()] = float(value)
        except ValueError:
            pass

    return overrides


def build_ranks(pid_override_str=None, qcov_override_str=None):
    pid_overrides = parse_threshold_arg(pid_override_str)
    qcov_overrides = parse_threshold_arg(qcov_override_str)

    ranks = []

    for rank in RANKS:
        ranks.append({
            **rank,
            "pidThreshold": pid_overrides.get(
                rank["name"],
                rank["pidThreshold"]
            ),
            "qcovThreshold": qcov_overrides.get(
                rank["name"],
                rank["qcovThreshold"]
            )
        })

    return ranks


def parse_lineage(lineage_str):
    return [
        x.strip()
        for x in lineage_str.split(";")
        if x.strip()
    ]


def lca(paths):
    """
    Lowest Common Ancestor of lineage paths
    """
    if not paths:
        return []

    base = paths[0]
    common_length = len(base)

    for current in paths[1:]:
        j = 0

        while (
            j < common_length
            and j < len(current)
            and base[j] == current[j]
        ):
            j += 1

        common_length = j

        if common_length == 0:
            break

    return base[:common_length]


def pad_lineage(lineage):
    """
    Pad lineage to number of ranks using rank prefixes.
    """
    padded = lineage.copy()

    while len(padded) < len(PREFIXES):
        padded.append(PREFIXES[len(padded)])

    return padded


def ranked_lca(hits, ranks):
    """
    Walk finest -> coarsest ranks.
    First rank with at least one hit passing PID and qcov:
      - keep hits passing thresholds
      - select highest PID
      - compute LCA
      - truncate to rank depth
    """
    for rank in ranks:

        passing = [
            hit for hit in hits
            if hit["pid"] >= rank["pidThreshold"]
            and hit["qcov"] >= rank["qcovThreshold"]
        ]

        if not passing:
            continue

        max_pid = max(hit["pid"] for hit in passing)

        top_hits = [
            hit for hit in passing
            if hit["pid"] == max_pid
        ]

        max_qcov = max(hit["qcov"] for hit in top_hits)

        common = lca([hit["lineage"] for hit in top_hits])

        return {
            "lineage": common[:rank["depth"]],
            "pid": max_pid,
            "qcov": max_qcov
        }

    return {
        "lineage": [],
        "pid": None,
        "qcov": None
    }


def main(input_file, pid_threshold_arg=None, qcov_threshold_arg=None):

    ranks = build_ranks(pid_threshold_arg, qcov_threshold_arg)

    groups = defaultdict(list)

    with open(input_file, "r") as f:
        for line in f:

            line = line.strip()

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                continue

            id_ = fields[0]
            lineage = fields[1]

            pid = float(fields[2]) if len(fields) > 2 else 0
            qcov = float(fields[3]) if len(fields) > 3 else 0

            groups[id_].append({
                "lineage": parse_lineage(lineage),
                "pid": pid,
                "qcov": qcov
            })

    for id_, hits in groups.items():

        result = ranked_lca(hits, ranks)

        padded = pad_lineage(result["lineage"])

        pid_out = result["pid"] if result["pid"] is not None else "NA"
        qcov_out = result["qcov"] if result["qcov"] is not None else "NA"

        print(
            f"{id_}\t"
            f"{';'.join(padded)};"
            f"\t{pid_out}\t{qcov_out}"
        )


if __name__ == "__main__":

    input_file = sys.argv[1]

    pid_overrides = sys.argv[2] if len(sys.argv) > 2 else None
    qcov_overrides = sys.argv[3] if len(sys.argv) > 3 else None

    main(
        input_file,
        pid_overrides,
        qcov_overrides
    )
