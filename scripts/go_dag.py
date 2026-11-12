"""Gene Ontology graph from an OBO file: parents, ancestors, depth, alternate IDs.

Obsolete terms stay in parse_obo's output but are not part of a Dag. A term can reach
its root by paths of different lengths, so there are two depths; edges to obsolete or
missing terms are ignored. Standard library only.

    terms, alt = parse_obo("go-basic.obo")
    dag = Dag(terms, alt, RELATION_SETS["is_a+part_of"])
    dag.ancestors(dag.resolve("GO:0044763"))
"""

RELATION_SETS = {
    "is_a": {"is_a"},
    "is_a+part_of": {"is_a", "part_of"},
    "is_a+part_of+regulates": {"is_a", "part_of", "regulates",
                               "positively_regulates", "negatively_regulates"},
}
DOMAINS = {"molecular_function": "F", "biological_process": "P", "cellular_component": "C"}


def obo_version(path):
    """The OBO file's data-version line, e.g. "releases/2026-07-26"."""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("data-version:"):
                return line.split(": ", 1)[1].strip()
            if line.startswith("["):
                break
    return "unknown"


def parse_obo(path):
    """Return ({id: {name, domain, obsolete, subsets, edges: [(relation, parent)]}}, {alt_id: id})."""
    terms, alt, cur = {}, {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("["):
                cur = {"name": None, "domain": None, "obsolete": False, "subsets": set(), "edges": []} if line == "[Term]" else None
            elif cur is not None and ": " in line:
                key, value = line.split(": ", 1)
                if key == "id":
                    terms[value], cur["id"] = cur, value
                elif key == "name":
                    cur["name"] = value
                elif key == "namespace":
                    cur["domain"] = DOMAINS[value]
                elif key == "is_obsolete":
                    cur["obsolete"] = value == "true"
                elif key == "subset":
                    cur["subsets"].add(value)
                elif key == "alt_id":
                    alt[value] = cur["id"]
                elif key == "is_a":
                    cur["edges"].append(("is_a", value.split()[0]))
                elif key == "relationship":
                    relation, parent = value.split()[:2]
                    cur["edges"].append((relation, parent))
    return terms, alt


class Dag:
    """The live GO graph restricted to one set of relations."""

    def __init__(self, terms, alt, relations):
        self.alt = alt
        self.live = {t for t, v in terms.items() if not v["obsolete"]}
        self.parents = {t: {p for r, p in terms[t]["edges"] if r in relations and p in self.live}
                        for t in self.live}
        self.has_child = {p for ps in self.parents.values() for p in ps}
        self.roots = {t for t, ps in self.parents.items() if not ps}
        self._ancestors, self._longest, self._shortest = {}, {}, {}

    def resolve(self, term):
        """Map a current or alternate ID to its live ID; None if obsolete or unknown."""
        term = self.alt.get(term, term)
        return term if term in self.live else None

    def ancestors(self, term):
        if term not in self._ancestors:
            found = set()
            for parent in self.parents[term]:
                found.add(parent)
                found |= self.ancestors(parent)
            self._ancestors[term] = frozenset(found)
        return self._ancestors[term]

    def longest_depth(self, term):
        """Edges on the longest path to a root."""
        if term not in self._longest:
            self._longest[term] = 1 + max(map(self.longest_depth, self.parents[term])) if self.parents[term] else 0
        return self._longest[term]

    def shortest_depth(self, term):
        """Edges on the shortest path to a root."""
        if term not in self._shortest:
            self._shortest[term] = 1 + min(map(self.shortest_depth, self.parents[term])) if self.parents[term] else 0
        return self._shortest[term]

    def is_leaf(self, term):
        return term not in self.has_child

    def at_level(self, term, level, depth="longest"):
        """Ancestors exactly `level` edges from the root; the term itself if it is no deeper.

        A term can have several, because paths through the DAG branch.
        """
        measure = self.longest_depth if depth == "longest" else self.shortest_depth
        if measure(term) <= level:
            return {term}
        return {a for a in self.ancestors(term) if measure(a) == level}

    def slim_terms(self, terms, name):
        """The live terms of a GO slim, without the domain roots, which would cover everything."""
        found = {t for t in self.live - self.roots if name in terms[t]["subsets"]}
        if not found:
            raise ValueError(f"no terms in subset {name}")
        return found

    def to_slim(self, term, slim):
        """The most specific slim terms at or above a term; empty if the slim doesn't cover it."""
        found = (self.ancestors(term) | {term}) & slim
        return found - set().union(*(self.ancestors(t) for t in found))
