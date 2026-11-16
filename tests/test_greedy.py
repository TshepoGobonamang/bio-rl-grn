from biorlgrn.demo import make_demo
from biorlgrn.io import read_prepared
from biorlgrn.ridge import GraphRidge
from biorlgrn.search import greedy_rank


def test_greedy_selects_and_logs(tmp_path):
    p = make_demo(tmp_path / "toy")
    X, tfs = read_prepared(p)
    engine = GraphRidge(X.iloc[:70], X.iloc[70:], alpha=0.1)
    ranking, trace = greedy_rank(engine, tfs, max_edges=4)
    assert len(ranking) == 12
    assert len(trace) <= 4
    assert trace.empty or trace.utility_gain.min() > 0
    assert ranking.score.is_monotonic_decreasing
