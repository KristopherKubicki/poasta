use poasta::aligner::alignment::{AlignedPair, Alignment};
use poasta::graphs::poa::POAGraph;
use poasta::io::graph::{graph_to_dot, graph_to_gfa};

fn build_branching_graph() -> (POAGraph<u16>, usize, usize, usize) {
    let mut graph: POAGraph<u16> = POAGraph::new();
    let seq1 = b"AC";
    graph
        .add_alignment_with_weights("s1", seq1, None, &vec![1; seq1.len()])
        .unwrap();
    let node_a = graph
        .all_nodes()
        .find(|n| graph.get_symbol(*n) == b'A')
        .unwrap();
    let seq2 = b"AG";
    let aln: Alignment<_> = vec![
        AlignedPair::new(Some(node_a), Some(0)),
        AlignedPair::new(None, Some(1)),
    ];
    graph
        .add_alignment_with_weights("s2", seq2, Some(&aln), &vec![1; seq2.len()])
        .unwrap();

    let node_c = graph
        .all_nodes()
        .find(|n| graph.get_symbol(*n) == b'C')
        .unwrap();
    let node_g = graph
        .all_nodes()
        .find(|n| graph.get_symbol(*n) == b'G')
        .unwrap();

    (
        graph,
        node_a.index(),
        node_c.index(),
        node_g.index(),
    )
}

#[test]
fn graph_to_gfa_basic() {
    let (graph, _, _, _) = build_branching_graph();

    let mut out = Vec::new();
    graph_to_gfa(&mut out, &graph).unwrap();
    let text = String::from_utf8(out).unwrap();

    assert!(text.contains("H\tVN:Z:1.1"));
    assert!(text.contains("S\ts0\tA"));
    assert!(text.contains("S\ts1\tC"));
    assert!(text.contains("S\ts2\tG"));
}

#[test]
fn graph_to_dot_basic() {
    let (graph, a_ix, c_ix, g_ix) = build_branching_graph();

    let mut out = String::new();
    graph_to_dot(&mut out, &graph).unwrap();

    assert!(out.contains("digraph"));
    assert!(out.contains("[label=\"A\""));
    assert!(out.contains("[label=\"C\""));
    assert!(out.contains("[label=\"G\""));
    assert!(out.contains(&format!("{} -> {}", a_ix, c_ix)));
    assert!(out.contains(&format!("{} -> {}", a_ix, g_ix)));
}
