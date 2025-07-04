use poasta::graphs::poa::POAGraph;
use poasta::io::graph::{graph_to_dot, graph_to_gfa};
use poasta::aligner::alignment::AlignedPair;

#[test]
fn poa_graph_output_basic() {
    let mut graph: POAGraph<u16> = POAGraph::new();
    let seq1 = b"AC";
    graph
        .add_alignment_with_weights("seq1", seq1, None, &vec![1; seq1.len()])
        .unwrap();

    let seq_nodes: Vec<_> = graph
        .topological_sorted
        .iter()
        .cloned()
        .filter(|n| *n != graph.start_node() && *n != graph.end_node())
        .collect();

    let alignment = vec![
        AlignedPair::new(Some(seq_nodes[0]), Some(0)),
        AlignedPair::new(None, Some(1)),
    ];
    let seq2 = b"AG";
    graph
        .add_alignment_with_weights("seq2", seq2, Some(&alignment), &vec![1; seq2.len()])
        .unwrap();

    let mut gfa_buf = Vec::new();
    graph_to_gfa(&mut gfa_buf, &graph).unwrap();
    let gfa_output = String::from_utf8(gfa_buf).unwrap();

    assert!(gfa_output.contains("S\ts0\tA"));
    assert!(gfa_output.contains("S\ts1\tC"));
    assert!(gfa_output.contains("S\ts2\tG"));
    assert!(gfa_output.contains("L\ts0\t+\ts1\t+\t0M"));
    assert!(gfa_output.contains("L\ts0\t+\ts2\t+\t0M"));

    let mut dot_buf = Vec::new();
    graph_to_dot(&mut dot_buf, &graph).unwrap();
    let dot_output = String::from_utf8(dot_buf).unwrap();

    assert!(dot_output.contains("digraph"));
    assert!(dot_output.contains("label=\"A\""));
    assert!(dot_output.contains("label=\"C\""));
    assert!(dot_output.contains("label=\"G\""));
}
