use poasta::graphs::poa::POAGraph;
use poasta::io::graph::{graph_to_dot, graph_to_gfa};

#[test]
fn test_graph_export_basic() {
    let mut graph: POAGraph<u32> = POAGraph::new();
    let seq = b"AC";
    let weights = vec![1usize; seq.len()];
    graph
        .add_alignment_with_weights("seq1", seq, None, &weights)
        .unwrap();

    let mut dot_buf = String::new();
    graph_to_dot(&mut dot_buf, &graph).unwrap();
    assert!(dot_buf.contains("2 [label=\"A\""));
    assert!(dot_buf.contains("3 [label=\"C\""));
    assert!(dot_buf.contains("2 -> 3"));

    let mut gfa_buf = Vec::new();
    graph_to_gfa(&mut gfa_buf, &graph).unwrap();
    let gfa_str = String::from_utf8(gfa_buf).unwrap();
    let s_records = gfa_str.lines().filter(|l| l.starts_with('S')).count();
    let l_records = gfa_str.lines().filter(|l| l.starts_with('L')).count();
    assert_eq!(s_records, 1);
    assert_eq!(l_records, 0);
}
