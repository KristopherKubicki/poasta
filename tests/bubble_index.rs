use poasta::graphs::poa::POAGraph;
use poasta::bubbles::index::BubbleIndex;
use poasta::graphs::AlignableRefGraph;
use poasta::aligner::alignment::{AlignedPair, Alignment};

#[test]
fn simple_bubble() {
    let mut graph: POAGraph<u32> = POAGraph::new();
    let seq1 = b"ABC";
    graph.add_alignment_with_weights("seq1", seq1, None, &vec![1; 3]).unwrap();

    // find nodes for A and C
    let node_a = graph.all_nodes().find(|n| graph.get_symbol(*n) == b'A').unwrap();
    let node_c = graph.all_nodes().find(|n| graph.get_symbol(*n) == b'C').unwrap();

    let seq2 = b"ADC";
    let aln: Alignment<_> = vec![
        AlignedPair::new(Some(node_a), Some(0)),
        AlignedPair::new(None, Some(1)),
        AlignedPair::new(Some(node_c), Some(2)),
    ];
    graph.add_alignment_with_weights("seq2", seq2, Some(&aln), &vec![1; 3]).unwrap();

    let node_b = graph.all_nodes().find(|n| graph.get_symbol(*n) == b'B').unwrap();
    let node_d = graph.all_nodes().find(|n| graph.get_symbol(*n) == b'D').unwrap();

    let index = BubbleIndex::new(&graph);

    assert!(index.is_entrance(node_a));
    assert!(index.is_exit(node_c));
    assert!(!index.is_entrance(node_b));
    assert!(!index.is_exit(node_b));

    let bubbles_a: Vec<_> = index
        .get_node_bubbles(node_a)
        .iter()
        .map(|b| (b.bubble_exit, b.min_dist_to_exit))
        .collect();
    let bubbles_b: Vec<_> = index
        .get_node_bubbles(node_b)
        .iter()
        .map(|b| (b.bubble_exit, b.min_dist_to_exit))
        .collect();
    let bubbles_d: Vec<_> = index
        .get_node_bubbles(node_d)
        .iter()
        .map(|b| (b.bubble_exit, b.min_dist_to_exit))
        .collect();
    let bubbles_c: Vec<_> = index
        .get_node_bubbles(node_c)
        .iter()
        .map(|b| (b.bubble_exit, b.min_dist_to_exit))
        .collect();

    assert_eq!(bubbles_a, vec![(node_c, 2)]);
    assert_eq!(bubbles_b, vec![(node_c, 1)]);
    assert_eq!(bubbles_d, vec![(node_c, 1)]);
    assert_eq!(bubbles_c, vec![(node_c, 0)]);
}
