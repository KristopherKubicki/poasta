use poasta::aligner::alignment::{Alignment, AlignedPair};
use poasta::graphs::poa::POAGraph;
use poasta::graphs::AlignableRefGraph;
use poasta::io::gaf::{alignment_to_gaf, NodeSegmentResolver};
use poasta::io::gfa::FieldValue;
use poasta::io::graph::GraphSegments;

#[test]
fn alignment_to_gaf_basic() {
    let mut graph = POAGraph::<u32>::new();
    let weights = vec![1; 2];
    let (s1_start, s1_end) = graph
        .add_nodes_for_sequence(b"AC", &weights, 0, 2)
        .unwrap();
    let (s2_start, s2_end) = graph
        .add_nodes_for_sequence(b"GT", &weights, 0, 2)
        .unwrap();
    graph.add_edge(s1_end, s2_start, 0, 1);

    let mut segments = GraphSegments::default();
    segments.names.push("s1".to_string());
    segments.start_nodes.push(s1_start);
    segments.end_nodes.push(s1_end);
    segments.segment_lengths.push(2);
    segments.names.push("s2".to_string());
    segments.start_nodes.push(s2_start);
    segments.end_nodes.push(s2_end);
    segments.segment_lengths.push(2);

    let resolver = NodeSegmentResolver::new(&graph, &segments);

    let sequence = b"ACGT";
    let alignment: Alignment<_> = vec![
        AlignedPair::new(Some(s1_start), Some(0)),
        AlignedPair::new(Some(s1_end), Some(1)),
        AlignedPair::new(Some(s2_start), Some(2)),
        AlignedPair::new(Some(s2_end), Some(3)),
    ];

    let gaf = alignment_to_gaf(&graph, &segments, "query", sequence, &alignment, &resolver).expect("gaf");
    assert_eq!(gaf.graph_path, ">s1>s2");
    assert_eq!(gaf.query_start, 0);
    assert_eq!(gaf.query_end, 3);
    assert_eq!(gaf.path_aln_start, 0);
    assert_eq!(gaf.path_aln_end, 3);
    assert_eq!(gaf.additional_fields.len(), 1);
    assert_eq!(gaf.additional_fields[0].tag, "cg");
    assert_eq!(gaf.additional_fields[0].value, FieldValue::String("4=".to_string()));
}

#[test]
fn alignment_to_gaf_empty_returns_none() {
    let graph = POAGraph::<u32>::new();
    let segments = GraphSegments::<u32>::default();
    let resolver = NodeSegmentResolver::new(&graph, &segments);
    let alignment: Alignment<_> = Vec::new();
    let result = alignment_to_gaf(&graph, &segments, "empty", b"", &alignment, &resolver);
    assert!(result.is_none());
}
