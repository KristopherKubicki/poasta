use std::time::Instant;
use poasta::io::graph::{load_graph_from_gfa, POAGraphFromGFA};
use poasta::graphs::AlignableRefGraph;

fn main() {
    let POAGraphFromGFA { graph, graph_segments } =
        load_graph_from_gfa::<u32>("tests/test.gfa").expect("load graph");

    let start = Instant::now();
    let max_index = graph.all_nodes().map(|n| n.index()).max().unwrap_or(0) + 1;
    let mut mapping = vec![None; max_index];

    for (segment_ix, _) in graph_segments.names.iter().enumerate() {
        let start_node = graph_segments.start_nodes[segment_ix];
        let end_node = graph_segments.end_nodes[segment_ix];
        mapping[start_node.index()] = Some((segment_ix, 0));
        if start_node == end_node {
            continue;
        }
        let mut curr = start_node;
        let mut pos = 1;
        while let Some(succ) = graph.successors(curr).next() {
            mapping[succ.index()] = Some((segment_ix, pos));
            curr = succ;
            pos += 1;
            if curr == end_node {
                break;
            }
        }
    }

    let elapsed = start.elapsed();
    println!("mapping built in {:?}", elapsed);
    // prevent optimizer
    std::hint::black_box(mapping);
}
