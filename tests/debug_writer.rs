use poasta::debug::{DebugOutputWriter, messages::DebugOutputMessage};
use tempfile::tempdir;

#[test]
fn debug_writer_creates_files() {
    let dir = tempdir().unwrap();
    let writer = DebugOutputWriter::init(dir.path());

    writer.log(DebugOutputMessage::NewSequence {
        seq_name: "seq1".to_string(),
        sequence: "ACGT".to_string(),
        max_rank: 1,
    });

    writer.log(DebugOutputMessage::IntermediateGraph {
        graph_dot: "digraph {}".to_string(),
    });

    writer.log(DebugOutputMessage::AstarData {
        visited_tsv: "0\t0\n".to_string(),
    });

    writer.log(DebugOutputMessage::Terminate);
    writer.join().unwrap();

    let graph_file = dir.path().join("graph_for_seq1.dot");
    assert!(graph_file.exists());

    let tsv_file = dir.path().join("astar_iterations/seq1.iter0.tsv");
    assert!(tsv_file.exists());
}
