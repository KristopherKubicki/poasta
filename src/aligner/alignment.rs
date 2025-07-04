use crate::graphs::{AlignableRefGraph, NodeIndexType};

#[derive(Clone, Debug)]
pub struct AlignedPair<N>
where
    N: NodeIndexType
{
    /// Represents the node rank in the graph
    pub rpos: Option<N>,

    /// Query sequence position
    pub qpos: Option<usize>
}

impl<N> AlignedPair<N>
where
    N: NodeIndexType
{
    pub fn new(rpos: Option<N>, qpos: Option<usize>) -> Self {
        Self { rpos, qpos }
    }

    pub fn is_aligned(&self) -> bool {
        matches!((self.rpos, self.qpos), (Some(_), Some(_)))
    }

    pub fn is_indel(&self) -> bool {
        !self.is_aligned()
    }
    
    pub fn is_deletion(&self) -> bool {
        self.rpos.is_none() && self.qpos.is_some()
    }
    
    pub fn is_insertion(&self) -> bool {
        self.rpos.is_some() && self.qpos.is_none()
    }
}

pub type Alignment<N> = Vec<AlignedPair<N>>;

pub fn print_alignment<G, N>(graph: &G, sequence: &[u8], aln: &Alignment<N>) -> String
where
    G: AlignableRefGraph<NodeIndex=N>,
    N: NodeIndexType
{
    let mut graph_chars = Vec::with_capacity(aln.len());
    let mut aln_chars = Vec::with_capacity(aln.len());
    let mut query_chars = Vec::with_capacity(aln.len());

    for pair in aln {
        if pair.is_aligned() {
            let node = graph.get_symbol_char(pair.rpos.unwrap());
            let qry = char::from(sequence[pair.qpos.unwrap()]);

            graph_chars.push(node);
            aln_chars.push(if node == qry { '|' } else { '·' });
            query_chars.push(qry);
        } else if let Some(nix) = pair.rpos {
            let node = graph.get_symbol_char(nix);
            graph_chars.push(node);
            aln_chars.push(' ');
            query_chars.push('-');
        } else if let Some(qpos) = pair.qpos {
            let qry = char::from(sequence[qpos]);
            graph_chars.push('-');
            aln_chars.push(' ');
            query_chars.push(qry);
        }
    }

    format!(
        "{}\n{}\n{}",
        graph_chars.into_iter().collect::<String>(),
        aln_chars.into_iter().collect::<String>(),
        query_chars.into_iter().collect::<String>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::algo::toposort;
    use petgraph::graph::{DiGraph, NodeIndex, NodeIndices, Neighbors};
    use petgraph::{Incoming, Outgoing};

    /// Simple directed graph storing a character per node.
    struct CharGraph {
        graph: DiGraph<char, ()>,
        end: NodeIndex,
    }

    impl CharGraph {
        fn linear(chars: &[char]) -> Self {
            let mut g: DiGraph<char, ()> = DiGraph::new();
            let nodes: Vec<_> = chars.iter().map(|&c| g.add_node(c)).collect();
            for i in 0..nodes.len() - 1 {
                g.add_edge(nodes[i], nodes[i + 1], ());
            }
            let end = g.add_node('-');
            g.add_edge(*nodes.last().unwrap(), end, ());
            Self { graph: g, end }
        }
    }

    impl AlignableRefGraph for CharGraph {
        type NodeIndex = NodeIndex;
        type NodeIterator<'a> = NodeIndices where Self: 'a;
        type PredecessorIterator<'a> = Neighbors<'a, (), u32> where Self: 'a;
        type SuccessorIterator<'a> = Neighbors<'a, (), u32> where Self: 'a;

        fn all_nodes(&self) -> Self::NodeIterator<'_> {
            self.graph.node_indices()
        }

        fn node_count(&self) -> usize {
            self.graph.node_count()
        }

        fn node_count_with_start_and_end(&self) -> usize {
            self.graph.node_count()
        }

        fn edge_count(&self) -> usize {
            self.graph.edge_count()
        }

        fn start_node(&self) -> Self::NodeIndex {
            NodeIndex::new(0)
        }

        fn end_node(&self) -> Self::NodeIndex {
            self.end
        }

        fn predecessors(&self, node: Self::NodeIndex) -> Self::PredecessorIterator<'_> {
            self.graph.neighbors_directed(node, Incoming)
        }

        fn successors(&self, node: Self::NodeIndex) -> Self::SuccessorIterator<'_> {
            self.graph.neighbors_directed(node, Outgoing)
        }

        fn in_degree(&self, node: Self::NodeIndex) -> usize {
            self.graph.neighbors_directed(node, Incoming).count()
        }

        fn out_degree(&self, node: Self::NodeIndex) -> usize {
            self.graph.neighbors_directed(node, Outgoing).count()
        }

        fn is_end(&self, node: Self::NodeIndex) -> bool {
            node == self.end
        }

        fn get_symbol_char(&self, node: Self::NodeIndex) -> char {
            self.graph[node]
        }

        fn is_symbol_equal(&self, node: Self::NodeIndex, symbol: u8) -> bool {
            self.graph[node] as u8 == symbol
        }

        fn get_node_ordering(&self) -> Vec<usize> {
            let toposorted = toposort(&self.graph, None).unwrap();
            let mut node_ordering = vec![0; toposorted.len()];
            for (rank, node) in toposorted.iter().enumerate() {
                node_ordering[node.index()] = rank;
            }
            node_ordering
        }
    }

    #[test]
    fn test_aligned_pair_methods() {
        let pair = AlignedPair::new(Some(NodeIndex::<u32>::new(1)), Some(2));
        assert!(pair.is_aligned());
        assert!(!pair.is_indel());

        let ins = AlignedPair::new(Some(NodeIndex::<u32>::new(1)), None);
        assert!(ins.is_indel());
        assert!(ins.is_insertion());
        assert!(!ins.is_deletion());

        let del: AlignedPair<NodeIndex> = AlignedPair::new(None, Some(0));
        assert!(del.is_indel());
        assert!(del.is_deletion());
        assert!(!del.is_insertion());
    }

    #[test]
    fn test_print_alignment_with_gaps_and_empty() {
        let graph = CharGraph::linear(&['A', 'B']);
        let seq = b"AC";
        let aln = vec![
            AlignedPair::new(Some(NodeIndex::<u32>::new(0)), Some(0)),
            AlignedPair::new(Some(NodeIndex::<u32>::new(1)), None),
            AlignedPair::new(None, Some(1)),
        ];

        let printed = print_alignment(&graph, seq, &aln);
        assert_eq!(printed, "AB-\n|  \nA-C");

        let empty_aln: Alignment<NodeIndex> = Vec::new();
        let empty = print_alignment(&graph, seq, &empty_aln);
        assert_eq!(empty, "\n\n");
    }
}
