use std::cell::{RefCell, RefMut};
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::astar::AstarVisited;
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::Score;
use crate::graphs::{AlignableRefGraph, NodeIndexType};

/// A node in the graph aligned to a certain offset in the query
/// that we will try to extend from
struct StackNode<N, O, I>(AlignmentGraphNode<N, O>, RefCell<I>)
where
    N: NodeIndexType,
    O: OffsetType;

impl<N, O, I> StackNode<N, O, I>
where
    N: NodeIndexType,
    O: OffsetType,
    I: Iterator,
{
    fn new(state: AlignmentGraphNode<N, O>, iter: I)  -> Self {
        Self(state, RefCell::new(iter))
    }

    #[inline(always)]
    fn offset(&self) -> O {
        self.0.offset()
    }

    #[inline(always)]
    fn aln_graph_node(&self) -> &AlignmentGraphNode<N, O> {
        &self.0
    }

    #[inline]
    fn children_iter(&self) -> RefMut<'_, I> {
        self.1.borrow_mut()
    }
}

#[derive(Debug)]
pub enum ExtendResult<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    /// Variant indicating that we reached the end node in the reference graph
    ///
    /// The associated values represent the edge traversed in the alignment graph
    /// and thus represents the parent alignment graph node and the
    /// just reached alignment graph node with N being the reference graph end node.
    RefGraphEnd(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),

    /// Variant indicating that we could not extend any further because all query sequence
    /// has been aligned.
    ///
    /// The associated value is the parent alignment graph node and the
    /// child reference graph node that we tried to align, but couldn't.
    QueryEnd(AlignmentGraphNode<N, O>, N),


    /// During extension, we encountered a reference graph node
    /// with a different symbol than the corresponding query symbol.
    ///
    /// The associated values indicate the edge in the alignment graph traversed, with the
    /// first ['AlignmentGraphNode'] the last node with matching symbols between the
    /// reference and the query, and the second ['AlignmentGraphNode'] the mismatching
    /// alignment graph node.
    ///
    /// We return the edge such that we can open indel states from the last matching
    /// alignment graph node (the first value), in addition to the mismatch state
    /// (the second value).
    Mismatch(AlignmentGraphNode<N, O>, AlignmentGraphNode<N, O>),
}

/// This struct represents the state for POASTA's
/// depth-first greedy alignment algorithm.
///
/// It traverses the reference graph in a depth-first
/// manner, and greedily aligns visited nodes with
/// the corresponding query offset if the characters
/// match.
///
/// This works because we set the alignment cost of
/// matching characters to 0.
pub struct DepthFirstGreedyAlignment<'a, G, O>
    where G: AlignableRefGraph + 'a,
          O: OffsetType,
{
    /// The alignment graph
    ref_graph: &'a G,

    /// Query sequence to align
    seq: &'a [u8],

    /// The current alignment score
    score: Score,

    /// Number of alignment states visited
    num_visited: usize,

    /// Number of alignment states pruned
    num_pruned: usize,

    /// Stack for depth-first alignment of matches between query and graph
    stack: Vec<StackNode<G::NodeIndex, O, G::SuccessorIterator<'a>>>,
}

impl<'a, G, O> DepthFirstGreedyAlignment<'a, G, O>
where
    G: AlignableRefGraph,
    O: OffsetType,
{
    pub fn new(
        ref_graph: &'a G,
        seq: &'a [u8],
        score: Score,
        start_node: &AlignmentGraphNode<G::NodeIndex, O>,
    ) -> Self {
        Self {
            ref_graph,
            seq,
            score,
            num_visited: 0,
            num_pruned: 0,
            stack: vec![StackNode::new(*start_node, ref_graph.successors(start_node.node()))],
        }
    }

    pub fn get_num_visited(&self) -> usize {
        self.num_visited
    }

    pub fn extend<V>(
        &mut self,
        astar_visited: &mut V,
    ) -> Option<ExtendResult<G::NodeIndex, O>>
    where
        V: AstarVisited<G::NodeIndex, O>,
    {
        // Special handling for ends-free alignment starting at offset 0
        if self.stack.len() == 1 && !self.seq.is_empty() {
            let initial_node = *self.stack[0].aln_graph_node();
            if initial_node.offset().as_usize() == 0 {
                if self.ref_graph.is_symbol_equal(initial_node.node(), self.seq[0]) {
                    let match_node = AlignmentGraphNode::new(initial_node.node(), initial_node.offset().increase_one());
                    if astar_visited.update_score_if_lower(&match_node, AlignState::Match, &initial_node, AlignState::Match, self.score) {
                        let child_succ = self.ref_graph.successors(match_node.node());
                        self.stack[0] = StackNode::new(match_node, child_succ);
                        astar_visited.dfa_match(self.score, &initial_node, &match_node);
                        self.num_visited += 1;
                        
                        
                        // If we've consumed the entire query, this might be an end state
                        // Create a special "end" result to signal completion
                        if match_node.offset().as_usize() == self.seq.len() {
                            // Create a virtual end node transition
                            let end_node = AlignmentGraphNode::new(match_node.node(), match_node.offset());
                            return Some(ExtendResult::RefGraphEnd(initial_node, end_node));
                        }
                    }
                }
            }
        }
        
        while !self.stack.is_empty() {
            let parent = self.stack.last().unwrap();

            match self.next_valid_successor(parent, astar_visited) {
                Successor::RefGraphEnd(aln_node) => return Some(ExtendResult::RefGraphEnd(
                    *parent.aln_graph_node(),
                    aln_node
                )),
                Successor::QueryEnd(ref_node) => return Some(ExtendResult::QueryEnd(
                    *parent.aln_graph_node(),
                    ref_node
                )),
                Successor::Match(child) => {
                    if astar_visited.prune(self.score, &child, AlignState::Match) {
                        // eprintln!("- dfa prune");
                        self.num_pruned += 1;
                        continue;
                    }

                    astar_visited.dfa_match(self.score, parent.aln_graph_node(), &child);
                    self.num_visited += 1;

                    let child_succ = self.ref_graph.successors(child.node());
                    self.stack.push(StackNode::new(child, child_succ));
                },
                Successor::Mismatch(child) => {
                    return Some(ExtendResult::Mismatch(
                        *parent.aln_graph_node(),
                        child
                    ))
                }
                Successor::SuccessorsExhausted => {
                    self.stack.pop();
                }
            }
        }

        None
    }

    fn next_valid_successor<V>(
        &self,
        parent: &StackNode<G::NodeIndex, O, G::SuccessorIterator<'a>>,
        // Separate mutable borrow for visited data
        // to prevent requiring mutable &self
        astar_visited: &mut V,
    ) -> Successor<G::NodeIndex, O>
    where
        V: AstarVisited<G::NodeIndex, O>,
    {
        while let Some(child) = parent.children_iter().next() {

            if child == self.ref_graph.end_node() {
                let aln_termination = AlignmentGraphNode::new(child, parent.offset());
                astar_visited.update_score_if_lower(&aln_termination, AlignState::Match,
                                                    parent.aln_graph_node(), AlignState::Match, self.score);

                return Successor::RefGraphEnd(aln_termination);
            }

            if parent.offset().as_usize() >= self.seq.len() {
                return Successor::QueryEnd(child)
            }

            let child_offset = parent.offset().increase_one();
            let child_node = AlignmentGraphNode::new(child, child_offset);

            if self.ref_graph.is_symbol_equal(child, self.seq[child_offset.as_usize()-1]) {
                if astar_visited
                    .update_score_if_lower(&child_node, AlignState::Match,
                                           parent.aln_graph_node(), AlignState::Match, self.score)
                {
                    return Successor::Match(child_node)
                }
            } else {
                return Successor::Mismatch(child_node)
            }
        }

        Successor::SuccessorsExhausted
    }
}

enum Successor<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    RefGraphEnd(AlignmentGraphNode<N, O>),
    QueryEnd(N),
    Match(AlignmentGraphNode<N, O>),
    Mismatch(AlignmentGraphNode<N, O>),
    SuccessorsExhausted,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aligner::Alignment;
    use crate::errors::PoastaError;
    use nonmax::NonMaxU32;
    use std::fmt::Write;
    use petgraph::algo::toposort;
    use petgraph::graph::{DiGraph, NodeIndex, NodeIndices, Neighbors};
    use petgraph::{Incoming, Outgoing};

    /// Minimal graph storing character per node for DFA tests.
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

        fn all_nodes(&self) -> Self::NodeIterator<'_> { self.graph.node_indices() }
        fn node_count(&self) -> usize { self.graph.node_count() }
        fn node_count_with_start_and_end(&self) -> usize { self.graph.node_count() }
        fn edge_count(&self) -> usize { self.graph.edge_count() }
        fn start_node(&self) -> Self::NodeIndex { NodeIndex::new(0) }
        fn end_node(&self) -> Self::NodeIndex { self.end }
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
        fn is_end(&self, node: Self::NodeIndex) -> bool { node == self.end }
        fn get_symbol_char(&self, node: Self::NodeIndex) -> char { self.graph[node] }
        fn is_symbol_equal(&self, node: Self::NodeIndex, s: u8) -> bool { self.graph[node] as u8 == s }
        fn get_node_ordering(&self) -> Vec<usize> {
            let toposorted = toposort(&self.graph, None).unwrap();
            let mut node_ordering = vec![0; toposorted.len()];
            for (rank, node) in toposorted.iter().enumerate() {
                node_ordering[node.index()] = rank;
            }
            node_ordering
        }
    }

    #[derive(Default)]
    struct DummyVisited;

    impl<N: NodeIndexType, O: OffsetType> AstarVisited<N, O> for DummyVisited {
        fn get_score(&self, _aln_node: &AlignmentGraphNode<N, O>, _state: AlignState) -> Score { Score::Unvisited }
        fn set_score(&mut self, _aln_node: &AlignmentGraphNode<N, O>, _state: AlignState, _score: Score) {}
        fn mark_reached(&mut self, _score: Score, _aln_node: &AlignmentGraphNode<N, O>, _state: AlignState) {}
        fn dfa_match(&mut self, _score: Score, _parent: &AlignmentGraphNode<N, O>, _child: &AlignmentGraphNode<N, O>) {}
        fn prune(&self, _score: Score, _aln_node: &AlignmentGraphNode<N, O>, _state: AlignState) -> bool { false }
        fn update_score_if_lower(&mut self, _aln_node: &AlignmentGraphNode<N, O>, _state: AlignState, _parent: &AlignmentGraphNode<N, O>, _parent_state: AlignState, _score: Score) -> bool { true }
        fn backtrace<G>(&self, _ref_graph: &G, _seq: &[u8], _aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
        where G: AlignableRefGraph<NodeIndex = N> { Vec::new() }
        fn write_tsv<W>(&self, _writer: &mut W) -> Result<(), PoastaError>
        where
            W: Write,
        {
            Ok(())
        }
    }

    #[test]
    fn test_extend_refgraph_end() {
        let graph = CharGraph::linear(&['A', 'B', 'C']);
        let seq = b"ABC";
        let start = AlignmentGraphNode::new(graph.start_node(), 0u8);
        let mut dfa = DepthFirstGreedyAlignment::new(&graph, seq, Score::Score(NonMaxU32::new(0).unwrap()), &start);
        let mut visited = DummyVisited::default();
        let res = dfa.extend(&mut visited);
        assert!(matches!(res, Some(ExtendResult::RefGraphEnd(_, _))));
        assert_eq!(dfa.get_num_visited(), 3);
    }

    #[test]
    fn test_extend_query_end() {
        let graph = CharGraph::linear(&['A', 'B', 'C']);
        let seq = b"AB";
        let start = AlignmentGraphNode::new(graph.start_node(), 0u8);
        let mut dfa = DepthFirstGreedyAlignment::new(&graph, seq, Score::Score(NonMaxU32::new(0).unwrap()), &start);
        let mut visited = DummyVisited::default();
        let res = dfa.extend(&mut visited);
        assert!(matches!(res, Some(ExtendResult::QueryEnd(_, _))));
    }

    #[test]
    fn test_extend_mismatch() {
        let graph = CharGraph::linear(&['A', 'B', 'C']);
        let seq = b"ABD";
        let start = AlignmentGraphNode::new(graph.start_node(), 0u8);
        let mut dfa = DepthFirstGreedyAlignment::new(&graph, seq, Score::Score(NonMaxU32::new(0).unwrap()), &start);
        let mut visited = DummyVisited::default();
        let res = dfa.extend(&mut visited);
        assert!(matches!(res, Some(ExtendResult::Mismatch(_, _))));
    }
}

