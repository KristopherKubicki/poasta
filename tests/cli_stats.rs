use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn poasta_stats_gfa() {
    Command::cargo_bin("poasta")
        .unwrap()
        .args(["stats", "tests/test.gfa"])
        .assert()
        .success()
        .stderr(predicate::str::contains("node_count:"))
        .stderr(predicate::str::contains("edge_count:"));
}
