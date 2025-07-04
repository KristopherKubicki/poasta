use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;

#[test]
fn poasta_align_and_view_produces_outputs() {
    let temp = assert_fs::TempDir::new().unwrap();
    let poasta_file = temp.child("out.poasta");
    let gfa_file = temp.child("out.gfa");

    Command::new("cargo")
        .args([
            "run",
            "--bin",
            "poasta",
            "--",
            "align",
            "-o",
            poasta_file.path().to_str().unwrap(),
            "tests/small_test.fa",
        ])
        .assert()
        .success();

    poasta_file.assert(predicate::path::is_file());
    assert!(fs::metadata(poasta_file.path()).unwrap().len() > 0);

    Command::new("cargo")
        .args([
            "run",
            "--bin",
            "poasta",
            "--",
            "view",
            "-O",
            "gfa",
            poasta_file.path().to_str().unwrap(),
            "-o",
            gfa_file.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    gfa_file.assert(predicate::path::is_file());
    assert!(fs::metadata(gfa_file.path()).unwrap().len() > 0);
}
