use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;

#[test]
fn lasagna_align_produces_gaf() {
    let temp = assert_fs::TempDir::new().unwrap();
    let gaf_file = temp.child("out.gaf");

    Command::new("cargo")
        .args([
            "run",
            "--bin",
            "lasagna",
            "--",
            "align",
            "tests/test.gfa",
            "tests/small_test.query.fa",
            "-o",
            gaf_file.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    gaf_file.assert(predicate::path::is_file());
    assert!(fs::metadata(gaf_file.path()).unwrap().len() > 0);
}
