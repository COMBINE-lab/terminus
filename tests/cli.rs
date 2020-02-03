extern crate assert_cmd;

use std::process::Command;

use assert_cmd::prelude::*;
//use predicates::prelude::*;

#[test]
fn terminus_group() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("terminus")?;
    cmd.arg("group");
    cmd.args(&["-d", "tests/data/salmon_quant/quant_1"]);
    cmd.args(&["-m", "0.05"]);
    cmd.args(&["--tolerance", "0.01"]);
    cmd.args(&["-o", "tests/terminus_out"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn terminus_collapse() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("terminus")?;
    cmd.arg("collapse");
    cmd.args(&["-d", "tests/data/salmon_quant/quant_1"]);
    cmd.args(&["-c", "0.25"]);
    cmd.args(&["-o", "tests/data/terminus_group"]);
    
    cmd.assert().success();

    Ok(())
}

