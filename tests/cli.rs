use std::process::Command;

use assert_cmd::prelude::*;
//use predicates::prelude::*;

#[test]
fn terminus_collapse() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("terminus")?;
    cmd.arg("collapse");
    cmd.args(&["--dirs", "tests/data/???"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn terminus_group() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("terminus")?;
    cmd.arg("group");
    cmd.args(&["--dir", "tests/data/???"]);

    cmd.assert().success();

    Ok(())
}
