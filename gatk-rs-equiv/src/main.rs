//! gatk-rs-equiv — scientific HaplotypeCaller equivalence vs GATK4 (hap.py / RTG).

use clap::Parser;
use gatk_rs_equiv::cli::{Cli, Command};
use gatk_rs_equiv::differential_fuzz;
use gatk_rs_equiv::report;
use gatk_rs_equiv::run;
use std::process::ExitCode;

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Run(args) => run::run(args),
        Command::Report(args) => report::report(args),
        Command::DifferentialFuzz(args) => differential_fuzz::run_from_cli(args),
    };
    match result {
        Ok(code) => ExitCode::from(code as u8),
        Err(err) => {
            eprintln!("[gatk-rs-equiv] ERROR: {err:#}");
            ExitCode::from(2)
        }
    }
}
