pub mod salmon_types;
mod util;

use std::collections::HashMap;
use std::fs::*;
use std::io::Write;
use std::io::{self, BufRead, BufReader};
use std::process;

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use ndarray::prelude::*;
use num_format::{Locale, ToFormattedString};
use petgraph as pg;
use petgraph::algo::{connected_components, tarjan_scc};
use petgraph::unionfind::UnionFind;
use rayon::prelude::*;
use std::path::PathBuf;

// Name of the program, to be used in diagnostic messages.
static PROGRAM_NAME: &str = "terminus";

/// Exit the program, printing an error message on stderr, and returning
/// a specific error code. The program name is prefixed onto the front of
/// the error message. If logging is enabled we also write the error
/// message to the log file.
#[allow(unused)]
fn exit_with_error(status: i32, message: &str) {
    //error!(target: "log_messages", "{}", message);
    writeln!(
        &mut std::io::stderr(),
        "{} ERROR: {}!",
        PROGRAM_NAME,
        message
    )
    .unwrap();
    std::process::exit(status);
}

fn do_group(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    //let mut groups: Vec<Vec<usize>> = Vec::new();
    //let dir_paths : Vec<_> = sub_m.values_of("dirs").unwrap().collect();

    let dname: String = sub_m.value_of("dir").unwrap().to_string();
    let prefix: String = sub_m.value_of("out").unwrap().to_string();
    create_dir_all(prefix.clone())?;

    // let mut unionfind_vec = Vec::<UnionFind<_>>::with_capacity(dir_paths.len()) ;
    // let mut gibbs_array_vec = Vec::<Array2<f64>>::with_capacity(dir_paths.len());
    // let mut meta_info_array = Vec::<salmon_types::MetaInfo>::with_capacity(dir_paths.len());
    // let mut num_global_targrts = 0u32 ;

    let seed = sub_m
        .value_of("seed")
        .unwrap()
        .parse::<u64>()
        .expect("generate random values from the seed");
    let min_spread = sub_m
        .value_of("min-spread")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert min-spread to float value");

    let tolerance = sub_m
        .value_of("tolerance")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert tolerance to float value");

    println!("------input configuration------");
    println!("seed : {}", seed);
    println!("min-spread : {}", min_spread);
    println!("tolerance : {}", tolerance);
    println!("dir : {}", dname);
    let compo: Vec<&str> = dname.rsplit('/').collect();
    //println!("{:?}",compo);
    let experiment_name = compo[0];
    let mut prefix_path = prefix;
    prefix_path.push_str("/");
    prefix_path.push_str(experiment_name);
    let file_list = salmon_types::FileList::new(dname.to_string());

    // create output directory
    println!("output folder: {}", prefix_path);
    println!("------------------------------");
    // create
    create_dir_all(prefix_path.clone())?;
    let file_list_out = salmon_types::FileList::new(prefix_path);

    // Load the gibbs samples
    let x = util::parse_json(&file_list.mi_file).unwrap();

    let mut gibbs_array =
        Array2::<f64>::zeros((x.num_valid_targets as usize, x.num_bootstraps as usize));
    util::read_gibbs_array(&file_list.bootstrap_file, &x, &mut gibbs_array);
    let mut gibbs_mat_mean = gibbs_array.mean_axis(Axis(1)).unwrap();

    println!("parsing eqfile {:?}", file_list.eq_file);
    let eq_class = util::parse_eq(&file_list.eq_file).unwrap();
    println!("length of eqclass {:?}", eq_class.neq);
    let mut eq_class_counts = vec![0 as u32; eq_class.neq];
    let mut i = 0 as usize;
    for eq in eq_class.classes.iter() {
        eq_class_counts[i] = eq.2;
        i += 1;
    }

    let inf_perc = 0.25f64;
    let p = util::get_infrv_percentile(&gibbs_array, inf_perc);
    println!("the {}% of infRV was : {}", inf_perc * 100., p);
    let thr = util::get_threhold(&gibbs_array, p, seed, &file_list_out);
    // thr = thr * 0.75;
    // thr = 0.645;
    println!("threshold: {}", thr);
    //let dpath = Path::new(file_list_out.delta_file.clone());
    let mut dfile =
        File::create(file_list_out.delta_file.clone()).expect("could not create collapse.log");
    let mut unionfind_struct = UnionFind::new(eq_class.ntarget);
    let mut gr = util::eq_experiment_to_graph(
        &eq_class,
        &mut gibbs_array,
        &eq_class_counts,
        tolerance,
        thr,
        p,
        min_spread,
        &mut dfile,
        &mut unionfind_struct,
    );
    util::verify_graph(&eq_class_counts, &mut gr);
    // Go over the graph and keep collapsing
    // edges until we hit a point where there
    // are no more edges to that satisfies the criteria
    // and we collapse

    // connected coponents
    let num_connected_components = connected_components(&gr);
    println!("#Connected components {:?}", num_connected_components);

    let mut num_collapses = 0 as usize;

    //let cpath = Path::new(file_list_out.collapsed_log_file.clone());
    let mut cfile = File::create(file_list_out.collapsed_log_file.clone())
        .expect("could not create collapse.log");

    let gcomp: Vec<petgraph::prelude::NodeIndex> = gr
        .node_indices()
        .map(|x| petgraph::graph::NodeIndex::new(x.index()))
        .collect();
    util::work_on_component(
        &eq_class_counts,
        &mut gibbs_array,
        &mut gibbs_mat_mean,
        &mut unionfind_struct,
        &mut gr,
        &gcomp,
        &mut num_collapses,
        thr,
        p,
        &mut cfile,
    );

    //write down the groups
    let mut groups = HashMap::new();
    //let mut grouped_set = HashSet::new();
    for i in 0..(x.num_valid_targets as usize) {
        let root = unionfind_struct.find(i);
        if root != i {
            groups.entry(root).or_insert_with(Vec::new).push(i);
        }
    }

    // println!("Number of collapsed transcripts from conn components with 2 {}", num_collapses_2.to_formatted_string(&Locale::en));
    println!(
        "Number of collapses {}",
        num_collapses.to_formatted_string(&Locale::en)
    );
    //let _res = util::write_modified_quants(&groups, &grouped_set, &file_list_out, &gibbs_array, &x, &rec, &collapsed_dim);
    let mut gfile = File::create(file_list_out.group_file).expect("could not create groups.txt");
    let _write = util::group_writer(&mut gfile, &groups);

    Ok(true)
}

fn do_collapse(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    let num_threads = sub_m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .expect("could not parse --threads option");
    // set the number of rayon threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap_or_else(|_| panic!("could not set number of rayon threads to {}", num_threads));

    let dir_paths: Vec<_> = sub_m.values_of("dirs").unwrap().collect();
    let prefix: String = sub_m.value_of("out").unwrap().to_string();

    let mut global_graph = pg::Graph::<usize, u32, petgraph::Undirected>::new_undirected();

    // add edges
    for (_, dname) in dir_paths.iter().enumerate() {
        let compo: Vec<&str> = dname.rsplit('/').collect();
        let experiment_name = compo[0];
        println!("experiment name {}", experiment_name);
        let mut prefix_path = prefix.clone();
        prefix_path.push_str("/");
        prefix_path.push_str(experiment_name);
        create_dir_all(prefix_path.clone())?;
        let file_list_out = salmon_types::FileList::new(prefix_path);

        //let groups = util::group_reader(&file_list_out.group_file);
        if global_graph.node_count() == 0 {
            let file_list = salmon_types::FileList::new((*dname).to_string());
            let x = util::parse_json(&file_list.mi_file).unwrap();

            for i in 0..x.num_valid_targets {
                let idx = global_graph.add_node(i as usize);
                // the index assigned by the graph should be the
                // order in which we add these
                debug_assert_eq!(i as usize, idx.index());
            }
        }

        let group_file = File::open(file_list_out.group_file).unwrap();
        let buf_reader = BufReader::new(group_file);
        for (_i, l) in buf_reader.lines().enumerate() {
            let s = l.unwrap();
            let v: Vec<_> = s.trim().rsplit(',').collect();
            let group: Vec<usize> = v.iter().map(|n| n.parse::<usize>().unwrap()).collect();

            assert_ne!(group.len(), 1 as usize);
            for g in group.iter() {
                if *g >= global_graph.node_count() {
                    process::exit(0x0100);
                }
            }

            for a in 0..group.len() {
                let na = group[a];
                for nb in group.iter().skip(a + 1) {
                    let va = pg::graph::NodeIndex::new(na);
                    let vb = pg::graph::NodeIndex::new(*nb);
                    let e = global_graph.find_edge(va, vb);
                    match e {
                        Some(ei) => {
                            let ew = global_graph.edge_weight_mut(ei).unwrap();
                            *ew += 1;
                        }
                        None => {
                            global_graph.add_edge(va, vb, 1);
                        }
                    }
                }
            } // all combinations
        } // all lines of the file
    } // all files

    // filter
    let consensus_thresh = sub_m
        .value_of("consensus-thresh")
        .unwrap()
        .parse::<f64>()
        .expect("could not parse --consensus-thresh option");
    let half_length = (dir_paths.len() as f64 * consensus_thresh).floor() as usize;

    let global_filtered_graph = global_graph.filter_map(
        |_, n| Some(*n),
        |_, &e| {
            if (e as usize) >= half_length {
                Some(e)
            } else {
                None
            }
        },
    );

    // components
    let mut comps: Vec<Vec<_>> = tarjan_scc(&global_filtered_graph);
    comps.sort_by(|v, w| v.len().cmp(&w.len()));
    // write a json file containing
    // new number of nodes, number of components, etc
    let (total_in_group, num_group) = comps.iter().fold((0, 0), |s, v| {
        if v.len() > 1 {
            (s.0 + v.len(), s.1 + 1)
        } else {
            (s.0, s.1)
        }
    });
    let term_info = salmon_types::TerminusInfo {
        num_nontrivial_groups: num_group,
        num_targets_in_nontrivial_group: total_in_group,
    };
    println!(
        "total number of transcripts in a non-trivial group : {}",
        total_in_group
    );
    println!("total consensus clusters : {}", num_group);

    // for each experiment call the writer
    println!("=============Writing collapsed output=============");
    dir_paths.into_par_iter().for_each(|dname| {
        // getting prepared for output
        let compo: Vec<&str> = dname.rsplit('/').collect();
        let experiment_name = compo[0];
        let mut prefix_path = prefix.clone();
        prefix_path.push_str("/");
        prefix_path.push_str(experiment_name);
        let file_list_out = salmon_types::FileList::new(prefix_path.to_string());
        let file_list = salmon_types::FileList::new(dname.to_string());
        // create output directory
        println!("Sample {}", experiment_name);

        let bpath = std::path::Path::new(&prefix_path);
        let term_info_path = bpath.join("terminus_info.json");
        let estr = format!("could not write terminus_info.json in {}", prefix_path);
        ::serde_json::to_writer(&File::create(term_info_path).expect(&estr), &term_info)
            .expect(&estr);
        // load old files
        // let x = util::parse_json(&file_list.mi_file).unwrap();
        // copy cmd_info.json from old location to new
        std::fs::copy(&file_list.cmd_file, &file_list_out.cmd_file)
            .expect("Could not copy cmd_info.json.");
        let x = util::parse_json(&file_list.mi_file).unwrap();
        let rec = util::parse_quant(&file_list.quant_file, &x).unwrap();
        let mut gibbs_array =
            Array2::<f64>::zeros((x.num_valid_targets as usize, x.num_bootstraps as usize));
        util::read_gibbs_array(&file_list.bootstrap_file, &x, &mut gibbs_array);

        //call the writer
        let _res =
            util::write_quants_from_components(&comps, &file_list_out, &gibbs_array, &x, &rec);
    });

    Ok(true)
}

fn alevin_processing(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    println!("This is about alevin");
    let dname: String = sub_m.value_of("dir").unwrap().to_string();
    let out_dname: String = sub_m.value_of("out").unwrap().to_string();
    let t2g_file = PathBuf::from(sub_m.value_of("transcript2gene").unwrap().to_string());
    let threshold = sub_m.value_of("threshold-fraction")
                         .unwrap()
                         .parse::<f32>()
                         .expect("could not convert threshold-fraction to float value");

    // read alevin output
    let mut alevin_exp = salmon_types::AlevinMetaData::new(dname);
    alevin_exp.load();
    // [debug] check if alevin is loaded correctly
    println!(
        "Number of cells {:?}, number of features: {:?}",
        alevin_exp.num_of_cells, alevin_exp.num_of_features
    );
    // read alevin matrix
    // extract the tier 3 fraction
    let mut tier_fraction_vec: Vec<f32> = vec![0.0 as f32; alevin_exp.num_of_features];
    // how many cells the gene has been expressed
    let mut gene_expression_count: Vec<f32> = vec![0.0 as f32; alevin_exp.num_of_features];
    util::matrix_reader(
        alevin_exp.tier_file.to_str().unwrap(),
        alevin_exp.quant_file.to_str().unwrap(),
        alevin_exp.num_of_cells,
        alevin_exp.num_of_features,
        &mut tier_fraction_vec,
        &mut gene_expression_count,
    )?;

    // read bfh file
    let bfh_classes =
        util::parse_bfh(&alevin_exp, &t2g_file).expect("Reading bfh class failed");

    // construct gene level graph
    let gr = util::bfh_to_graph(&bfh_classes, &tier_fraction_vec, &alevin_exp, 
        threshold);

    let mut num_connected_components = 0;

    let mut comps: Vec<Vec<_>> = tarjan_scc(&gr);
    comps.sort_by(|v, w| v.len().cmp(&w.len()));

    create_dir_all(out_dname.clone())?;
    let out_prefix = PathBuf::from(out_dname);
    let group_file = out_prefix.join("group.txt");
    let mut gfile = File::create(group_file)?;
    for (_i, comp) in comps.iter().enumerate() {
        if comp.len() != 1 {
            num_connected_components += 1;
            let mut member_names: Vec<String> = Vec::with_capacity(comp.len());
            for g in comp.iter() {
                let node_ind = g.index();
                member_names.push(alevin_exp.feature_vector[node_ind].clone());
            }
            writeln!(gfile, "{}", member_names.join(","))?;
        }
    }
    println!("#Connected components {:?}", num_connected_components);

    Ok(true)
}

// The entry point of the program
// that reads from the files and build
// graphs or later produces the collapsed
// files.

fn main() -> io::Result<()> {
    let matches = App::new("Terminus")
	.setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1.0")
        .author("Sarkar et al.")
        .about("Data-driven grouping of transcripts to reduce inferential uncertainty")
        .subcommand(
            SubCommand::with_name("group")
            .about("perform per-sample grouping of transcripts; required prior to consensus collapse.")
            .arg(
                Arg::with_name("dir")
                    .long("dir")
                    .short("d")
                    .required(true)
                    .takes_value(true)
                    .help("directory to read input from")
            )
            .arg(
                Arg::with_name("min-spread")
                    .long("min-spread")
                    .short("m")
                    .takes_value(true)
                    .default_value("0.1")
                    .help("the minimum spread a transcript must exhibit to enable an \
                          attached edge to be a collapse candidate")
            )
            .arg(
                Arg::with_name("tolerance")
                    .long("tolerance")
                    .takes_value(true)
                    .default_value("0.001")
                    .help("The allowable difference between the weights of transcripts \
                          in same equivalence classes to treat them as identical")
            )
            .arg(
                Arg::with_name("seed")
                    .long("seed")
                    .takes_value(true)
                    .default_value("10")
                    .help("The allowable difference between the weights of transcripts \
                          in same equivalence classes to treat them as identical")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dir")
                    .help("prefix where output would be written")
            )
        )
        .subcommand(
            SubCommand::with_name("collapse")
            .about("analyze a collection of per-sample groups, and produce a consensus grouping.")
            .arg(
                Arg::with_name("dirs")
                    .long("dirs")
                    .short("d")
                    .required(true)
                    .multiple(true)
                    .takes_value(true)
                    .help("direcotories to read the group files from")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dirs")
                    .help("prefix where output would be written")
            )
            .arg(
                Arg::with_name("threads")
                    .long("threads")
                    .short("t")
                    .takes_value(true)
                    .default_value("8")
                    .help("number of threads to use when writing down the collapsed quantifications")
            )
            .arg(
                Arg::with_name("consensus-thresh")
                .long("consensus-thresh")
                .short("c")
                .takes_value(true)
                .default_value("0.5")
                .help("threshold for edge consensus")
            )
        )
        .subcommand(
            SubCommand::with_name("alevin")
            .about("Receive the alevin generated matrix and bootstraps / clusters")
            .arg(
                Arg::with_name("dir")
                    .long("dir")
                    .short("d")
                    .required(true)
                    .takes_value(true)
                    .help("directory to read input from")
            )
            .arg(
                Arg::with_name("transcript2gene")
                    .long("t2g")
                    .short("t")
                    .required(true)
                    .takes_value(true)
                    .help("directory to read input from")
            )
            .arg(
                Arg::with_name("threshold-fraction")
                    .long("threshold-fraction")
                    .short("h")
                    .takes_value(true)
                    .default_value("0.5")
                    .help("Ratio of number of cells where the gene \
                          has to be in tier 3 in order to be considered")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dir")
                    .help("prefix where output would be written")
            )
        ).get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        ("group", Some(sub_m)) => {
            do_group(&sub_m).expect("Grouping failed");
        }
        ("collapse", Some(sub_m)) => {
            do_collapse(&sub_m).expect("Grouping failed");
        }
        ("alevin", Some(sub_m)) => {
            alevin_processing(&sub_m).expect("Alevin ");
        }
        _ => unreachable!(),
    }

    Ok(())
}
